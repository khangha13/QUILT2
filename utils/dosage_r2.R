#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  opts <- list(
    imputed_ds = NULL,
    truth_gt = NULL,
    imputed_gp = NULL,
    samples = NULL,
    out_prefix = NULL,
    plots = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--imputed-ds", "--truth-gt", "--imputed-gp", "--samples", "--out-prefix")) {
      if (i == length(args)) stop(key, " requires a value")
      val <- args[[i + 1]]
      switch(key,
             "--imputed-ds" = opts$imputed_ds <- val,
             "--truth-gt" = opts$truth_gt <- val,
             "--imputed-gp" = opts$imputed_gp <- val,
             "--samples" = opts$samples <- val,
             "--out-prefix" = opts$out_prefix <- val)
      i <- i + 2
    } else if (key == "--plots") {
      opts$plots <- TRUE
      i <- i + 1
    } else if (key %in% c("--help", "-h")) {
      cat("Usage: dosage_r2.R --imputed-ds <tsv> --truth-gt <tsv> --samples <file> --out-prefix <prefix> [--imputed-gp <tsv>] [--plots]\n")
      quit(status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
  }
  required <- c("imputed_ds", "truth_gt", "samples", "out_prefix")
  missing <- required[sapply(required, function(k) is.null(opts[[k]]))]
  if (length(missing) > 0) stop("Missing required args: ", paste(missing, collapse = ", "))
  opts
}

opts <- parse_args(args)

sample_ids <- readLines(opts$samples)
if (length(sample_ids) == 0) stop("Sample list is empty")

load_ds_table <- function(path) {
  dt <- fread(path, na.strings = c(".", "NA"), showProgress = FALSE)
  if (ncol(dt) < 6) stop("Expected metadata + at least one sample column in ", path)
  meta <- dt[, c("CHROM", "POS", "REF", "ALT", "ID"), with = FALSE]
  mat <- as.matrix(dt[, -(1:5)])
  storage.mode(mat) <- "double"
  list(meta = meta, mat = mat, samples = colnames(dt)[-(1:5)])
}

load_gt_table <- function(path) {
  dt <- fread(path, na.strings = c(".", "NA"), showProgress = FALSE)
  if (ncol(dt) < 6) stop("Expected metadata + at least one sample column in ", path)
  meta <- dt[, c("CHROM", "POS", "REF", "ALT", "ID"), with = FALSE]
  mat <- as.matrix(dt[, -(1:5)])
  list(meta = meta, mat = mat, samples = colnames(dt)[-(1:5)])
}

gp_to_ds <- function(gp_mat) {
  flat <- as.vector(gp_mat)
  split <- tstrsplit(flat, ",", fixed = TRUE)
  if (length(split) < 3) {
    return(matrix(NA_real_, nrow = nrow(gp_mat), ncol = ncol(gp_mat)))
  }
  p1 <- suppressWarnings(as.numeric(split[[2]]))
  p2 <- suppressWarnings(as.numeric(split[[3]]))
  ds <- p1 + 2 * p2
  dim(ds) <- dim(gp_mat)
  ds
}

gt_to_dosage <- function(gt_mat) {
  flat <- gsub("\\|", "/", gt_mat)
  parts <- tstrsplit(flat, "/", fixed = TRUE)
  a1 <- suppressWarnings(as.numeric(parts[[1]]))
  a2 <- suppressWarnings(as.numeric(parts[[2]]))
  dosage <- a1 + a2
  dosage[is.na(a1) | is.na(a2)] <- NA_real_
  dim(dosage) <- dim(gt_mat)
  dosage
}

calc_variant_metrics <- function(ds_row, truth_row) {
  keep <- !(is.na(ds_row) | is.na(truth_row))
  n_nonmiss <- sum(keep)
  maf <- NA_real_
  if (n_nonmiss > 0) {
    af <- mean(truth_row[keep] / 2)
    maf <- min(af, 1 - af)
  }
  if (n_nonmiss < 2) {
    return(list(r2 = NA_real_, concordance = NA_real_, maf = maf, n = n_nonmiss))
  }
  r2 <- suppressWarnings(cor(ds_row[keep], truth_row[keep])^2)
  ghat <- ifelse(ds_row[keep] >= 1.5, 2,
                 ifelse(ds_row[keep] <= 0.5, 0, 1))
  conc <- mean(ghat == round(truth_row[keep]))
  list(r2 = r2, concordance = conc, maf = maf, n = n_nonmiss)
}

imputed_ds <- load_ds_table(opts$imputed_ds)
truth_gt <- load_gt_table(opts$truth_gt)

if (!identical(imputed_ds$meta, truth_gt$meta)) {
  stop("Imputed and truth tables have different variant ordering or metadata")
}
if (!identical(imputed_ds$samples, truth_gt$samples)) {
  stop("Imputed and truth tables have different sample ordering")
}
if (!identical(imputed_ds$samples, sample_ids)) {
  stop("Extracted sample order does not match provided sample list")
}

ds_mat <- imputed_ds$mat
if (all(is.na(ds_mat)) && !is.null(opts$imputed_gp)) {
  gp_tbl <- load_ds_table(opts$imputed_gp)
  if (!identical(gp_tbl$meta, imputed_ds$meta) || !identical(gp_tbl$samples, imputed_ds$samples)) {
    stop("GP table does not match DS table in ordering")
  }
  ds_from_gp <- gp_to_ds(gp_tbl$mat)
  ds_mat <- ds_from_gp
} else if (!is.null(opts$imputed_gp)) {
  gp_tbl <- load_ds_table(opts$imputed_gp)
  if (identical(gp_tbl$meta, imputed_ds$meta) && identical(gp_tbl$samples, imputed_ds$samples)) {
    ds_from_gp <- gp_to_ds(gp_tbl$mat)
    replace_idx <- is.na(ds_mat)
    ds_mat[replace_idx] <- ds_from_gp[replace_idx]
  }
}

if (all(is.na(ds_mat))) {
  stop("No usable dosage values were found (DS and GP missing/empty)")
}

truth_dosage <- gt_to_dosage(truth_gt$mat)

metrics_list <- vector("list", nrow(imputed_ds$meta))
for (i in seq_len(nrow(imputed_ds$meta))) {
  res <- calc_variant_metrics(ds_mat[i, ], truth_dosage[i, ])
  metrics_list[[i]] <- cbind(imputed_ds$meta[i], data.table(
    maf = res$maf,
    n_non_missing = res$n,
    r2 = res$r2,
    concordance = res$concordance
  ))
}
metrics_dt <- rbindlist(metrics_list)

summary_dt <- data.table(
  metric = c("r2_mean", "r2_median", "concordance_mean", "variants_total", "variants_with_r2"),
  value = c(
    mean(metrics_dt$r2, na.rm = TRUE),
    median(metrics_dt$r2, na.rm = TRUE),
    mean(metrics_dt$concordance, na.rm = TRUE),
    nrow(metrics_dt),
    sum(!is.na(metrics_dt$r2))
  )
)

maf_bins <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
metrics_dt[, maf_bin := cut(maf, breaks = maf_bins, include.lowest = TRUE, right = FALSE)]
maf_bins_dt <- metrics_dt[!is.na(maf_bin), .(
  variants = .N,
  r2_mean = mean(r2, na.rm = TRUE),
  concordance_mean = mean(concordance, na.rm = TRUE)
), by = maf_bin]

out_prefix <- opts$out_prefix
fwrite(metrics_dt, sprintf("%s.metrics.tsv", out_prefix), sep = "\t")
fwrite(summary_dt, sprintf("%s.summary.tsv", out_prefix), sep = "\t", col.names = TRUE)
fwrite(maf_bins_dt, sprintf("%s.maf_bins.tsv", out_prefix), sep = "\t")

open_png <- function(path, width = 1600, height = 1200, res = 200) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = path, width = width, height = height, units = "px", res = res)
  } else {
    grDevices::png(filename = path, width = width / res, height = height / res, units = "in", res = res)
  }
}

if (isTRUE(opts$plots)) {
  if (nrow(metrics_dt) > 0) {
    open_png(sprintf("%s.r2_hist.png", out_prefix))
    print(ggplot(metrics_dt[!is.na(r2)], aes(x = r2)) +
            geom_histogram(bins = 60, fill = "#1f77b4", color = "white") +
            theme_minimal(base_size = 12) +
            labs(title = "Imputation r2 distribution", x = expression(r^2), y = "Variants"))
    dev.off()

    open_png(sprintf("%s.concordance_hist.png", out_prefix))
    print(ggplot(metrics_dt[!is.na(concordance)], aes(x = concordance)) +
            geom_histogram(bins = 40, fill = "#2ca02c", color = "white") +
            theme_minimal(base_size = 12) +
            labs(title = "Genotype concordance", x = "Concordance", y = "Variants"))
    dev.off()

    open_png(sprintf("%s.r2_vs_maf.png", out_prefix))
    print(ggplot(metrics_dt[!is.na(r2) & !is.na(maf)], aes(x = maf, y = r2)) +
            geom_hex(bins = 40) +
            scale_fill_viridis_c(option = "plasma") +
            theme_minimal(base_size = 12) +
            labs(title = "r2 vs MAF (truth-derived)", x = "Minor allele frequency", y = expression(r^2)))
    dev.off()
  }
}

cat(sprintf("Metrics written with prefix: %s\n", out_prefix))
