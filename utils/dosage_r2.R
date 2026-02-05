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
    imputed_gt = NULL,
    truth_gt = NULL,
    imputed_gp = NULL,
    samples = NULL,
    out_prefix = NULL,
    plots = FALSE,
    use_vcfpp = FALSE,
    vcfpp_imputed = NULL,
    vcfpp_truth = NULL,
    write_parquet = TRUE
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--imputed-ds", "--imputed-gt", "--truth-gt", "--imputed-gp", "--samples", "--out-prefix", "--vcfpp-imputed", "--vcfpp-truth", "--write-parquet")) {
      if (i == length(args)) stop(key, " requires a value")
      val <- args[[i + 1]]
      switch(key,
             "--imputed-ds" = opts$imputed_ds <- val,
             "--imputed-gt" = opts$imputed_gt <- val,
             "--truth-gt" = opts$truth_gt <- val,
             "--imputed-gp" = opts$imputed_gp <- val,
             "--samples" = opts$samples <- val,
             "--out-prefix" = opts$out_prefix <- val,
             "--vcfpp-imputed" = opts$vcfpp_imputed <- val,
             "--vcfpp-truth" = opts$vcfpp_truth <- val,
             "--write-parquet" = opts$write_parquet <- val)
      i <- i + 2
    } else if (key == "--plots") {
      opts$plots <- TRUE
      i <- i + 1
    } else if (key == "--use-vcfpp") {
      opts$use_vcfpp <- TRUE
      i <- i + 1
    } else if (key %in% c("--help", "-h")) {
      cat("Usage: dosage_r2.R --imputed-gt <tsv> --truth-gt <tsv> --samples <file> --out-prefix <prefix> [--imputed-ds <tsv>] [--imputed-gp <tsv>] [--use-vcfpp --vcfpp-imputed <vcf> --vcfpp-truth <vcf>] [--write-parquet <true|false>] [--plots]\n")
      quit(status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
  }
  if (is.null(opts$imputed_gt) && is.null(opts$imputed_ds)) {
    stop("Provide --imputed-gt (preferred) or --imputed-ds")
  }
  required <- c("truth_gt", "samples", "out_prefix")
  missing <- required[sapply(required, function(k) is.null(opts[[k]]))]
  if (length(missing) > 0) stop("Missing required args: ", paste(missing, collapse = ", "))
  if (isTRUE(opts$use_vcfpp) && (is.null(opts$vcfpp_imputed) || is.null(opts$vcfpp_truth))) {
    stop("When --use-vcfpp is set, provide --vcfpp-imputed and --vcfpp-truth VCFs")
  }
  if (is.character(opts$write_parquet)) {
    opts$write_parquet <- tolower(opts$write_parquet) %in% c("true", "1", "yes", "y")
  }
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

variant_key <- function(meta_dt) {
  # Stable join key independent of VCF ID field.
  paste(meta_dt$CHROM, meta_dt$POS, meta_dt$REF, meta_dt$ALT, sep = ":")
}

align_to_truth <- function(truth_meta, imputed_meta, imputed_mat, label = "imputed") {
  tk <- variant_key(truth_meta)
  ik <- variant_key(imputed_meta)

  if (anyDuplicated(tk)) stop("Truth table has duplicated variants by (CHROM,POS,REF,ALT)")
  if (anyDuplicated(ik)) stop(label, " table has duplicated variants by (CHROM,POS,REF,ALT)")

  idx <- match(tk, ik)
  if (any(is.na(idx))) {
    missing_keys <- tk[is.na(idx)]
    stop(label, " table is missing ", length(missing_keys),
         " truth variants by (CHROM,POS,REF,ALT). Example missing key: ",
         missing_keys[[1]])
  }

  # Reorder to match truth. Use truth metadata (including ID) for downstream outputs.
  list(meta = truth_meta, mat = imputed_mat[idx, , drop = FALSE])
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

open_png <- function(path, width = 1600, height = 1200, res = 200) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = path, width = width, height = height, units = "px", res = res)
  } else {
    grDevices::png(filename = path, width = width / res, height = height / res, units = "in", res = res)
  }
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

write_concordance_parquet <- function(meta_dt, imputed_ds, truth_ds, sample_ids, out_prefix) {
  if (!requireNamespace("arrow", quietly = TRUE)) {
    warning("arrow package not installed; skipping concordance parquet")
    return(invisible(NULL))
  }
  imp_gt <- round(imputed_ds)
  truth_gt <- round(truth_ds)
  conc_mat <- matrix(NA_real_, nrow = nrow(imp_gt), ncol = ncol(imp_gt))
  keep <- !(is.na(imp_gt) | is.na(truth_gt))
  conc_mat[keep] <- as.numeric(imp_gt[keep] == truth_gt[keep])
  colnames(conc_mat) <- sample_ids
  conc_dt <- cbind(as.data.table(meta_dt), as.data.table(conc_mat))
  arrow::write_parquet(conc_dt, sprintf("%s.concordance.parquet", out_prefix))
}

run_vcfpp <- function(imputed_vcf, truth_vcf, out_prefix) {
  if (!requireNamespace("vcfppR", quietly = TRUE)) {
    warning("vcfppR not installed; skipping vcfpp evaluation")
    return(invisible(NULL))
  }
  res <- tryCatch(
    vcfppR::vcfcomp(test = imputed_vcf, truth = truth_vcf, stats = "r2", formats = c("GT")),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    warning("vcfppR::vcfcomp failed: ", conditionMessage(res))
    return(invisible(NULL))
  }
  saveRDS(res, sprintf("%s.vcfpp.rds", out_prefix))
  dt <- NULL
  if (is.data.frame(res)) {
    dt <- as.data.table(res)
  } else if (is.list(res) && !is.null(res$stats)) {
    dt <- tryCatch(as.data.table(res$stats), error = function(e) NULL)
  }
  if (!is.null(dt)) {
    fwrite(dt, sprintf("%s.vcfpp.r2.tsv", out_prefix), sep = "\t")
  }
  pngfile <- sprintf("%s.vcfpp.r2.png", out_prefix)
  open_png(pngfile)
  tryCatch({
    vcfppR::vcfplot(res, col = 2, cex = 1.5, lwd = 2, type = "b")
  }, error = function(e) {
    warning("vcfppR::vcfplot failed: ", conditionMessage(e))
  })
  grDevices::dev.off()
}

truth_gt <- load_gt_table(opts$truth_gt)
if (!identical(truth_gt$samples, sample_ids)) {
  stop("Extracted sample order does not match provided sample list")
}

variant_meta <- truth_gt$meta
sample_order <- truth_gt$samples

if (!is.null(opts$imputed_gt)) {
  imputed_gt <- load_gt_table(opts$imputed_gt)
  aligned <- align_to_truth(variant_meta, imputed_gt$meta, imputed_gt$mat, label = "Imputed GT")
  imputed_gt$meta <- aligned$meta
  imputed_gt$mat <- aligned$mat
  if (!identical(imputed_gt$samples, sample_order)) {
    stop("Imputed GT table sample order differs from truth")
  }
  ds_mat <- gt_to_dosage(imputed_gt$mat)
} else {
  imputed_ds <- load_ds_table(opts$imputed_ds)
  aligned <- align_to_truth(variant_meta, imputed_ds$meta, imputed_ds$mat, label = "Imputed DS")
  imputed_ds$meta <- aligned$meta
  imputed_ds$mat <- aligned$mat
  if (!identical(imputed_ds$samples, sample_order)) {
    stop("Imputed and truth tables have different sample ordering")
  }
  ds_mat <- imputed_ds$mat

  if (all(is.na(ds_mat)) && !is.null(opts$imputed_gp)) {
    gp_tbl <- load_gt_table(opts$imputed_gp)
    if (!identical(gp_tbl$samples, imputed_ds$samples)) stop("GP table sample order does not match DS table")
    gp_aligned <- align_to_truth(variant_meta, gp_tbl$meta, gp_tbl$mat, label = "Imputed GP")
    gp_tbl$meta <- gp_aligned$meta
    gp_tbl$mat <- gp_aligned$mat
    ds_from_gp <- gp_to_ds(gp_tbl$mat)
    ds_mat <- ds_from_gp
  } else if (!is.null(opts$imputed_gp)) {
    gp_tbl <- load_gt_table(opts$imputed_gp)
    if (!identical(gp_tbl$samples, imputed_ds$samples)) stop("GP table sample order does not match DS table")
    gp_aligned <- align_to_truth(variant_meta, gp_tbl$meta, gp_tbl$mat, label = "Imputed GP")
    gp_tbl$meta <- gp_aligned$meta
    gp_tbl$mat <- gp_aligned$mat
    ds_from_gp <- gp_to_ds(gp_tbl$mat)
    replace_idx <- is.na(ds_mat)
    ds_mat[replace_idx] <- ds_from_gp[replace_idx]
  }

  if (all(is.na(ds_mat))) {
    stop("No usable dosage values were found (DS and GP missing/empty)")
  }
}

truth_dosage <- gt_to_dosage(truth_gt$mat)

metrics_list <- vector("list", nrow(variant_meta))
for (i in seq_len(nrow(variant_meta))) {
  res <- calc_variant_metrics(ds_mat[i, ], truth_dosage[i, ])
  metrics_list[[i]] <- cbind(variant_meta[i], data.table(
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

if (isTRUE(opts$write_parquet)) {
  write_concordance_parquet(variant_meta, ds_mat, truth_dosage, sample_order, out_prefix)
}

if (isTRUE(opts$use_vcfpp)) {
  run_vcfpp(opts$vcfpp_imputed, opts$vcfpp_truth, out_prefix)
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
