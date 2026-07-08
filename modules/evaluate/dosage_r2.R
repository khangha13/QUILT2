#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  opts <- list(
    imputed_ds = NULL,
    imputed_gt = NULL,
    truth_gt = NULL,
    imputed_gp = NULL,
    samples = NULL,
    eval_dir = NULL,
    per_sample_out = NULL,
    parquet_out = NULL,
    use_vcfpp = FALSE,
    vcfpp_imputed = NULL,
    vcfpp_truth = NULL
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--imputed-ds", "--imputed-gt", "--truth-gt", "--imputed-gp", "--samples",
                   "--eval-dir", "--per-sample-out", "--parquet-out",
                   "--vcfpp-imputed", "--vcfpp-truth")) {
      if (i == length(args)) stop(key, " requires a value")
      val <- args[[i + 1]]
      switch(key,
             "--imputed-ds" = opts$imputed_ds <- val,
             "--imputed-gt" = opts$imputed_gt <- val,
             "--truth-gt" = opts$truth_gt <- val,
             "--imputed-gp" = opts$imputed_gp <- val,
             "--samples" = opts$samples <- val,
             "--eval-dir" = opts$eval_dir <- val,
             "--per-sample-out" = opts$per_sample_out <- val,
             "--parquet-out" = opts$parquet_out <- val,
             "--vcfpp-imputed" = opts$vcfpp_imputed <- val,
             "--vcfpp-truth" = opts$vcfpp_truth <- val)
      i <- i + 2
    } else if (key == "--use-vcfpp") {
      opts$use_vcfpp <- TRUE
      i <- i + 1
    } else if (key %in% c("--help", "-h")) {
      cat(paste(
        "Usage: dosage_r2.R --imputed-gt <tsv> --truth-gt <tsv> --samples <file>",
        "--eval-dir <dir> --per-sample-out <path>",
        "[--parquet-out <path>] [--imputed-ds <tsv>] [--imputed-gp <tsv>]",
        "[--use-vcfpp --vcfpp-imputed <vcf> --vcfpp-truth <vcf>]\n",
        sep = " "
      ))
      quit(status = 0)
    } else {
      stop("Unknown argument: ", key)
    }
  }
  if (is.null(opts$imputed_gt) && is.null(opts$imputed_ds)) {
    stop("Provide --imputed-gt (preferred) or --imputed-ds")
  }
  required <- c("truth_gt", "samples", "eval_dir", "per_sample_out")
  missing <- required[sapply(required, function(k) is.null(opts[[k]]))]
  if (length(missing) > 0) stop("Missing required args: ", paste(missing, collapse = ", "))
  if (isTRUE(opts$use_vcfpp) && (is.null(opts$vcfpp_imputed) || is.null(opts$vcfpp_truth))) {
    stop("When --use-vcfpp is set, provide --vcfpp-imputed and --vcfpp-truth VCFs")
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
  # Position-only key. Genotypes have been translated to A/B format by
  # dosage_r2.sh; REF/ALT columns retain the original allele strings.
  # Matching on (CHROM, POS) is sufficient.
  paste(meta_dt$CHROM, meta_dt$POS, sep = ":")
}

align_to_truth <- function(truth_meta, imputed_meta, imputed_mat, label = "imputed") {
  tk <- variant_key(truth_meta)
  ik <- variant_key(imputed_meta)

  if (anyDuplicated(tk)) stop("Truth table has duplicated variants by (CHROM,POS)")
  if (anyDuplicated(ik)) stop(label, " table has duplicated variants by (CHROM,POS)")

  idx <- match(tk, ik)
  if (any(is.na(idx))) {
    missing_keys <- tk[is.na(idx)]
    stop(label, " table is missing ", length(missing_keys),
         " truth variants by (CHROM,POS). Example missing key: ",
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
  # Convert A/B genotype strings to numeric dosage.
  # Input genotypes are in A/B format produced by dosage_r2.sh:
  #   A/A → 0  (homozygous AT-group / array REF)
  #   A/B → 1  (heterozygous)
  #   B/A → 1  (heterozygous, reversed order)
  #   B/B → 2  (homozygous CG-group / array ALT)
  #   ./. → NA (missing)
  # Phased separators (|) are normalized to / before matching.
  flat <- gsub("\\|", "/", gt_mat)
  dosage <- rep(NA_real_, length(flat))
  dosage[flat == "A/A"] <- 0
  dosage[flat == "A/B" | flat == "B/A"] <- 1
  dosage[flat == "B/B"] <- 2
  # Everything else (including "./." and any unexpected strings) stays NA.
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
    return(list(r = NA_real_, r2 = NA_real_, concordance = NA_real_, maf = maf, n = n_nonmiss))
  }
  # r is the signed Pearson correlation; r2 = r^2 (direction-aware).
  # A negative r flags variants where imputed and truth dosages are anti-correlated —
  # indicating a possible translation or allele-encoding error for that variant.
  r  <- suppressWarnings(cor(ds_row[keep], truth_row[keep]))
  r2 <- if (!is.na(r)) r^2 else NA_real_
  ghat <- ifelse(ds_row[keep] >= 1.5, 2,
                 ifelse(ds_row[keep] <= 0.5, 0, 1))
  conc <- mean(ghat == round(truth_row[keep]))
  list(r = r, r2 = r2, concordance = conc, maf = maf, n = n_nonmiss)
}

write_concordance_parquet <- function(meta_dt, imputed_ds, truth_ds, sample_ids, parquet_out) {
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
  arrow::write_parquet(conc_dt, parquet_out)
}

run_vcfpp <- function(imputed_vcf, truth_vcf, eval_dir) {
  if (!requireNamespace("vcfppR", quietly = TRUE)) {
    warning("vcfppR not installed; skipping vcfpp evaluation")
    return(invisible(NULL))
  }
  qc_dir <- file.path(eval_dir, "intermediate", "qc")
  res <- tryCatch(
    vcfppR::vcfcomp(test = imputed_vcf, truth = truth_vcf, stats = "r2", formats = c("GT")),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    warning("vcfppR::vcfcomp failed: ", conditionMessage(res))
    return(invisible(NULL))
  }
  saveRDS(res, file.path(qc_dir, "vcfpp.rds"))
  dt <- NULL
  if (is.data.frame(res)) {
    dt <- as.data.table(res)
  } else if (is.list(res) && !is.null(res$stats)) {
    dt <- tryCatch(as.data.table(res$stats), error = function(e) NULL)
  }
  if (!is.null(dt)) {
    fwrite(dt, file.path(qc_dir, "vcfpp.r2.tsv"), sep = "\t")
  }
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

# =============================================================================
# PER-VARIANT METRICS
# =============================================================================
# Both imputed and truth genotypes are translated using the same nucleotide-based
# A/B rule (dosage_r2.sh step 5): decode GT index → nucleotide, then assign
# A (AT-group, dosage 0) or B (CG-group, dosage 2). Because the rule is
# identical for both VCFs, a biologically identical genotype always yields the
# same dosage value regardless of which allele is listed as REF or ALT in
# either file. No post-hoc direction flipping is needed or appropriate.
#
# r²          = cor(imputed_dosage, truth_dosage)²   per variant
# concordance = mean(round(imputed) == round(truth))  per variant

metrics_list <- vector("list", nrow(variant_meta))
for (i in seq_len(nrow(variant_meta))) {
  res <- calc_variant_metrics(ds_mat[i, ], truth_dosage[i, ])
  metrics_list[[i]] <- cbind(variant_meta[i], data.table(
    maf           = res$maf,
    n_non_missing = res$n,
    r             = res$r,
    r2            = res$r2,
    concordance   = res$concordance
  ))
}
metrics_dt <- rbindlist(metrics_list)

# --- Per-sample metrics ---
# For each sample, compute overall r² across all variants, plus r² within each
# 0.1 MAF bin. This helps identify individual samples that may be pulling the
# dataset r² down, and whether the issue is concentrated in specific MAF ranges.

calc_sample_metrics <- function(ds_col, truth_col) {
  keep <- !(is.na(ds_col) | is.na(truth_col))
  n_nonmiss <- sum(keep)
  if (n_nonmiss < 2) return(list(r = NA_real_, r2 = NA_real_, n = n_nonmiss))
  r  <- suppressWarnings(cor(ds_col[keep], truth_col[keep]))
  r2 <- if (!is.na(r)) r^2 else NA_real_
  list(r = r, r2 = r2, n = n_nonmiss)
}

# MAF bins for per-sample breakdown: [0.0,0.1), [0.1,0.2), ..., [0.4,0.5]
sample_maf_breaks <- seq(0, 0.5, by = 0.1)
sample_maf_breaks[length(sample_maf_breaks)] <- sample_maf_breaks[length(sample_maf_breaks)] + 1e-8
sample_maf_labels <- c("[0.0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)", "[0.4,0.5]")
# Assign each variant to a 0.1 MAF bin (using truth-derived MAF from metrics_dt).
variant_maf <- metrics_dt$maf
variant_maf_bin <- cut(variant_maf, breaks = sample_maf_breaks,
                       include.lowest = TRUE, right = FALSE, labels = sample_maf_labels)
maf_bin_labels <- levels(variant_maf_bin)

sample_metrics_list <- vector("list", length(sample_order))
for (j in seq_along(sample_order)) {
  # Overall r and r².
  res <- calc_sample_metrics(ds_mat[, j], truth_dosage[, j])
  row <- data.table(sample     = sample_order[j],
                    r_overall  = res$r,
                    r2_overall = res$r2,
                    n_variants = res$n)

  # r and r² per 0.1 MAF bin.
  for (b in maf_bin_labels) {
    idx <- which(variant_maf_bin == b)
    if (length(idx) >= 2) {
      bres <- calc_sample_metrics(ds_mat[idx, j], truth_dosage[idx, j])
      set(row, j = paste0("r_maf_",  b), value = bres$r)
      set(row, j = paste0("r2_maf_", b), value = bres$r2)
      set(row, j = paste0("n_maf_",  b), value = bres$n)
    } else {
      set(row, j = paste0("r_maf_",  b), value = NA_real_)
      set(row, j = paste0("r2_maf_", b), value = NA_real_)
      set(row, j = paste0("n_maf_",  b), value = length(idx))
    }
  }

  sample_metrics_list[[j]] <- row
}
sample_metrics_dt <- rbindlist(sample_metrics_list)
setorder(sample_metrics_dt, r2_overall)  # worst samples first

fwrite(sample_metrics_dt, opts$per_sample_out, sep = "\t")

if (!is.null(opts$parquet_out)) {
  write_concordance_parquet(variant_meta, ds_mat, truth_dosage, sample_order, opts$parquet_out)
}

if (isTRUE(opts$use_vcfpp)) {
  run_vcfpp(opts$vcfpp_imputed, opts$vcfpp_truth, opts$eval_dir)
}

cat(sprintf("Per-sample metrics written to: %s\n", opts$per_sample_out))
