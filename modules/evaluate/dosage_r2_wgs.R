#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("R package 'data.table' is required for WGS GT evaluation.", call. = FALSE)
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("R package 'arrow' is required for WGS GT evaluation.", call. = FALSE)
  }
})

library(data.table)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Expected --key value arguments; found: ", key, call. = FALSE)
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

require_options <- function(opt, names) {
  missing <- names[!vapply(names, function(x) !is.null(opt[[x]]) && nzchar(opt[[x]]), logical(1))]
  if (length(missing)) stop("Missing options: ", paste(missing, collapse = ", "), call. = FALSE)
}

as_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  x <- as.character(x)
  x[x %chin% c("", ".")] <- NA_character_
  suppressWarnings(as.numeric(x))
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  fwrite(x, tmp, sep = "\t", na = "NA", quote = FALSE)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not promote output: ", path, call. = FALSE)
}

atomic_parquet <- function(x, dataset_dir, chromosome) {
  partition_dir <- file.path(dataset_dir, paste0("CHROM=", chromosome))
  dir.create(partition_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(partition_dir, "part-000.parquet")
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  payload <- copy(x)
  if ("CHROM" %chin% names(payload)) payload[, CHROM := NULL]
  arrow::write_parquet(payload, tmp, compression = "snappy")
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) stop("Could not promote Parquet output: ", path, call. = FALSE)
}

cor_from_sums <- function(n, sx, sy, sxx, syy, sxy, min_pairs) {
  numerator <- n * sxy - sx * sy
  vx <- n * sxx - sx * sx
  vy <- n * syy - sy * sy
  valid <- n >= min_pairs & is.finite(vx) & is.finite(vy) & vx > 0 & vy > 0
  r <- rep(NA_real_, length(n))
  r[valid] <- numerator[valid] / sqrt(vx[valid] * vy[valid])
  r[valid] <- pmax(-1, pmin(1, r[valid]))
  list(r = r, r2 = r * r)
}

maf_labels <- c("[0.0,0.1)", "[0.1,0.2)", "[0.2,0.3)", "[0.3,0.4)", "[0.4,0.5]")
maf_breaks <- c(0, 0.1, 0.2, 0.3, 0.4, 0.50000001)

decode_gt <- function(gt) {
  gt <- as.character(gt)
  missing <- is.na(gt) | gt %chin% c("", ".", "./.", ".|.")
  valid_pattern <- "^(0[/|]0|0[/|]1|1[/|]0|1[/|]1)$"
  invalid <- !missing & !grepl(valid_pattern, gt)
  dosage <- rep(NA_real_, length(gt))
  dosage[gt %chin% c("0/0", "0|0")] <- 0
  dosage[gt %chin% c("0/1", "1/0", "0|1", "1|0")] <- 1
  dosage[gt %chin% c("1/1", "1|1")] <- 2
  list(dosage = dosage, missing = missing, invalid = invalid)
}

read_raw <- function(path) {
  fread(path, sep = "\t", header = TRUE, na.strings = "NA", check.names = FALSE)
}

run_chromosome <- function(opt) {
  require_options(opt, c(
    "imputed-common-raw", "truth-common-raw", "imputed-only-raw", "truth-only-raw",
    "samples", "out-dir", "chromosome", "filter-enabled",
    "min-qual", "min-qd", "max-sor", "max-fs", "min-mq",
    "min-mq-rank-sum", "min-read-pos-rank-sum", "min-gq", "min-dp"
  ))

  threads <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  if (is.na(threads) || threads < 1L) threads <- 1L
  setDTthreads(threads)

  chromosome <- opt[["chromosome"]]
  filter_enabled <- identical(opt[["filter-enabled"]], "true")
  thresholds <- list(
    min_qual = as.numeric(opt[["min-qual"]]), min_qd = as.numeric(opt[["min-qd"]]),
    max_sor = as.numeric(opt[["max-sor"]]), max_fs = as.numeric(opt[["max-fs"]]),
    min_mq = as.numeric(opt[["min-mq"]]),
    min_mq_rank_sum = as.numeric(opt[["min-mq-rank-sum"]]),
    min_read_pos_rank_sum = as.numeric(opt[["min-read-pos-rank-sum"]]),
    min_gq = as.numeric(opt[["min-gq"]]), min_dp = as.numeric(opt[["min-dp"]])
  )

  samples <- readLines(opt[["samples"]], warn = FALSE)
  samples <- samples[nzchar(samples)]
  if (!length(samples) || anyDuplicated(samples)) {
    stop("The sample file must contain at least one unique sample ID.", call. = FALSE)
  }

  imputed <- read_raw(opt[["imputed-common-raw"]])
  truth <- read_raw(opt[["truth-common-raw"]])
  imputed_only <- read_raw(opt[["imputed-only-raw"]])
  truth_only <- read_raw(opt[["truth-only-raw"]])
  variant_cols <- c("CHROM", "POS", "REF", "ALT", "ID")
  key_cols <- c("CHROM", "POS", "REF", "ALT")
  truth_site_cols <- c("QUAL", "QD", "SOR", "FS", "MQ", "MQRankSum", "ReadPosRankSum")

  missing_imputed <- setdiff(c(variant_cols, samples), names(imputed))
  required_truth <- c(variant_cols, truth_site_cols, paste0(samples, ".GT"))
  if (filter_enabled) required_truth <- c(required_truth, paste0(samples, ".GQ"), paste0(samples, ".DP"))
  missing_truth <- setdiff(required_truth, names(truth))
  missing_imputed_only <- setdiff(variant_cols, names(imputed_only))
  missing_truth_only <- setdiff(variant_cols, names(truth_only))
  if (length(missing_imputed)) stop("Missing shared imputed columns: ", paste(missing_imputed, collapse = ", "), call. = FALSE)
  if (length(missing_truth)) stop("Missing shared truth columns: ", paste(missing_truth, collapse = ", "), call. = FALSE)
  if (length(missing_imputed_only)) stop("Missing imputed-only columns: ", paste(missing_imputed_only, collapse = ", "), call. = FALSE)
  if (length(missing_truth_only)) stop("Missing truth-only columns: ", paste(missing_truth_only, collapse = ", "), call. = FALSE)
  if (!nrow(imputed) || !nrow(truth)) stop("No exact shared variants remain after bcftools isec.", call. = FALSE)
  if (uniqueN(imputed, by = key_cols) != nrow(imputed)) stop("Duplicate exact variants in shared imputed input.", call. = FALSE)
  if (uniqueN(truth, by = key_cols) != nrow(truth)) stop("Duplicate exact variants in shared truth input.", call. = FALSE)

  setorderv(imputed, key_cols)
  setorderv(truth, key_cols)
  same_keys <- nrow(imputed) == nrow(truth) &&
    all(vapply(key_cols, function(column) identical(imputed[[column]], truth[[column]]), logical(1)))
  if (!same_keys) {
    stop("bcftools isec shared records are not aligned by CHROM, POS, REF, and ALT.", call. = FALSE)
  }
  setorderv(imputed_only, key_cols)
  setorderv(truth_only, key_cols)

  qual <- as_num(truth$QUAL)
  qd <- as_num(truth$QD)
  sor <- as_num(truth$SOR)
  fs <- as_num(truth$FS)
  mq <- as_num(truth$MQ)
  mq_rank <- as_num(truth$MQRankSum)
  read_pos_rank <- as_num(truth$ReadPosRankSum)
  site_fail <- list(
    missing_QUAL = is.na(qual),
    QUAL_below_min = !is.na(qual) & qual < thresholds$min_qual,
    missing_QD = is.na(qd),
    QD_below_min = !is.na(qd) & qd < thresholds$min_qd,
    SOR_above_max = !is.na(sor) & sor > thresholds$max_sor,
    FS_above_max = !is.na(fs) & fs > thresholds$max_fs,
    missing_MQ = is.na(mq),
    MQ_below_min = !is.na(mq) & mq < thresholds$min_mq,
    MQRankSum_below_min = !is.na(mq_rank) & mq_rank < thresholds$min_mq_rank_sum,
    ReadPosRankSum_below_min = !is.na(read_pos_rank) & read_pos_rank < thresholds$min_read_pos_rank_sum
  )
  if (filter_enabled) {
    site_pass <- !Reduce(`|`, site_fail)
  } else {
    site_pass <- rep(TRUE, nrow(truth))
    site_fail <- lapply(site_fail, function(x) rep(FALSE, length(x)))
  }
  site_reasons <- rep("PASS", nrow(truth))
  for (reason in names(site_fail)) {
    failed <- which(site_fail[[reason]])
    if (length(failed)) {
      site_reasons[failed] <- ifelse(
        site_reasons[failed] == "PASS", reason, paste0(site_reasons[failed], ";", reason)
      )
    }
  }

  matched_variants <- imputed[, .(CHROM, POS, REF, ALT, ID)]
  matched_variants[, `:=`(
    ID_truth = truth$ID,
    site_filter_pass = site_pass,
    site_filter_reasons = site_reasons
  )]
  site_filtered_variants <- matched_variants[site_filter_pass == FALSE]

  imeta <- rbindlist(list(
    imputed[, .(CHROM, POS, REF, ALT, ID_imputed = ID)],
    imputed_only[, .(CHROM, POS, REF, ALT, ID_imputed = ID)]
  ), use.names = TRUE)
  tmeta <- rbindlist(list(
    truth[, .(CHROM, POS, REF, ALT, ID_truth = ID)],
    truth_only[, .(CHROM, POS, REF, ALT, ID_truth = ID)]
  ), use.names = TRUE)
  allele_pairs <- merge(
    imeta, tmeta, by = c("CHROM", "POS"), all = FALSE,
    suffixes = c("_imputed", "_truth"), allow.cartesian = TRUE
  )
  allele_mismatches <- allele_pairs[REF_imputed != REF_truth | ALT_imputed != ALT_truth]
  setorder(allele_mismatches, CHROM, POS, REF_imputed, ALT_imputed, REF_truth, ALT_truth)

  retained <- which(site_pass)
  meta <- imputed[retained, ..variant_cols]
  meta[, ID_truth := truth$ID[retained]]

  truth_gt_count <- matrix(NA_real_, nrow = length(retained), ncol = length(samples), dimnames = list(NULL, samples))
  imputed_gt_count <- matrix(NA_real_, nrow = length(retained), ncol = length(samples), dimnames = list(NULL, samples))
  masking_rows <- vector("list", length(samples) * 9L)
  masking_index <- 0L

  for (j in seq_along(samples)) {
    sample <- samples[[j]]
    truth_decoded <- decode_gt(truth[[paste0(sample, ".GT")]][retained])
    truth_invalid_gt <- truth_decoded$missing | truth_decoded$invalid
    truth_count_sample <- truth_decoded$dosage

    if (filter_enabled) {
      gq <- as_num(truth[[paste0(sample, ".GQ")]][retained])
      dp <- as_num(truth[[paste0(sample, ".DP")]][retained])
      missing_gq <- is.na(gq)
      low_gq <- !missing_gq & gq < thresholds$min_gq
      missing_dp <- is.na(dp)
      low_dp <- !missing_dp & dp < thresholds$min_dp
    } else {
      missing_gq <- low_gq <- missing_dp <- low_dp <- rep(FALSE, length(retained))
    }
    truth_mask <- truth_invalid_gt | missing_gq | low_gq | missing_dp | low_dp
    truth_count_sample[truth_mask] <- NA_real_
    truth_gt_count[, j] <- truth_count_sample

    imputed_decoded <- decode_gt(imputed[[sample]][retained])
    imputed_count_sample <- imputed_decoded$dosage
    imputed_count_sample[imputed_decoded$missing | imputed_decoded$invalid] <- NA_real_
    imputed_gt_count[, j] <- imputed_count_sample

    usable <- is.finite(truth_count_sample) & is.finite(imputed_count_sample)
    reason_vectors <- list(
      imputed_missing_GT = imputed_decoded$missing,
      imputed_invalid_GT = imputed_decoded$invalid,
      truth_invalid_GT = truth_invalid_gt,
      truth_missing_GQ = missing_gq,
      truth_GQ_below_min = low_gq,
      truth_missing_DP = missing_dp,
      truth_DP_below_min = low_dp,
      truth_masked_any = truth_mask,
      usable_pair = usable
    )
    for (reason in names(reason_vectors)) {
      masking_index <- masking_index + 1L
      masking_rows[[masking_index]] <- data.table(
        sample = sample, chromosome = chromosome, reason = reason,
        count = sum(reason_vectors[[reason]])
      )
    }
  }
  masking_summary <- rbindlist(masking_rows[seq_len(masking_index)])

  truth_valid <- is.finite(truth_gt_count)
  truth_count <- rowSums(truth_valid)
  truth_zero <- truth_gt_count
  truth_zero[!truth_valid] <- 0
  truth_af <- rep(NA_real_, nrow(meta))
  has_truth <- truth_count > 0
  truth_af[has_truth] <- rowSums(truth_zero)[has_truth] / (2 * truth_count[has_truth])
  truth_maf <- ifelse(is.na(truth_af), NA_real_, pmin(truth_af, 1 - truth_af))
  maf_bin <- cut(truth_maf, breaks = maf_breaks, labels = maf_labels, include.lowest = TRUE, right = FALSE)

  pair_ok <- is.finite(imputed_gt_count) & is.finite(truth_gt_count)
  x <- imputed_gt_count
  y <- truth_gt_count
  x[!pair_ok] <- 0
  y[!pair_ok] <- 0
  n_pairs <- rowSums(pair_ok)
  sx <- rowSums(x)
  sy <- rowSums(y)
  sxx <- rowSums(x * x)
  syy <- rowSums(y * y)
  sxy <- rowSums(x * y)
  variant_cor <- cor_from_sums(n_pairs, sx, sy, sxx, syy, sxy, 3L)
  n_concordant <- rowSums(pair_ok & imputed_gt_count == truth_gt_count, na.rm = TRUE)
  concordance <- rep(NA_real_, nrow(meta))
  concordance[n_pairs > 0] <- n_concordant[n_pairs > 0] / n_pairs[n_pairs > 0]
  variant_metrics <- copy(meta)
  variant_metrics[, `:=`(
    n_pairs = as.integer(n_pairs),
    truth_alt_frequency = truth_af,
    truth_maf = truth_maf,
    r = variant_cor$r,
    r2 = variant_cor$r2,
    concordance = concordance
  )]

  stats_rows <- vector("list", length(samples) * (length(maf_labels) + 1L))
  stats_index <- 0L
  for (j in seq_along(samples)) {
    for (bin in c("ALL", maf_labels)) {
      in_bin <- if (bin == "ALL") rep(TRUE, nrow(meta)) else !is.na(maf_bin) & maf_bin == bin
      ok <- pair_ok[, j] & in_bin
      xv <- imputed_gt_count[ok, j]
      yv <- truth_gt_count[ok, j]
      stats_index <- stats_index + 1L
      stats_rows[[stats_index]] <- data.table(
        sample = samples[[j]], chromosome = chromosome, maf_bin = bin,
        n = length(xv), sx = sum(xv), sy = sum(yv),
        sxx = sum(xv * xv), syy = sum(yv * yv), sxy = sum(xv * yv)
      )
    }
  }
  sample_stats <- rbindlist(stats_rows)

  variant_counts <- data.table(
    section = "variants", chromosome = chromosome,
    reason = c(
      "imputed_structural", "truth_structural", "exact_matches", "retained_exact_matches",
      "imputed_only", "truth_only", "allele_mismatch_pairs"
    ),
    count = c(
      nrow(imputed) + nrow(imputed_only), nrow(truth) + nrow(truth_only),
      nrow(imputed), length(retained), nrow(imputed_only), nrow(truth_only),
      nrow(allele_mismatches)
    ),
    value = NA_character_
  )
  site_counts <- data.table(
    section = "site_filter", chromosome = chromosome,
    reason = c("sites_total", "sites_passed", "sites_failed", names(site_fail)),
    count = c(nrow(truth), sum(site_pass), sum(!site_pass), vapply(site_fail, sum, numeric(1))),
    value = NA_character_
  )
  filter_summary <- rbindlist(list(variant_counts, site_counts), use.names = TRUE)

  out_dir <- opt[["out-dir"]]
  atomic_parquet(variant_metrics, file.path(out_dir, "metrics", "per_variant_metrics"), chromosome)
  atomic_parquet(site_filtered_variants, file.path(out_dir, "metrics", "site_filtered_variants"), chromosome)
  atomic_parquet(imputed_only[, ..variant_cols], file.path(out_dir, "metrics", "imputed_only_variants"), chromosome)
  atomic_parquet(truth_only[, ..variant_cols], file.path(out_dir, "metrics", "truth_only_variants"), chromosome)
  atomic_parquet(allele_mismatches, file.path(out_dir, "metrics", "allele_mismatches"), chromosome)
  atomic_fwrite(sample_stats, file.path(out_dir, "intermediate", "chromosome_stats", paste0(chromosome, ".sample_stats.tsv")))
  atomic_fwrite(filter_summary, file.path(out_dir, "qc", "chromosomes", paste0(chromosome, ".filter_summary.tsv")))
  atomic_fwrite(masking_summary, file.path(out_dir, "qc", "chromosomes", paste0(chromosome, ".genotype_masking_summary.tsv")))

  cat(sprintf("[INFO] WGS GT chromosome metrics complete: %s, %d retained variants, %d samples\n", chromosome, nrow(meta), length(samples)))
}

run_finalize <- function(opt) {
  require_options(opt, c(
    "out-dir", "chromosome-manifest", "filter-enabled", "min-qual", "min-qd",
    "max-sor", "max-fs", "min-mq", "min-mq-rank-sum",
    "min-read-pos-rank-sum", "min-gq", "min-dp"
  ))
  out_dir <- opt[["out-dir"]]
  tasks <- fread(opt[["chromosome-manifest"]], sep = "\t", header = TRUE)
  if (!nrow(tasks) || anyDuplicated(tasks$chromosome)) stop("Invalid chromosome task manifest.", call. = FALSE)

  stats <- rbindlist(lapply(tasks$chromosome, function(chromosome) {
    fread(file.path(out_dir, "intermediate", "chromosome_stats", paste0(chromosome, ".sample_stats.tsv")))
  }))
  sum_cols <- c("n", "sx", "sy", "sxx", "syy", "sxy")
  combined <- stats[, lapply(.SD, sum), by = .(sample, maf_bin), .SDcols = sum_cols]
  cors <- cor_from_sums(combined$n, combined$sx, combined$sy, combined$sxx, combined$syy, combined$sxy, 2L)
  combined[, `:=`(r = cors$r, r2 = cors$r2)]

  sample_order <- unique(stats$sample)
  sample_rows <- lapply(sample_order, function(sample_id) {
    x <- combined[sample == sample_id]
    overall <- x[maf_bin == "ALL"]
    row <- data.table(
      sample = sample_id,
      r_overall = overall$r,
      r2_overall = overall$r2,
      n_variants = as.integer(overall$n)
    )
    for (bin in maf_labels) {
      bx <- x[maf_bin == bin]
      set(row, j = paste0("r_maf_", bin), value = bx$r)
      set(row, j = paste0("r2_maf_", bin), value = bx$r2)
      set(row, j = paste0("n_maf_", bin), value = as.integer(bx$n))
    }
    row
  })
  sample_metrics <- rbindlist(sample_rows)
  setorder(sample_metrics, r2_overall)
  atomic_fwrite(sample_metrics, file.path(out_dir, "per_sample_metrics.tsv"))

  filter_base <- rbindlist(lapply(tasks$chromosome, function(chromosome) {
    fread(file.path(out_dir, "qc", "chromosomes", paste0(chromosome, ".filter_summary.tsv")))
  }), fill = TRUE)
  filter_all <- filter_base[, .(count = sum(count), value = NA_character_), by = .(section, reason)]
  filter_all[, chromosome := "ALL"]
  setcolorder(filter_all, names(filter_base))
  config_values <- c(
    filter_enabled = opt[["filter-enabled"]],
    min_qual = opt[["min-qual"]],
    min_qd = opt[["min-qd"]], max_sor = opt[["max-sor"]], max_fs = opt[["max-fs"]],
    min_mq = opt[["min-mq"]], min_mq_rank_sum = opt[["min-mq-rank-sum"]],
    min_read_pos_rank_sum = opt[["min-read-pos-rank-sum"]], min_gq = opt[["min-gq"]],
    min_dp = opt[["min-dp"]]
  )
  config_rows <- data.table(
    section = "config", chromosome = "ALL", reason = names(config_values),
    count = NA_real_, value = unname(config_values)
  )
  filter_summary <- rbindlist(list(config_rows, filter_base, filter_all), use.names = TRUE, fill = TRUE)
  atomic_fwrite(filter_summary, file.path(out_dir, "qc", "filter_summary.tsv"))

  masking_base <- rbindlist(lapply(tasks$chromosome, function(chromosome) {
    fread(file.path(out_dir, "qc", "chromosomes", paste0(chromosome, ".genotype_masking_summary.tsv")))
  }))
  masking_sample_all <- masking_base[, .(count = sum(count)), by = .(sample, reason)][, chromosome := "ALL"]
  setcolorder(masking_sample_all, names(masking_base))
  masking_chr_all <- masking_base[, .(count = sum(count)), by = .(chromosome, reason)][, sample := "ALL"]
  setcolorder(masking_chr_all, names(masking_base))
  masking_all <- masking_base[, .(count = sum(count)), by = reason][, `:=`(sample = "ALL", chromosome = "ALL")]
  setcolorder(masking_all, names(masking_base))
  masking_summary <- rbindlist(list(masking_base, masking_sample_all, masking_chr_all, masking_all))
  atomic_fwrite(masking_summary, file.path(out_dir, "qc", "genotype_masking_summary.tsv"))

  cat(sprintf("[INFO] WGS GT finalization complete: %d chromosomes, %d samples\n", nrow(tasks), length(sample_order)))
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
mode <- opt[["mode"]]
if (is.null(mode) || !mode %chin% c("chromosome", "finalize")) {
  stop("--mode must be chromosome or finalize.", call. = FALSE)
}
if (mode == "chromosome") run_chromosome(opt) else run_finalize(opt)
