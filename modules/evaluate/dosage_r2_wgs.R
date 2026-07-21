#!/usr/bin/env Rscript

# Compute DS-versus-GT dosage metrics from normalized, biallelic SNP tables.
# Input extraction and VCF normalization are handled by dosage_r2_wgs.sh.

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) {
      stop("Invalid argument list near: ", key, call. = FALSE)
    }
    out[[substring(key, 3L)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

required_args <- c(
  "imputed-raw", "truth-raw", "samples", "out-dir", "filter-enabled",
  "min-qual", "min-qd", "max-sor", "max-fs", "min-mq",
  "min-mq-rank-sum", "min-read-pos-rank-sum", "min-gq", "min-dp"
)
opt <- parse_args(commandArgs(trailingOnly = TRUE))
missing_args <- required_args[!required_args %in% names(opt)]
if (length(missing_args)) {
  stop("Missing arguments: ", paste(missing_args, collapse = ", "), call. = FALSE)
}

filter_enabled <- identical(opt[["filter-enabled"]], "true")
thresholds <- c(
  min_qual = as.numeric(opt[["min-qual"]]),
  min_qd = as.numeric(opt[["min-qd"]]),
  max_sor = as.numeric(opt[["max-sor"]]),
  max_fs = as.numeric(opt[["max-fs"]]),
  min_mq = as.numeric(opt[["min-mq"]]),
  min_mq_rank_sum = as.numeric(opt[["min-mq-rank-sum"]]),
  min_read_pos_rank_sum = as.numeric(opt[["min-read-pos-rank-sum"]]),
  min_gq = as.numeric(opt[["min-gq"]]),
  min_dp = as.numeric(opt[["min-dp"]])
)

samples <- readLines(opt$samples, warn = FALSE)
samples <- samples[nzchar(samples)]
if (!length(samples) || anyDuplicated(samples)) {
  stop("The sample file must contain at least one unique sample ID.", call. = FALSE)
}

read_tsv <- function(path) {
  read.delim(
    path, header = TRUE, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, na.strings = character()
  )
}

imputed <- read_tsv(opt[["imputed-raw"]])
truth <- read_tsv(opt[["truth-raw"]])
variant_cols <- c("CHROM", "POS", "REF", "ALT", "ID")
truth_site_cols <- c("QUAL", "QD", "SOR", "FS", "MQ", "MQRankSum", "ReadPosRankSum")

missing_imputed <- setdiff(c(variant_cols, samples), names(imputed))
missing_truth <- setdiff(c(variant_cols, truth_site_cols, paste0(samples, ".GT")), names(truth))
if (length(missing_imputed)) {
  stop("Missing imputed columns: ", paste(missing_imputed, collapse = ", "), call. = FALSE)
}
if (length(missing_truth)) {
  stop("Missing truth columns: ", paste(missing_truth, collapse = ", "), call. = FALSE)
}

make_key <- function(x) paste(x$CHROM, x$POS, x$REF, x$ALT, sep = ":")
make_pos_key <- function(x) paste(x$CHROM, x$POS, sep = ":")
imputed$key <- make_key(imputed)
truth$key <- make_key(truth)
if (anyDuplicated(imputed$key)) stop("Duplicate exact variants in normalized imputed input.", call. = FALSE)
if (anyDuplicated(truth$key)) stop("Duplicate exact variants in normalized truth input.", call. = FALSE)

chrom_number <- function(x) suppressWarnings(as.integer(sub("^Chr", "", x)))
variant_order <- function(x) order(chrom_number(x$CHROM), as.numeric(x$POS), x$REF, x$ALT, na.last = TRUE)
imputed <- imputed[variant_order(imputed), , drop = FALSE]
truth <- truth[variant_order(truth), , drop = FALSE]

as_num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("", "."), NA, x)))
qual <- as_num(truth$QUAL)
qd <- as_num(truth$QD)
sor <- as_num(truth$SOR)
fs <- as_num(truth$FS)
mq <- as_num(truth$MQ)
mq_rank <- as_num(truth$MQRankSum)
read_pos_rank <- as_num(truth$ReadPosRankSum)

site_fail <- list(
  missing_QUAL = is.na(qual),
  QUAL_below_min = !is.na(qual) & qual < thresholds[["min_qual"]],
  missing_QD = is.na(qd),
  QD_below_min = !is.na(qd) & qd < thresholds[["min_qd"]],
  SOR_above_max = !is.na(sor) & sor > thresholds[["max_sor"]],
  FS_above_max = !is.na(fs) & fs > thresholds[["max_fs"]],
  missing_MQ = is.na(mq),
  MQ_below_min = !is.na(mq) & mq < thresholds[["min_mq"]],
  MQRankSum_below_min = !is.na(mq_rank) & mq_rank < thresholds[["min_mq_rank_sum"]],
  ReadPosRankSum_below_min = !is.na(read_pos_rank) & read_pos_rank < thresholds[["min_read_pos_rank_sum"]]
)
if (filter_enabled) {
  site_pass <- !Reduce(`|`, site_fail)
} else {
  site_pass <- rep(TRUE, nrow(truth))
  site_fail <- lapply(site_fail, function(x) rep(FALSE, length(x)))
}
truth$site_filter_pass <- site_pass
truth$site_filter_reasons <- vapply(seq_len(nrow(truth)), function(i) {
  failed <- names(site_fail)[vapply(site_fail, `[[`, logical(1), i)]
  if (length(failed)) paste(failed, collapse = ";") else "PASS"
}, character(1))

imputed_pos <- unique(imputed[c("CHROM", "POS")])
truth_pos <- unique(truth[c("CHROM", "POS")])
shared_pos <- merge(imputed_pos, truth_pos, by = c("CHROM", "POS"))
imputed_shared <- merge(shared_pos, imputed[c("CHROM", "POS", "REF", "ALT", "ID")], by = c("CHROM", "POS"))
truth_shared <- merge(shared_pos, truth[c("CHROM", "POS", "REF", "ALT", "ID")], by = c("CHROM", "POS"))
allele_pairs <- merge(imputed_shared, truth_shared, by = c("CHROM", "POS"), suffixes = c("_imputed", "_truth"))
allele_mismatches <- allele_pairs[
  allele_pairs$REF_imputed != allele_pairs$REF_truth | allele_pairs$ALT_imputed != allele_pairs$ALT_truth,
  , drop = FALSE
]

exact_keys <- intersect(imputed$key, truth$key)
matched_i <- match(exact_keys, imputed$key)
matched_t <- match(exact_keys, truth$key)
matched_variants <- data.frame(
  CHROM = imputed$CHROM[matched_i], POS = imputed$POS[matched_i],
  REF = imputed$REF[matched_i], ALT = imputed$ALT[matched_i],
  ID_imputed = imputed$ID[matched_i], ID_truth = truth$ID[matched_t],
  site_filter_pass = truth$site_filter_pass[matched_t],
  site_filter_reasons = truth$site_filter_reasons[matched_t],
  stringsAsFactors = FALSE
)
if (nrow(matched_variants)) matched_variants <- matched_variants[variant_order(matched_variants), , drop = FALSE]

imputed_only <- imputed[!imputed$key %in% truth$key, variant_cols, drop = FALSE]
truth_only <- truth[!truth$key %in% imputed$key, variant_cols, drop = FALSE]

retained_keys <- truth$key[truth$site_filter_pass & truth$key %in% imputed$key]
retained_i <- match(retained_keys, imputed$key)
retained_t <- match(retained_keys, truth$key)
if (length(retained_keys)) {
  ord <- variant_order(imputed[retained_i, , drop = FALSE])
  retained_keys <- retained_keys[ord]
  retained_i <- retained_i[ord]
  retained_t <- retained_t[ord]
}
meta <- imputed[retained_i, variant_cols, drop = FALSE]

truth_dosage <- matrix(NA_real_, nrow = length(retained_keys), ncol = length(samples), dimnames = list(NULL, samples))
imputed_dosage <- truth_dosage
mask_rows <- list()

valid_gt_pattern <- "^(0[/|]0|0[/|]1|1[/|]0|1[/|]1)$"
gt_to_dosage <- function(gt) {
  out <- rep(NA_real_, length(gt))
  out[gt %in% c("0/0", "0|0")] <- 0
  out[gt %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  out[gt %in% c("1/1", "1|1")] <- 2
  out
}

for (j in seq_along(samples)) {
  sample <- samples[[j]]
  gt <- truth[[paste0(sample, ".GT")]][retained_t]
  invalid_gt <- is.na(gt) | gt %in% c("", ".", "./.", ".|.") | !grepl(valid_gt_pattern, gt)

  if (filter_enabled) {
    gq_col <- paste0(sample, ".GQ")
    dp_col <- paste0(sample, ".DP")
    if (!gq_col %in% names(truth) || !dp_col %in% names(truth)) {
      stop("Truth GQ/DP columns are required when filtering is enabled.", call. = FALSE)
    }
    gq <- as_num(truth[[gq_col]][retained_t])
    dp <- as_num(truth[[dp_col]][retained_t])
    missing_gq <- is.na(gq)
    low_gq <- !missing_gq & gq < thresholds[["min_gq"]]
    missing_dp <- is.na(dp)
    low_dp <- !missing_dp & dp < thresholds[["min_dp"]]
  } else {
    missing_gq <- low_gq <- missing_dp <- low_dp <- rep(FALSE, length(retained_t))
  }

  truth_mask <- invalid_gt | missing_gq | low_gq | missing_dp | low_dp
  td <- gt_to_dosage(gt)
  td[truth_mask] <- NA_real_
  truth_dosage[, j] <- td

  ds_raw <- imputed[[sample]][retained_i]
  ds_missing <- is.na(ds_raw) | ds_raw %in% c("", ".")
  ds <- as_num(ds_raw)
  ds_invalid <- !ds_missing & (is.na(ds) | !is.finite(ds) | ds < 0 | ds > 2)
  ds[ds_missing | ds_invalid] <- NA_real_
  imputed_dosage[, j] <- ds

  usable <- !is.na(td) & !is.na(ds)
  reason_vectors <- list(
    truth_invalid_GT = invalid_gt,
    truth_missing_GQ = missing_gq,
    truth_GQ_below_min = low_gq,
    truth_missing_DP = missing_dp,
    truth_DP_below_min = low_dp,
    truth_masked_any = truth_mask,
    imputed_missing_DS = ds_missing,
    imputed_invalid_DS = ds_invalid,
    usable_pair = usable
  )
  chromosomes <- unique(c(as.character(meta$CHROM), "ALL"))
  for (chrom in chromosomes) {
    in_chrom <- if (chrom == "ALL") rep(TRUE, nrow(meta)) else meta$CHROM == chrom
    for (reason in names(reason_vectors)) {
      mask_rows[[length(mask_rows) + 1L]] <- data.frame(
        sample = sample, chromosome = chrom, reason = reason,
        count = sum(reason_vectors[[reason]] & in_chrom), stringsAsFactors = FALSE
      )
    }
  }
}

masking_summary <- if (length(mask_rows)) do.call(rbind, mask_rows) else data.frame(
  sample = character(), chromosome = character(), reason = character(), count = integer()
)
if (nrow(masking_summary)) {
  non_all_rows <- masking_summary[masking_summary$chromosome != "ALL", , drop = FALSE]
  if (nrow(non_all_rows)) {
    aggregate_rows <- aggregate(count ~ chromosome + reason, data = non_all_rows, sum)
    aggregate_rows$sample <- "ALL"
    aggregate_rows <- aggregate_rows[c("sample", "chromosome", "reason", "count")]
  } else {
    aggregate_rows <- masking_summary[FALSE, , drop = FALSE]
  }
  all_rows <- aggregate(count ~ reason, data = masking_summary[masking_summary$chromosome == "ALL", ], sum)
  all_rows$chromosome <- "ALL"
  all_rows$sample <- "ALL"
  all_rows <- all_rows[c("sample", "chromosome", "reason", "count")]
  masking_summary <- rbind(masking_summary, aggregate_rows, all_rows)
}

safe_cor <- function(x, y, min_pairs) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < min_pairs || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) {
    return(c(n = n, r = NA_real_, r2 = NA_real_))
  }
  r <- suppressWarnings(cor(x[ok], y[ok], method = "pearson"))
  c(n = n, r = r, r2 = r * r)
}

variant_metrics <- data.frame(
  CHROM = meta$CHROM, POS = meta$POS, REF = meta$REF, ALT = meta$ALT, ID = meta$ID,
  n_pairs = integer(nrow(meta)), truth_alt_frequency = numeric(nrow(meta)),
  truth_maf = numeric(nrow(meta)), pearson_r = numeric(nrow(meta)), dosage_r2 = numeric(nrow(meta)),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(meta))) {
  metric <- safe_cor(imputed_dosage[i, ], truth_dosage[i, ], 3L)
  truth_values <- truth_dosage[i, is.finite(truth_dosage[i, ])]
  af <- if (length(truth_values)) mean(truth_values) / 2 else NA_real_
  variant_metrics$n_pairs[i] <- metric[["n"]]
  variant_metrics$truth_alt_frequency[i] <- af
  variant_metrics$truth_maf[i] <- if (is.na(af)) NA_real_ else min(af, 1 - af)
  variant_metrics$pearson_r[i] <- metric[["r"]]
  variant_metrics$dosage_r2[i] <- metric[["r2"]]
}

sample_metrics <- data.frame(
  sample = samples, n_pairs = integer(length(samples)), pearson_r = numeric(length(samples)),
  dosage_r2 = numeric(length(samples)), stringsAsFactors = FALSE
)
for (j in seq_along(samples)) {
  metric <- safe_cor(imputed_dosage[, j], truth_dosage[, j], 2L)
  sample_metrics$n_pairs[j] <- metric[["n"]]
  sample_metrics$pearson_r[j] <- metric[["r"]]
  sample_metrics$dosage_r2[j] <- metric[["r2"]]
}

config_values <- c(
  filter_enabled = opt[["filter-enabled"]], min_qual = opt[["min-qual"]],
  min_qd = opt[["min-qd"]], max_sor = opt[["max-sor"]], max_fs = opt[["max-fs"]],
  min_mq = opt[["min-mq"]], min_mq_rank_sum = opt[["min-mq-rank-sum"]],
  min_read_pos_rank_sum = opt[["min-read-pos-rank-sum"]], min_gq = opt[["min-gq"]],
  min_dp = opt[["min-dp"]]
)
filter_rows <- lapply(names(config_values), function(name) data.frame(
  section = "config", chromosome = "ALL", reason = name, count = NA_integer_,
  value = config_values[[name]], stringsAsFactors = FALSE
))
filter_rows[[length(filter_rows) + 1L]] <- data.frame(
  section = "variants", chromosome = "ALL", reason = c(
    "imputed_structural", "truth_structural", "exact_matches", "retained_exact_matches",
    "imputed_only", "truth_only", "allele_mismatch_pairs"
  ),
  count = c(nrow(imputed), nrow(truth), length(exact_keys), length(retained_keys),
            nrow(imputed_only), nrow(truth_only), nrow(allele_mismatches)),
  value = NA_character_, stringsAsFactors = FALSE
)
site_chromosomes <- unique(c(as.character(truth$CHROM), "ALL"))
for (chrom in site_chromosomes) {
  in_chrom <- if (chrom == "ALL") rep(TRUE, nrow(truth)) else truth$CHROM == chrom
  filter_rows[[length(filter_rows) + 1L]] <- data.frame(
    section = "site_filter", chromosome = chrom,
    reason = c("sites_total", "sites_passed", "sites_failed"),
    count = c(sum(in_chrom), sum(site_pass & in_chrom), sum(!site_pass & in_chrom)),
    value = NA_character_, stringsAsFactors = FALSE
  )
  for (reason in names(site_fail)) {
    filter_rows[[length(filter_rows) + 1L]] <- data.frame(
      section = "site_filter", chromosome = chrom, reason = reason,
      count = sum(site_fail[[reason]] & in_chrom), value = NA_character_, stringsAsFactors = FALSE
    )
  }
}
filter_summary <- do.call(rbind, filter_rows)

dir.create(opt[["out-dir"]], recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt[["out-dir"]], "intermediate"), recursive = TRUE, showWarnings = FALSE)
write_out <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_out(variant_metrics, file.path(opt[["out-dir"]], "per_variant_metrics.tsv"))
write_out(sample_metrics, file.path(opt[["out-dir"]], "per_sample_metrics.tsv"))
write_out(filter_summary, file.path(opt[["out-dir"]], "filter_summary.tsv"))
write_out(masking_summary, file.path(opt[["out-dir"]], "genotype_masking_summary.tsv"))
write_out(matched_variants, file.path(opt[["out-dir"]], "matched_variants.tsv"))
write_out(imputed_only[variant_order(imputed_only), , drop = FALSE], file.path(opt[["out-dir"]], "imputed_only_variants.tsv"))
write_out(truth_only[variant_order(truth_only), , drop = FALSE], file.path(opt[["out-dir"]], "truth_only_variants.tsv"))
write_out(allele_mismatches, file.path(opt[["out-dir"]], "allele_mismatches.tsv"))
write_out(cbind(meta, as.data.frame(imputed_dosage, check.names = FALSE)), file.path(opt[["out-dir"]], "intermediate", "imputed_ds.tsv"))
write_out(cbind(meta, as.data.frame(truth_dosage, check.names = FALSE)), file.path(opt[["out-dir"]], "intermediate", "truth_gt_dosage.tsv"))

cat(sprintf("[INFO] WGS metrics complete: %d retained variants, %d samples\n", nrow(meta), length(samples)))
