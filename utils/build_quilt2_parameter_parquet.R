#!/usr/bin/env Rscript
# Agreed compromise: score reported GT against truth regardless of GP confidence.
# Write the audited locus/call schemas. GP and maxGP are retained for auditing;
# no GP-argmax correctness or calibration is derived. FORMAT/HD is omitted.
# Final loci require evaluable GT-truth comparisons for all 25 targets in all six conditions.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

fail <- function(...) stop(paste0(...), call. = FALSE)

analysis_compromise <- paste(
  "Only reported GT correctness is evaluated, regardless of GP confidence.",
  "GP bins, GP-argmax correctness, and probability calibration are not assessed.",
  "Missing, malformed, or tied GP does not exclude an otherwise evaluable GT call.",
  "The final mask requires valid aligned truth and reported GT for every target in every run."
)

# Base order follows Locus_mock and Calls_mock in quilt2_parquet_schema_mock.xlsx.
# Calls additionally retain same-run INFO_SCORE, HWE, EAF, ERC, EAC, PAF and
# INFO_SCORE/HWE validity flags for call-level analysis without a locus join.
# Additional emitted INFO/FORMAT fields are retained at their own resolution,
# except FORMAT/HD. Internal validity checks and provenance are not row columns.
locus_columns <- c(
  "source_id", "chrom", "pos", "ref", "alt", "treatment", "panel", "truth_source",
  "info_score", "eaf", "hwe", "erc", "eac", "paf",
  "truth_vcf_id", "truth_ref", "truth_alt", "truth_qual", "allele_alignment_status",
  "info_score_valid", "hwe_valid", "n_correct_source", "percent_correct_source",
  "n_correct_panel", "percent_correct_panel", "n_correct_pooled", "percent_correct_pooled",
  "mean_info_score", "truth_alt_allele_count", "truth_maf"
)
call_columns <- c(
  "source_id", "chrom", "pos", "ref", "alt", "treatment", "panel", "truth_source",
  "sample_id", "array_group", "fmt__gt", "fmt__gp", "fmt__ds", "truth_gt_raw",
  "truth_gq", "truth_dp", "imputed_gt_dosage", "truth_dosage_aligned", "gt_correct",
  "gp_valid", "maxGP", "info_score", "info_score_valid",
  "hwe", "eaf", "erc", "eac", "paf", "hwe_valid"
)
join_columns <- c("chrom", "pos", "ref", "alt", "source_id")
schema_version <- "quilt2-parameter-validation-v8"

parse_args <- function(args) {
  value_options <- c(
    "--staged-manifest", "--truth-manifest", "--eval-mask-manifest", "--eval-mask",
    "--sample-map", "--array-truth", "--wgs-truth",
    "--locus-output", "--calls-output", "--mask-out", "--summary-out", "--min-gq", "--min-dp",
    "--gp-sum-tolerance"
  )
  opts <- list(
    staged_manifest = NULL,
    truth_manifest = NULL,
    eval_mask_manifest = NULL,
    eval_mask = NULL,
    sample_map = NULL,
    array_truth = NULL,
    wgs_truth = NULL,
    locus_output = NULL,
    calls_output = NULL,
    mask_out = NULL,
    summary_out = NULL,
    min_gq = 60,
    min_dp = 10,
    gp_sum_tolerance = 0.001
  )
  key_map <- c(
    "--staged-manifest" = "staged_manifest",
    "--truth-manifest" = "truth_manifest",
    "--eval-mask-manifest" = "eval_mask_manifest",
    "--eval-mask" = "eval_mask",
    "--sample-map" = "sample_map",
    "--array-truth" = "array_truth",
    "--wgs-truth" = "wgs_truth",
    "--locus-output" = "locus_output",
    "--calls-output" = "calls_output",
    "--mask-out" = "mask_out",
    "--summary-out" = "summary_out",
    "--min-gq" = "min_gq",
    "--min-dp" = "min_dp",
    "--gp-sum-tolerance" = "gp_sum_tolerance"
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      cat(paste(
        "Usage: build_quilt2_parameter_parquet.R",
        "--staged-manifest FILE --truth-manifest FILE",
        "--eval-mask-manifest FILE --eval-mask FILE --sample-map FILE",
        "--array-truth VCF --wgs-truth VCF --locus-output FILE --calls-output FILE",
        "--mask-out FILE --summary-out FILE",
        "[--min-gq 60] [--min-dp 10] [--gp-sum-tolerance 0.001]\n"
      ))
      quit(status = 0)
    }
    if (!key %in% value_options || i == length(args)) {
      fail("Unknown or incomplete option: ", key)
    }
    opts[[key_map[[key]]]] <- args[[i + 1L]]
    i <- i + 2L
  }

  required <- c(
    "staged_manifest", "truth_manifest", "eval_mask_manifest", "eval_mask",
    "sample_map", "array_truth", "wgs_truth",
    "locus_output", "calls_output", "mask_out", "summary_out"
  )
  missing <- required[vapply(required, function(x) is.null(opts[[x]]), logical(1))]
  if (length(missing)) fail("Missing required option(s): ", paste(missing, collapse = ", "))

  opts$min_gq <- suppressWarnings(as.numeric(opts$min_gq))
  opts$min_dp <- suppressWarnings(as.numeric(opts$min_dp))
  opts$gp_sum_tolerance <- suppressWarnings(as.numeric(opts$gp_sum_tolerance))
  if (!is.finite(opts$min_gq) || opts$min_gq < 0) fail("--min-gq must be non-negative")
  if (!is.finite(opts$min_dp) || opts$min_dp < 0) fail("--min-dp must be non-negative")
  if (!is.finite(opts$gp_sum_tolerance) || opts$gp_sum_tolerance <= 0) {
    fail("--gp-sum-tolerance must be positive")
  }
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

for (path in c(
  opts$staged_manifest, opts$truth_manifest, opts$eval_mask_manifest, opts$eval_mask,
  opts$sample_map, opts$array_truth, opts$wgs_truth
)) {
  if (!file.exists(path)) fail("Input not found: ", path)
}
if (opts$locus_output == opts$calls_output) fail("Locus and call outputs must be different files")
for (path in c(opts$locus_output, opts$calls_output)) {
  if (file.exists(path)) fail("Refusing to overwrite existing output: ", path)
}

read_tabular <- function(path, label) {
  result <- fread(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    na.strings = character(),
    showProgress = FALSE
  )
  if (!nrow(result)) fail(label, " is empty: ", path)
  result
}

manifest <- read_tabular(opts$staged_manifest, "Staged manifest")
manifest_required <- c(
  "source_id", "logical_run_id", "treatment", "panel", "truth_source",
  "input_type", "vcf", "source_path", "source_sample_count", "selected_sample_count",
  "source_samples_sha256", "source_signature_sha256", "candidate_vcf_sha256",
  "candidate_header_sha256"
)
manifest_missing <- setdiff(manifest_required, names(manifest))
if (length(manifest_missing)) {
  fail("Staged manifest is missing: ", paste(manifest_missing, collapse = ", "))
}
manifest <- manifest[, ..manifest_required]
for (column in manifest_required) set(manifest, j = column, value = trimws(as.character(manifest[[column]])))
if (anyNA(manifest) || any(manifest == "")) fail("Staged manifest contains blank values")
if (anyDuplicated(manifest$source_id)) fail("Staged manifest source_id values are not unique")
if (any(!manifest$treatment %in% c("Filtered", "No_filter"))) fail("Unexpected treatment")
if (any(!manifest$panel %in% c("Liao", "NCBI", "Combined"))) fail("Unexpected panel")
if (any(!manifest$truth_source %in% c("array", "wgs"))) fail("Unexpected truth_source")
if (any(!manifest$input_type %in% c("vcf", "chromosome_vcfs"))) {
  fail("Unexpected input_type; expected vcf or chromosome_vcfs")
}
if (nrow(manifest) != 12L) fail("Expected 12 staged sources; found ", nrow(manifest))
hash_columns <- c(
  "source_samples_sha256", "source_signature_sha256", "candidate_vcf_sha256",
  "candidate_header_sha256"
)
if (any(!vapply(
  unlist(manifest[, ..hash_columns], use.names = FALSE),
  function(value) grepl("^[0-9a-fA-F]{64}$", value),
  logical(1)
))) fail("Staged manifest contains an invalid SHA-256 value")

expected_conditions <- CJ(
  treatment = c("Filtered", "No_filter"),
  panel = c("Liao", "NCBI", "Combined"),
  truth_source = c("array", "wgs"),
  unique = TRUE
)
condition_key <- function(x) paste(x$treatment, x$panel, x$truth_source, sep = "\r")
if (!setequal(condition_key(manifest), condition_key(expected_conditions))) {
  fail("Staged manifest must contain each treatment x panel x truth_source combination exactly once")
}
if (anyDuplicated(condition_key(manifest))) fail("Staged manifest contains duplicate conditions")
logical_runs <- unique(manifest[, .(logical_run_id, treatment, panel)])
logical_run_counts <- logical_runs[, .N, by = .(treatment, panel)]
if (nrow(logical_runs) != 6L || any(logical_run_counts$N != 1L)) {
  fail("Each treatment x panel condition must have exactly one logical_run_id shared by Array and WGS")
}
if (anyDuplicated(logical_runs$logical_run_id)) fail("logical_run_id values must identify one condition each")
if (any(!file.exists(manifest$vcf))) fail("One or more staged VCFs are missing")
manifest[, `:=`(
  source_sample_count = suppressWarnings(as.integer(source_sample_count)),
  selected_sample_count = suppressWarnings(as.integer(selected_sample_count))
)]
if (anyNA(manifest$source_sample_count) || anyNA(manifest$selected_sample_count) ||
    any(manifest$source_sample_count < manifest$selected_sample_count)) {
  fail("Staged manifest contains invalid source/selected sample counts")
}
expected_selected <- fifelse(manifest$truth_source == "array", 18L, 7L)
if (any(manifest$selected_sample_count != expected_selected)) {
  fail("Staged manifest selected_sample_count must be 18 for Array and 7 for WGS")
}

eval_mask_manifest <- read_tabular(opts$eval_mask_manifest, "Evaluation mask manifest")
eval_mask_required <- c(
  "logical_run_id", "treatment", "panel", "truth_source", "mask_input_type", "mask_input_path",
  "source_sha256", "position_count", "positions_sha256"
)
eval_mask_missing <- setdiff(eval_mask_required, names(eval_mask_manifest))
if (length(eval_mask_missing)) {
  fail("Evaluation mask manifest is missing: ", paste(eval_mask_missing, collapse = ", "))
}
eval_mask_manifest <- eval_mask_manifest[, ..eval_mask_required]
for (column in eval_mask_required) {
  set(eval_mask_manifest, j = column, value = trimws(as.character(eval_mask_manifest[[column]])))
}
if (anyNA(eval_mask_manifest) || any(eval_mask_manifest == "")) {
  fail("Evaluation mask manifest contains blank values")
}
if (nrow(eval_mask_manifest) != 12L) fail("Expected six Array and six WGS mask sources")
expected_mask_types <- fifelse(
  eval_mask_manifest$truth_source == "array", "concordance_parquet", "per_variant_metrics"
)
if (any(eval_mask_manifest$mask_input_type != expected_mask_types)) {
  fail("Expected Array concordance_parquet and WGS per_variant_metrics mask sources")
}
if (any(!file.exists(eval_mask_manifest$mask_input_path))) {
  fail("One or more original evaluation mask artifacts are missing")
}
eval_mask_manifest[, position_count := suppressWarnings(as.integer(position_count))]
if (anyNA(eval_mask_manifest$position_count) || any(eval_mask_manifest$position_count <= 0L)) {
  fail("Evaluation mask manifest contains an invalid position_count")
}
eval_mask_hash_columns <- c("source_sha256", "positions_sha256")
if (any(!vapply(
  unlist(eval_mask_manifest[, ..eval_mask_hash_columns], use.names = FALSE),
  function(value) grepl("^[0-9a-fA-F]{64}$", value),
  logical(1)
))) fail("Evaluation mask manifest contains an invalid SHA-256 value")
mask_condition_key <- function(x) paste(x$logical_run_id, condition_key(x), sep = "\r")
if (anyDuplicated(mask_condition_key(eval_mask_manifest)) ||
    !setequal(mask_condition_key(eval_mask_manifest), mask_condition_key(manifest))) {
  fail("Evaluation mask sources must match all twelve treatment-panel-truth-source runs")
}

eval_mask <- fread(
  opts$eval_mask,
  sep = "\t",
  header = FALSE,
  quote = "",
  na.strings = character(),
  showProgress = FALSE
)
if (ncol(eval_mask) != 2L || !nrow(eval_mask)) {
  fail("Evaluation mask must contain two non-empty columns: CHROM and POS")
}
setnames(eval_mask, c("chrom", "pos"))
# Coordinates come directly from the six Array and six WGS evaluation Parquets.
eval_mask[, position_id := paste(chrom, pos, sep = ":")]
if (any(eval_mask_manifest$position_count < nrow(eval_mask))) {
  fail("A source position count is smaller than the twelve-run mask intersection")
}
eval_mask_sha256 <- digest::digest(file = opts$eval_mask, algo = "sha256")
eval_mask_manifest_sha256 <- digest::digest(file = opts$eval_mask_manifest, algo = "sha256")

truth_manifest <- read_tabular(opts$truth_manifest, "Truth manifest")
truth_manifest_required <- c(
  "truth_source", "vcf", "source_path", "source_sample_count", "selected_sample_count",
  "source_samples_sha256", "source_signature_sha256", "candidate_vcf_sha256",
  "candidate_header_sha256"
)
truth_manifest_missing <- setdiff(truth_manifest_required, names(truth_manifest))
if (length(truth_manifest_missing)) {
  fail("Truth manifest is missing: ", paste(truth_manifest_missing, collapse = ", "))
}
truth_manifest <- truth_manifest[, ..truth_manifest_required]
for (column in truth_manifest_required) {
  set(truth_manifest, j = column, value = trimws(as.character(truth_manifest[[column]])))
}
if (anyNA(truth_manifest) || any(truth_manifest == "")) fail("Truth manifest contains blank values")
if (!setequal(truth_manifest$truth_source, c("array", "wgs")) || nrow(truth_manifest) != 2L) {
  fail("Truth manifest must contain exactly one array and one WGS source")
}
if (anyDuplicated(truth_manifest$truth_source)) fail("Truth manifest contains duplicate truth sources")
if (any(!file.exists(truth_manifest$vcf))) fail("One or more staged truth VCFs are missing")
truth_manifest[, `:=`(
  source_sample_count = suppressWarnings(as.integer(source_sample_count)),
  selected_sample_count = suppressWarnings(as.integer(selected_sample_count))
)]
if (anyNA(truth_manifest$source_sample_count) || anyNA(truth_manifest$selected_sample_count) ||
    any(truth_manifest$source_sample_count < truth_manifest$selected_sample_count)) {
  fail("Truth manifest contains invalid source/selected sample counts")
}
expected_truth_selected <- fifelse(truth_manifest$truth_source == "array", 18L, 7L)
if (any(truth_manifest$selected_sample_count != expected_truth_selected)) {
  fail("Truth manifest selected_sample_count must be 18 for Array and 7 for WGS")
}
if (any(!vapply(
  unlist(truth_manifest[, ..hash_columns], use.names = FALSE),
  function(value) grepl("^[0-9a-fA-F]{64}$", value),
  logical(1)
))) fail("Truth manifest contains an invalid SHA-256 value")
truth_paths <- setNames(truth_manifest$vcf, truth_manifest$truth_source)
if (!identical(normalizePath(opts$array_truth), normalizePath(truth_paths[["array"]])) ||
    !identical(normalizePath(opts$wgs_truth), normalizePath(truth_paths[["wgs"]]))) {
  fail("--array-truth/--wgs-truth must match the candidate VCFs in --truth-manifest")
}

sample_map <- read_tabular(opts$sample_map, "Sample map")
sample_required <- c("sample_id", "truth_source", "array_group")
sample_missing <- setdiff(sample_required, names(sample_map))
if (length(sample_missing)) fail("Sample map is missing: ", paste(sample_missing, collapse = ", "))
sample_map <- sample_map[, ..sample_required]
for (column in sample_required) set(sample_map, j = column, value = trimws(as.character(sample_map[[column]])))
if (anyNA(sample_map) || any(sample_map == "")) fail("Sample map contains blank values")
if (anyDuplicated(sample_map$sample_id)) fail("Sample map contains duplicate sample IDs")
if (any(!sample_map$truth_source %in% c("array", "wgs"))) fail("Unexpected sample truth_source")
if (nrow(sample_map[truth_source == "array"]) != 18L || nrow(sample_map[truth_source == "wgs"]) != 7L) {
  fail("Sample map must contain 18 array and 7 WGS samples")
}
if (any(!sample_map[truth_source == "array", array_group] %in% c("Group I", "Group II"))) {
  fail("Array samples must belong to Group I or Group II")
}
if (any(sample_map[truth_source == "wgs", array_group] != "WGS")) {
  fail("WGS sample array_group must be WGS")
}

run_bcftools <- function(arguments, stdout = FALSE) {
  result <- system2("bcftools", arguments, stdout = stdout, stderr = "")
  status <- attr(result, "status")
  if (!is.null(status) && status != 0L) {
    fail("bcftools failed: bcftools ", paste(arguments, collapse = " "))
  }
  result
}

read_vcf <- function(path, label) {
  header <- run_bcftools(c("view", "-h", shQuote(path)), stdout = TRUE)
  header_line <- grep("^#CHROM\\t", header, value = TRUE)
  if (length(header_line) != 1L) fail(label, " does not contain one #CHROM header")
  header_columns <- strsplit(sub("^#", "", header_line), "\t", fixed = FALSE)[[1]]
  if (length(header_columns) < 10L) fail(label, " does not contain sample columns")
  samples <- header_columns[-seq_len(9L)]

  body_path <- tempfile(fileext = ".vcf")
  on.exit(unlink(body_path), add = TRUE)
  run_bcftools(c("view", "-H", "-o", shQuote(body_path), shQuote(path)))
  if (!file.exists(body_path) || file.info(body_path)$size == 0) {
    fail(label, " contains no candidate variants")
  }
  records <- fread(
    body_path,
    sep = "\t",
    header = FALSE,
    quote = "",
    na.strings = character(),
    showProgress = FALSE
  )
  expected_columns <- 9L + length(samples)
  if (ncol(records) != expected_columns) {
    fail(label, " has ", ncol(records), " body columns; expected ", expected_columns)
  }
  setnames(
    records,
    c(
      "chrom", "pos", "vcf_id", "ref", "alt", "qual_raw", "filter",
      "info_raw", "format_keys", samples
    )
  )
  records[, qual_raw := as.character(qual_raw)]
  records[, pos := suppressWarnings(as.integer(pos))]
  if (anyNA(records$pos)) fail(label, " contains a non-integer position")
  records[, qual := suppressWarnings(as.numeric(fifelse(
    qual_raw == ".", NA_character_, qual_raw
  )))]
  records[, `:=`(
    position_id = paste(chrom, pos, sep = ":"),
    locus_id = paste(chrom, pos, ref, alt, sep = ":")
  )]

  definition_rows <- function(kind) {
    prefix <- paste0("^##", kind, "=<ID=")
    lines <- header[grepl(prefix, header)]
    if (!length(lines)) return(data.table(kind = character(), id = character(), number = character(), type = character()))
    parsed <- lapply(lines, function(line) {
      match <- regexec(
        paste0("^##", kind, "=<ID=([^,>]+),Number=([^,>]+),Type=([^,>]+)"),
        line
      )
      values <- regmatches(line, match)[[1]]
      if (length(values) != 4L) return(NULL)
      data.table(kind = kind, id = values[[2]], number = values[[3]], type = values[[4]])
    })
    rbindlist(parsed, fill = TRUE)
  }

  list(
    path = path,
    label = label,
    header = header,
    samples = samples,
    records = records,
    definitions = rbindlist(list(definition_rows("INFO"), definition_rows("FORMAT")))
  )
}

sanitize_field_id <- function(id) {
  result <- tolower(gsub("[^A-Za-z0-9]+", "_", id))
  result <- gsub("^_+|_+$", "", result)
  if (any(!nzchar(result))) fail("Cannot sanitize one or more VCF field IDs")
  result
}

validate_field_ids <- function(definitions, kind) {
  subset <- unique(definitions[definitions$kind == kind, .(id)])
  if (!nrow(subset)) return(invisible(NULL))
  subset[, column := sanitize_field_id(id)]
  collisions <- subset[, .N, by = column][N > 1L]
  if (nrow(collisions)) {
    fail(kind, " field IDs collide after column-name normalization: ", paste(collisions$column, collapse = ", "))
  }
}

unique_position_records <- function(vcf, label) {
  records <- vcf$records
  simple_snp <- nchar(records$ref) == 1L & nchar(records$alt) == 1L &
    !grepl(",", records$alt, fixed = TRUE) &
    records$ref %chin% c("A", "C", "G", "T") & records$alt %chin% c("A", "C", "G", "T")
  records <- records[simple_snp]
  position_counts <- records[, .N, by = position_id]
  unique_positions <- position_counts[N == 1L, position_id]
  excluded <- nrow(position_counts[N != 1L])
  records <- records[position_id %chin% unique_positions]
  if (!nrow(records)) fail(label, " has no unique biallelic SNP positions")
  list(records = records, duplicate_position_count = excluded)
}

cat("Reading staged candidate VCFs\n")
run_vcfs <- setNames(vector("list", nrow(manifest)), manifest$source_id)
duplicate_counts <- data.table(source_id = character(), duplicate_position_count = integer())
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i]
  object <- read_vcf(row$vcf, paste0("source ", row$source_id))
  expected_samples <- sample_map[truth_source == row$truth_source, sample_id]
  if (!setequal(object$samples, expected_samples)) {
    fail("Source ", row$source_id, " does not contain exactly its expected target sample set")
  }
  object$records <- object$records[, c(
    "chrom", "pos", "vcf_id", "ref", "alt", "qual", "qual_raw", "filter",
    "info_raw", "format_keys", "position_id", "locus_id", expected_samples
  ), with = FALSE]
  unique_result <- unique_position_records(object, paste0("source ", row$source_id))
  object$records <- unique_result$records
  run_vcfs[[row$source_id]] <- object
  duplicate_counts <- rbind(
    duplicate_counts,
    data.table(source_id = row$source_id, duplicate_position_count = unique_result$duplicate_position_count)
  )
}
all_run_definitions <- rbindlist(lapply(run_vcfs, `[[`, "definitions"), fill = TRUE)
validate_field_ids(all_run_definitions, "INFO")
validate_field_ids(all_run_definitions, "FORMAT")

array_truth <- read_vcf(opts$array_truth, "array truth")
wgs_truth <- read_vcf(opts$wgs_truth, "WGS truth")
array_samples <- sample_map[truth_source == "array", sample_id]
wgs_samples <- sample_map[truth_source == "wgs", sample_id]
if (!setequal(array_truth$samples, array_samples)) fail("Array truth sample set does not match sample map")
if (!setequal(wgs_truth$samples, wgs_samples)) fail("WGS truth sample set does not match sample map")
array_truth$records <- array_truth$records[, c(
  "chrom", "pos", "vcf_id", "ref", "alt", "qual", "qual_raw", "filter",
  "info_raw", "format_keys", "position_id", "locus_id", array_samples
), with = FALSE]
wgs_truth$records <- wgs_truth$records[, c(
  "chrom", "pos", "vcf_id", "ref", "alt", "qual", "qual_raw", "filter",
  "info_raw", "format_keys", "position_id", "locus_id", wgs_samples
), with = FALSE]
array_unique <- unique_position_records(array_truth, "array truth")
wgs_unique <- unique_position_records(wgs_truth, "WGS truth")
array_truth$records <- array_unique$records
wgs_truth$records <- wgs_unique$records

missing_array_truth_positions <- setdiff(
  eval_mask$position_id,
  array_truth$records$position_id
)
if (length(missing_array_truth_positions)) {
  fail(
    "Staged Array truth is missing ", length(missing_array_truth_positions),
    " twelve-run evaluator-mask positions"
  )
}
for (source_id in manifest[truth_source == "array", source_id]) {
  missing_source_positions <- setdiff(
    eval_mask$position_id,
    run_vcfs[[source_id]]$records$position_id
  )
  if (length(missing_source_positions)) {
    fail(
      "Array source ", source_id, " is missing ", length(missing_source_positions),
      " positions retained by all twelve completed evaluations"
    )
  }
}

cat("Constructing the exact-allele candidate mask within the twelve-run evaluator mask\n")
exact_key_sets <- lapply(run_vcfs, function(x) x$records$locus_id)
exact_key_sets[["wgs_truth"]] <- wgs_truth$records$locus_id
common_keys <- Reduce(intersect, exact_key_sets)
array_positions <- array_truth$records$position_id
eligible_positions <- intersect(eval_mask$position_id, array_positions)
common_keys <- common_keys[
  sub("^([^:]+:[0-9]+):.*$", "\\1", common_keys) %chin% eligible_positions
]
if (!length(common_keys)) fail("No loci remain in the cross-run/truth intersection")

common_mask <- unique(rbindlist(lapply(run_vcfs, function(x) {
  x$records[locus_id %chin% common_keys, .(locus_id, position_id, chrom, pos, ref, alt)]
})))[locus_id %chin% common_keys]
if (anyDuplicated(common_mask$locus_id) || nrow(common_mask) != length(common_keys)) {
  fail("Common-mask key construction is inconsistent")
}
common_mask[, chromosome_number := suppressWarnings(as.integer(sub("^Chr", "", chrom)))]
setorder(common_mask, chromosome_number, pos, ref, alt)
common_mask[, chromosome_number := NULL]
common_keys <- common_mask$locus_id
exact_allele_locus_count <- nrow(common_mask)
cat("Exact-allele candidate loci before sample completeness: ", exact_allele_locus_count, "\n", sep = "")

extract_format_tag <- function(format_keys, sample_values, tag) {
  result <- rep(NA_character_, length(sample_values))
  schemas <- unique(format_keys)
  schemas <- schemas[!is.na(schemas) & schemas != "" & schemas != "."]
  for (schema in schemas) {
    keys <- strsplit(schema, ":", fixed = TRUE)[[1]]
    tag_position <- match(tag, keys)
    if (is.na(tag_position)) next
    rows <- which(format_keys == schema)
    pieces <- tstrsplit(sample_values[rows], ":", fixed = TRUE, fill = NA_character_)
    if (tag_position <= length(pieces)) result[rows] <- pieces[[tag_position]]
  }
  result[result %in% c("", ".")] <- NA_character_
  result
}

extract_info_tag <- function(info_values, tag) {
  result <- rep(NA_character_, length(info_values))
  pieces <- strsplit(info_values, ";", fixed = TRUE)
  prefix <- paste0(tag, "=")
  for (i in seq_along(pieces)) {
    current <- pieces[[i]]
    exact <- which(current == tag)
    if (length(exact)) {
      result[[i]] <- "TRUE"
      next
    }
    matched <- which(startsWith(current, prefix))
    if (length(matched)) result[[i]] <- substring(current[[matched[[1]]]], nchar(prefix) + 1L)
  }
  result[result %in% c("", ".")] <- NA_character_
  result
}

add_info_columns <- function(records, definitions) {
  info_defs <- unique(definitions[kind == "INFO", .(id, number, type)])
  observed_tokens <- unlist(strsplit(
    records[!is.na(info_raw) & info_raw != ".", info_raw],
    ";",
    fixed = TRUE
  ))
  observed_ids <- unique(sub("=.*$", "", observed_tokens))
  observed_ids <- observed_ids[nzchar(observed_ids)]
  missing_defs <- setdiff(observed_ids, info_defs$id)
  if (length(missing_defs)) {
    info_defs <- rbind(
      info_defs,
      data.table(id = missing_defs, number = ".", type = "String")
    )
  }
  validate_field_ids(data.table(kind = "INFO", id = info_defs$id), "INFO")
  for (i in seq_len(nrow(info_defs))) {
    definition <- info_defs[i]
    column <- paste0("info__", sanitize_field_id(definition$id))
    value <- extract_info_tag(records$info_raw, definition$id)
    if ((definition$number == "1" || tolower(definition$id) %in% c("info_score", "hwe")) &&
        definition$type %in% c("Integer", "Float")) {
      numeric_value <- suppressWarnings(as.numeric(value))
      invalid <- !is.na(value) & is.na(numeric_value)
      if (any(invalid)) fail("Non-numeric value in numeric INFO/", definition$id)
      value <- numeric_value
    } else if (definition$type == "Flag") {
      value <- !is.na(value)
    }
    set(records, j = column, value = value)
  }
  records
}

add_format_columns <- function(calls, definitions) {
  declared <- unique(definitions[kind == "FORMAT", .(id, number, type)])
  observed <- unique(unlist(strsplit(
    calls[!is.na(format_keys) & format_keys != ".", format_keys],
    ":",
    fixed = TRUE
  )))
  observed <- observed[nzchar(observed)]
  missing_defs <- setdiff(observed, declared$id)
  if (length(missing_defs)) {
    declared <- rbind(
      declared,
      data.table(id = missing_defs, number = ".", type = "String")
    )
  }
  validate_field_ids(data.table(kind = "FORMAT", id = declared$id), "FORMAT")
  declared <- declared[toupper(id) != "HD"]
  for (i in seq_len(nrow(declared))) {
    definition <- declared[i]
    column <- paste0("fmt__", sanitize_field_id(definition$id))
    value <- extract_format_tag(calls$format_keys, calls$sample_format_raw, definition$id)
    if (definition$number == "1" && definition$type %in% c("Integer", "Float")) {
      numeric_value <- suppressWarnings(as.numeric(value))
      invalid <- !is.na(value) & is.na(numeric_value)
      if (any(invalid)) fail("Non-numeric value in numeric FORMAT/", definition$id)
      value <- numeric_value
    }
    set(calls, j = column, value = value)
  }
  calls
}

gt_alt_dosage <- function(gt) {
  normalized <- gsub("|", "/", gt, fixed = TRUE)
  result <- rep(NA_integer_, length(normalized))
  result[normalized == "0/0"] <- 0L
  result[normalized %in% c("0/1", "1/0")] <- 1L
  result[normalized == "1/1"] <- 2L
  result
}

nucleotide_group <- function(base) {
  fifelse(base %chin% c("A", "T"), "A", fifelse(base %chin% c("C", "G"), "B", NA_character_))
}

# Check the wide VCF records before expanding all INFO/FORMAT fields. Do not
# use correctness, GP, DS, INFO_SCORE, HWE, or evaluator n_pairs for selection.
# These rules mirror the gt_valid/truth_valid checks on the final call table.
complete_genotype_keys <- function(vcf, require_wgs_quality = FALSE) {
  records <- vcf$records[position_id %chin% common_mask$position_id]
  complete <- rep(TRUE, nrow(records))
  for (sample_id in vcf$samples) {
    sample_values <- records[[sample_id]]
    valid <- !is.na(gt_alt_dosage(extract_format_tag(records$format_keys, sample_values, "GT")))
    if (require_wgs_quality) {
      gq <- suppressWarnings(as.numeric(extract_format_tag(records$format_keys, sample_values, "GQ")))
      dp <- suppressWarnings(as.numeric(extract_format_tag(records$format_keys, sample_values, "DP")))
      valid <- valid & !is.na(gq) & !is.na(dp) & gq >= opts$min_gq & dp >= opts$min_dp
    }
    complete <- complete & valid
  }
  records[complete, .(locus_id, position_id)]
}

cat("Requiring valid GT and aligned truth for all 25 samples across all six conditions (150 comparisons per locus)\n")
complete_keys <- common_keys
for (source_id in names(run_vcfs)) {
  complete_keys <- intersect(complete_keys, complete_genotype_keys(run_vcfs[[source_id]])$locus_id)
}
complete_keys <- intersect(
  complete_keys,
  complete_genotype_keys(wgs_truth, require_wgs_quality = TRUE)$locus_id
)
complete_array_positions <- complete_genotype_keys(array_truth)$position_id
# Array A/B-to-ALT alignment is unresolved for A/T and C/G allele pairs.
common_mask <- common_mask[
  locus_id %chin% complete_keys & position_id %chin% complete_array_positions &
    nucleotide_group(ref) != nucleotide_group(alt)
]
completeness_excluded_locus_count <- exact_allele_locus_count - nrow(common_mask)
cat("Loci excluded by all-samples/all-runs completeness: ", completeness_excluded_locus_count, "\n", sep = "")
cat("Final all-samples/all-runs mask loci (150/150 valid comparisons): ", nrow(common_mask), "\n", sep = "")
if (!nrow(common_mask)) fail("No loci have valid GT and aligned truth for every sample in every run")
common_keys <- common_mask$locus_id
mask_payload <- paste0(paste(common_keys, collapse = "\n"), "\n")
mask_sha256 <- digest::digest(mask_payload, algo = "sha256", serialize = FALSE)
fwrite(common_mask, opts$mask_out, sep = "\t", quote = FALSE)

make_truth_long <- function(vcf, truth_source, expected_samples) {
  records <- vcf$records
  if (truth_source == "wgs") {
    records <- records[locus_id %chin% common_keys]
  } else {
    records <- records[position_id %chin% common_mask$position_id]
  }
  truth <- melt(
    records,
    id.vars = c(
      "chrom", "pos", "vcf_id", "ref", "alt", "qual", "qual_raw", "filter",
      "info_raw", "format_keys", "position_id", "locus_id"
    ),
    measure.vars = expected_samples,
    variable.name = "sample_id",
    value.name = "truth_sample_format_raw",
    variable.factor = FALSE
  )
  setnames(
    truth,
    c("vcf_id", "ref", "alt", "qual", "qual_raw", "filter", "info_raw", "format_keys", "locus_id"),
    c(
      "truth_vcf_id", "truth_ref", "truth_alt", "truth_qual", "truth_qual_raw", "truth_filter",
      "truth_info_raw", "truth_format_keys", "truth_locus_id"
    )
  )
  truth[, truth_gt_raw := extract_format_tag(truth_format_keys, truth_sample_format_raw, "GT")]
  truth[, truth_gq := suppressWarnings(as.numeric(extract_format_tag(
    truth_format_keys, truth_sample_format_raw, "GQ"
  )))]
  truth[, truth_dp := suppressWarnings(as.numeric(extract_format_tag(
    truth_format_keys, truth_sample_format_raw, "DP"
  )))]
  truth
}

array_truth_long <- make_truth_long(array_truth, "array", array_samples)
wgs_truth_long <- make_truth_long(wgs_truth, "wgs", wgs_samples)

cat("Extracting locus INFO and sample FORMAT fields (excluding FORMAT/HD)\n")
call_tables <- vector("list", nrow(manifest))
all_headers <- list()
all_definitions <- list()
for (i in seq_len(nrow(manifest))) {
  config <- manifest[i]
  object <- run_vcfs[[config$source_id]]
  records <- object$records[locus_id %chin% common_keys]
  records <- add_info_columns(records, object$definitions)
  expected_samples <- sample_map[truth_source == config$truth_source, sample_id]
  calls <- melt(
    records,
    id.vars = setdiff(names(records), expected_samples),
    measure.vars = expected_samples,
    variable.name = "sample_id",
    value.name = "sample_format_raw",
    variable.factor = FALSE
  )
  calls <- add_format_columns(calls, object$definitions)
  calls[, `:=`(
    source_id = config$source_id,
    treatment = config$treatment,
    panel = config$panel,
    truth_source = config$truth_source
  )]
  calls <- merge(
    calls,
    sample_map,
    by = c("sample_id", "truth_source"),
    all.x = TRUE,
    sort = FALSE
  )
  if (anyNA(calls$array_group)) fail("Source ", config$source_id, " contains an unmapped sample")

  if (config$truth_source == "array") {
    calls <- merge(
      calls,
      array_truth_long,
      by = c("chrom", "pos", "position_id", "sample_id"),
      all.x = TRUE,
      sort = FALSE
    )
  } else {
    calls <- merge(
      calls,
      wgs_truth_long,
      by.x = c("chrom", "pos", "position_id", "locus_id", "sample_id"),
      by.y = c("chrom", "pos", "position_id", "truth_locus_id", "sample_id"),
      all.x = TRUE,
      sort = FALSE
    )
  }

  all_headers[[config$source_id]] <- object$header
  all_definitions[[config$source_id]] <- object$definitions
  call_tables[[i]] <- calls
}
calls <- rbindlist(call_tables, fill = TRUE, use.names = TRUE)
rm(call_tables, run_vcfs)
invisible(gc())
emitted_info_columns <- sort(grep("^info__", names(calls), value = TRUE))
emitted_format_columns <- sort(grep("^fmt__", names(calls), value = TRUE))

if (nrow(calls) != nrow(common_mask) * 6L * 25L) {
  fail(
    "Call-table row count is ", nrow(calls), "; expected ",
    nrow(common_mask) * 6L * 25L
  )
}

if (!"fmt__gt" %in% names(calls)) fail("QUILT2 candidate VCFs do not expose FORMAT/GT")
if (!"fmt__gp" %in% names(calls)) calls[, fmt__gp := NA_character_]
if (!"info__info_score" %in% names(calls)) calls[, info__info_score := NA_real_]
if (!"info__hwe" %in% names(calls)) calls[, info__hwe := NA_real_]
calls[, info__info_score := suppressWarnings(as.numeric(info__info_score))]
calls[, info__hwe := suppressWarnings(as.numeric(info__hwe))]
calls[, info_score_valid := is.finite(info__info_score) & info__info_score >= 0 & info__info_score <= 1]
calls[, hwe_valid := is.finite(info__hwe) & info__hwe >= 0 & info__hwe <= 1]

cat("Analysis compromise: ", analysis_compromise, "\n", sep = "")
cat("Calculating aligned dosages and reported GT correctness\n")
calls[, imputed_gt_dosage := gt_alt_dosage(fmt__gt)]
calls[, gt_valid := !is.na(imputed_gt_dosage)]
calls[, `:=`(
  ref_group = nucleotide_group(ref),
  alt_group = nucleotide_group(alt)
)]
calls[, array_strand_ambiguous := truth_source == "array" &
  (is.na(ref_group) | is.na(alt_group) | ref_group == alt_group)]

truth_index_dosage <- gt_alt_dosage(calls$truth_gt_raw)
calls[, truth_dosage_aligned := as.integer(NA)]
wgs_rows <- which(calls$truth_source == "wgs")
calls[wgs_rows, truth_dosage_aligned := truth_index_dosage[wgs_rows]]
array_rows <- which(calls$truth_source == "array")
array_b_dosage <- truth_index_dosage[array_rows]
array_alt_group <- calls$alt_group[array_rows]
array_aligned <- fifelse(
  array_alt_group == "B",
  array_b_dosage,
  fifelse(array_alt_group == "A", 2L - array_b_dosage, NA_integer_)
)
calls[array_rows, truth_dosage_aligned := array_aligned]

calls[, allele_alignment_status := fifelse(
  truth_source == "wgs",
  "exact_ref_alt",
  fifelse(array_strand_ambiguous, "array_strand_ambiguous", "array_ab_to_alt")
)]
calls[, truth_exclusion_reason := NA_character_]
calls[is.na(truth_gt_raw), truth_exclusion_reason := "missing_truth_gt"]
calls[!is.na(truth_gt_raw) & is.na(truth_dosage_aligned), truth_exclusion_reason := "invalid_truth_gt"]
calls[truth_source == "array" & array_strand_ambiguous, truth_exclusion_reason := "array_strand_ambiguous"]
calls[
  truth_source == "wgs" & is.na(truth_exclusion_reason) & (is.na(truth_gq) | is.na(truth_dp)),
  truth_exclusion_reason := "missing_truth_gq_or_dp"
]
calls[
  truth_source == "wgs" & is.na(truth_exclusion_reason) & truth_gq < opts$min_gq,
  truth_exclusion_reason := "truth_gq_below_threshold"
]
calls[
  truth_source == "wgs" & is.na(truth_exclusion_reason) & truth_dp < opts$min_dp,
  truth_exclusion_reason := "truth_dp_below_threshold"
]
calls[, truth_valid := is.na(truth_exclusion_reason) & !is.na(truth_dosage_aligned)]
# GT correctness depends only on the reported GT and valid aligned truth, never GP.
calls[, gt_comparable := truth_valid & gt_valid]
calls[, gt_correct := fifelse(
  gt_comparable,
  imputed_gt_dosage == truth_dosage_aligned,
  NA
)]
complete_call_counts <- calls[, .(
  n_calls = .N,
  n_valid = sum(gt_comparable)
), by = locus_id]
if (anyNA(complete_call_counts) ||
    any(complete_call_counts$n_calls != 150L | complete_call_counts$n_valid != 150L) ||
    anyDuplicated(calls, by = c("locus_id", "treatment", "panel", "sample_id"))) {
  fail("Final mask must contain exactly 150 distinct valid GT-truth comparisons per locus")
}

# Parse GP temporarily to calculate gp_valid and maxGP. Only the original vector,
# validity flag and maximum are exported; these never select or score GT calls.
# Tied probabilities retain their shared maximum.
gp_pieces <- tstrsplit(calls$fmt__gp, ",", fixed = TRUE, fill = NA_character_)
gp_component_count <- fifelse(
  is.na(calls$fmt__gp) | calls$fmt__gp == "",
  0L,
  nchar(calls$fmt__gp) - nchar(gsub(",", "", calls$fmt__gp, fixed = TRUE)) + 1L
)
for (index in 1:3) {
  value <- if (length(gp_pieces) >= index) suppressWarnings(as.numeric(gp_pieces[[index]])) else rep(NA_real_, nrow(calls))
  set(calls, j = paste0("fmt__gp_", index - 1L), value = value)
}
gp_matrix <- as.matrix(calls[, .(fmt__gp_0, fmt__gp_1, fmt__gp_2)])
gp_finite <- gp_component_count == 3L & rowSums(is.finite(gp_matrix)) == 3L
gp_in_range <- gp_finite & rowSums(gp_matrix >= 0 & gp_matrix <= 1) == 3L
gp_sum <- rowSums(gp_matrix, na.rm = FALSE)
gp_valid <- gp_in_range & abs(gp_sum - 1) <= opts$gp_sum_tolerance
calls[, `:=`(
  gp_valid = gp_valid,
  maxGP = fifelse(gp_valid, pmax(fmt__gp_0, fmt__gp_1, fmt__gp_2), NA_real_)
)]

cat("Calculating fixed-denominator site metrics\n")
source_expected <- sample_map[, .N, by = truth_source]
setnames(source_expected, "N", "n_expected_source")
source_metrics <- calls[, .(
  n_valid_source = sum(!is.na(gt_correct)),
  n_correct_source = sum(gt_correct %in% TRUE)
), by = .(locus_id, treatment, panel, truth_source)]
source_metrics <- merge(source_metrics, source_expected, by = "truth_source", all.x = TRUE)
source_metrics[, percent_correct_source := fifelse(
  n_valid_source == n_expected_source,
  100 * n_correct_source / n_expected_source,
  NA_real_
)]
calls <- merge(
  calls,
  source_metrics,
  by = c("locus_id", "treatment", "panel", "truth_source"),
  all.x = TRUE,
  sort = FALSE
)

panel_metrics <- calls[, .(
  n_valid_panel = sum(!is.na(gt_correct)),
  n_correct_panel = sum(gt_correct %in% TRUE)
), by = .(locus_id, treatment, panel)]
panel_metrics[, `:=`(
  n_expected_panel = 25L,
  site_cc_eligible = n_valid_panel == 25L
)]
panel_metrics[, percent_correct_panel := fifelse(
  site_cc_eligible,
  100 * n_correct_panel / n_expected_panel,
  NA_real_
)]
calls <- merge(calls, panel_metrics, by = c("locus_id", "treatment", "panel"), all.x = TRUE, sort = FALSE)

pooled_metrics <- calls[, .(
  n_valid_pooled = sum(!is.na(gt_correct)),
  n_correct_pooled = sum(gt_correct %in% TRUE)
), by = .(locus_id, treatment)]
pooled_metrics[, `:=`(
  n_expected_pooled = 75L,
  pooled_site_eligible = n_valid_pooled == 75L
)]
pooled_metrics[, percent_correct_pooled := fifelse(
  pooled_site_eligible,
  100 * n_correct_pooled / n_expected_pooled,
  NA_real_
)]

site_info <- unique(calls[, .(
  locus_id, treatment, panel, truth_source, info__info_score, info_score_valid
)])
if (nrow(site_info) != nrow(common_mask) * 2L * 3L * 2L) {
  fail("Expected one INFO_SCORE value per locus x treatment x panel x truth source")
}
site_info <- merge(site_info, source_expected, by = "truth_source", all.x = TRUE)
weighted_info <- site_info[, .(
  info_score_source_count = sum(info_score_valid),
  mean_info_score = if (all(info_score_valid)) {
    sum(info__info_score * n_expected_source) / sum(n_expected_source)
  } else {
    NA_real_
  }
), by = .(locus_id, treatment)]
pooled_metrics <- merge(pooled_metrics, weighted_info, by = c("locus_id", "treatment"), all.x = TRUE)
pooled_metrics[, info_plot_eligible := pooled_site_eligible &
  info_score_source_count == 6L & !is.na(mean_info_score)]
calls <- merge(calls, pooled_metrics, by = c("locus_id", "treatment"), all.x = TRUE, sort = FALSE)

truth_consistency <- calls[, .(
  valid_values = list(unique(truth_dosage_aligned[truth_valid]))
), by = .(locus_id, sample_id)]
if (truth_consistency[, any(lengths(valid_values) > 1L)]) {
  fail("Aligned truth dosage is inconsistent across source runs")
}
truth_once <- calls[, .(
  truth_valid_once = any(truth_valid),
  truth_dosage_once = {
    values <- unique(truth_dosage_aligned[truth_valid])
    if (length(values)) values[[1]] else NA_integer_
  }
), by = .(locus_id, sample_id)]
truth_frequency <- truth_once[, .(
  truth_called_samples = sum(truth_valid_once),
  truth_alt_allele_count = sum(truth_dosage_once[truth_valid_once])
), by = locus_id]
truth_frequency[, truth_maf_eligible := truth_called_samples == 25L]
truth_frequency[, truth_maf := fifelse(
  truth_maf_eligible,
  pmin(truth_alt_allele_count / 50, 1 - truth_alt_allele_count / 50),
  NA_real_
)]
calls <- merge(calls, truth_frequency, by = "locus_id", all.x = TRUE, sort = FALSE)

calls[, chromosome_number := suppressWarnings(as.integer(sub("^Chr", "", chrom)))]
setorder(calls, chromosome_number, pos, treatment, panel, truth_source, sample_id)
calls[, c("chromosome_number", "ref_group", "alt_group") := NULL]

# Keep fixed VCF fields in metadata when constant within every physical run.
# A field that varies within any run remains a locus column, as agreed in the audit.
constant_vcf_fields <- list()
for (column in c("vcf_id", "qual", "filter", "truth_filter")) {
  values <- unique(calls[, c("source_id", column), with = FALSE])
  if (all(values[, .N, by = source_id]$N == 1L)) {
    constant_vcf_fields[[column]] <- values
  } else {
    locus_columns <- c(locus_columns, column)
  }
}

# INFO columns lose only the info__ prefix; preserve their existing payload types.
exported_info_columns <- sub("^info__", "", emitted_info_columns)
if (length(intersect(exported_info_columns, setdiff(names(calls), emitted_info_columns)))) {
  fail("An INFO tag conflicts with an existing column after removing the info__ prefix")
}
setnames(calls, emitted_info_columns, exported_info_columns)
# Missing optional fields keep a stable base schema and never change the mask.
for (column in setdiff(c("eaf", "erc", "eac", "paf", "fmt__gp"), names(calls))) {
  set(calls, j = column, value = NA_character_)
}
for (column in setdiff(c("info_score", "hwe", "fmt__ds"), names(calls))) {
  set(calls, j = column, value = NA_real_)
}
locus_columns <- c(locus_columns, setdiff(exported_info_columns, locus_columns))
call_columns <- c(call_columns, setdiff(emitted_format_columns, call_columns))
loci <- unique(calls[, ..locus_columns])
sample_calls <- calls[, ..call_columns]
if (nrow(loci) != nrow(common_mask) * 12L || anyDuplicated(loci, by = join_columns)) {
  fail("Locus table must contain one row per locus x physical source (12 rows per locus)")
}
if (anyDuplicated(sample_calls, by = c(join_columns, "sample_id"))) {
  fail("Call table contains duplicate locus x source x sample rows")
}

field_dictionary <- rbindlist(c(
  lapply(seq_len(nrow(manifest)), function(i) {
    definitions <- all_definitions[[manifest$source_id[[i]]]]
    if (!nrow(definitions)) return(NULL)
    definitions[, source_id := manifest$source_id[[i]]]
    definitions
  }),
  list(
    transform(array_truth$definitions, source_id = "array_truth"),
    transform(wgs_truth$definitions, source_id = "wgs_truth")
  )
), fill = TRUE)

metadata <- list(
  extraction_schema_version = schema_version,
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  accuracy_target = "reported_gt",
  analysis_compromise = analysis_compromise,
  eval_mask_position_count = as.character(nrow(eval_mask)),
  eval_mask_sha256 = eval_mask_sha256,
  eval_mask_manifest_sha256 = eval_mask_manifest_sha256,
  eval_mask_definition = paste(
    "CHROM/POS presence intersection of six Array concordance Parquets and six WGS",
    "per_variant_metrics Parquet datasets; no call-validity or correctness selection.",
    "The same position mask is applied to all Array and WGS sources and truth."
  ),
  eval_mask_source_hash_definition = paste(
    "Array source_sha256 hashes the Parquet file bytes. WGS source_sha256 hashes",
    "a TSV of SHA-256 and absolute path for every *.parquet part beneath the dataset,",
    "ordered by path using LC_ALL=C; each TSV row ends with a newline.",
    "positions_sha256 hashes the selected-chromosome coordinate TSV."
  ),
  common_locus_count = as.character(nrow(common_mask)),
  exact_allele_locus_count = as.character(exact_allele_locus_count),
  completeness_excluded_locus_count = as.character(completeness_excluded_locus_count),
  mask_completeness = "all_samples_all_runs",
  mask_expected_comparisons_per_locus = "150",
  final_mask_definition = paste(
    "Unique allele-compatible SNPs within the twelve-run position intersection,",
    "with valid reported GT in every physical source and valid aligned truth for all",
    "18 Array and seven WGS samples: 150 comparisons across six conditions per locus.",
    "WGS truth GQ/DP thresholds and unambiguous Array A/B alignment are required.",
    "Correctness, GP, DS, INFO_SCORE, and HWE do not determine mask membership."
  ),
  final_mask_attrition_from_eval = as.character(nrow(eval_mask) - nrow(common_mask)),
  expected_call_rows = as.character(nrow(common_mask) * 6L * 25L),
  expected_locus_rows = as.character(nrow(common_mask) * 12L),
  join_columns_json = jsonlite::toJSON(join_columns),
  fixed_denominators_json = jsonlite::toJSON(list(
    source = list(array = 18L, wgs = 7L), panel = 25L, pooled = 75L,
    truth_samples = 25L, truth_alleles = 50L
  ), auto_unbox = TRUE),
  constant_vcf_fields_json = jsonlite::toJSON(constant_vcf_fields, dataframe = "rows", auto_unbox = TRUE, na = "null"),
  output_schema_definition = paste(
    "Locus_mock and Calls_mock base order, with info_score, info_score_valid, hwe, eaf,",
    "erc, eac, paf and hwe_valid appended to calls;",
    "additional emitted INFO/FORMAT tags retained",
    "at their own resolution except FORMAT/HD. INFO prefix removed, FORMAT prefix retained.",
    "vcf_id, qual, filter and truth_filter are metadata if constant within every source,",
    "otherwise locus columns. Raw containers, GP components, duplicate IDs and constant",
    "completeness flags are not exported as columns. Source VCFs remain unchanged."
  ),
  call_site_annotations_definition = paste(
    "Original INFO_SCORE, HWE, EAF, ERC, EAC, PAF and INFO_SCORE/HWE validity flags",
    "copied from the same locus x physical source",
    "to every selected sample row, and also retained in the locus table. Shared within",
    "locus x treatment x panel x truth_source, not necessarily across Array and WGS.",
    "No averaging, allele reorientation, or sample-specific recalculation is applied;",
    "payload types and missing values match the locus table. Missing/invalid annotations",
    "do not exclude otherwise evaluable GT calls. HWE zero is retained as a valid p-value."
  ),
  mask_sha256 = mask_sha256,
  wgs_truth_min_gq = as.character(opts$min_gq),
  wgs_truth_min_dp = as.character(opts$min_dp),
  gp_sum_tolerance = as.character(opts$gp_sum_tolerance),
  gp_valid_definition = paste(
    "Three finite probabilities in [0,1] with sum within gp_sum_tolerance of one;",
    "ties permitted. Audit flag only, never used to filter GT accuracy."
  ),
  maxGP_definition = paste(
    "Maximum of the three original GP probabilities when gp_valid is TRUE; otherwise null.",
    "Ties retain their shared maximum. Not an argmax genotype index or a GT accuracy filter."
  ),
  info_score_aggregation = paste(
    "Target-count-weighted mean of the Array-run and WGS-run INFO_SCORE values",
    "within each panel, then equal mean across Liao, NCBI, and Combined"
  ),
  source_manifest_json = jsonlite::toJSON(manifest, dataframe = "rows", auto_unbox = TRUE),
  eval_mask_manifest_json = jsonlite::toJSON(
    eval_mask_manifest,
    dataframe = "rows",
    auto_unbox = TRUE
  ),
  truth_manifest_json = jsonlite::toJSON(truth_manifest, dataframe = "rows", auto_unbox = TRUE),
  sample_map_json = jsonlite::toJSON(sample_map, dataframe = "rows", auto_unbox = TRUE),
  field_dictionary_json = jsonlite::toJSON(field_dictionary, dataframe = "rows", auto_unbox = TRUE),
  exported_info_columns_json = jsonlite::toJSON(exported_info_columns),
  exported_format_columns_json = jsonlite::toJSON(emitted_format_columns),
  omitted_format_tags_json = jsonlite::toJSON("HD"),
  source_headers_json = jsonlite::toJSON(c(
    all_headers,
    list(array_truth = array_truth$header, wgs_truth = wgs_truth$header)
  ), auto_unbox = TRUE),
  duplicate_position_counts_json = jsonlite::toJSON(rbind(
    duplicate_counts,
    data.table(
      source_id = c("array_truth", "wgs_truth"),
      duplicate_position_count = c(
        array_unique$duplicate_position_count,
        wgs_unique$duplicate_position_count
      )
    )
  ), dataframe = "rows", auto_unbox = TRUE)
)

metadata$extraction_id <- digest::digest(metadata, algo = "sha256")
write_output <- function(data, path, role, keys) {
  table_metadata <- c(metadata, list(
    table_role = role,
    primary_key_json = jsonlite::toJSON(keys),
    output_columns_json = jsonlite::toJSON(names(data))
  ))
  output_table <- Table$create(as.data.frame(data))$ReplaceSchemaMetadata(table_metadata)
  write_parquet(output_table, path, compression = "zstd", use_dictionary = TRUE, write_statistics = TRUE)
  reader <- ParquetFileReader$create(path)
  schema <- reader$GetSchema()
  if (reader$num_rows != nrow(data) || !identical(names(schema), names(data)) ||
      !identical(schema$metadata$extraction_id, metadata$extraction_id) ||
      !identical(schema$metadata$table_role, role)) {
    fail("Written ", role, " Parquet failed row-count/schema/metadata verification")
  }
  cat("Created ", role, " Parquet: ", nrow(data), " rows x ", ncol(data), " columns\n", sep = "")
}
cat("Writing separate locus and sample-call Parquets\n")
write_output(loci, opts$locus_output, "loci", join_columns)
write_output(sample_calls, opts$calls_output, "calls", c(join_columns, "sample_id"))

summary <- data.table(
  metric = c(
    "schema_version", "eval_mask_positions", "eval_mask_sha256",
    "exact_allele_loci", "completeness_excluded_loci", "common_loci",
    "mask_completeness", "mask_expected_comparisons_per_locus",
    "final_mask_attrition_from_eval", "logical_conditions",
    "physical_sources", "samples", "array_samples", "wgs_samples", "call_rows", "mask_sha256",
    "locus_rows", "locus_columns", "call_columns", "extraction_id",
    "emitted_info_columns", "emitted_format_columns", "site_cc_eligible_rows",
    "pooled_site_eligible_rows", "gt_comparable_calls", "accuracy_target", "analysis_compromise"
  ),
  value = c(
    schema_version,
    nrow(eval_mask),
    eval_mask_sha256,
    exact_allele_locus_count,
    completeness_excluded_locus_count,
    nrow(common_mask),
    "all_samples_all_runs",
    150L,
    nrow(eval_mask) - nrow(common_mask),
    6L,
    nrow(manifest),
    nrow(sample_map),
    nrow(sample_map[truth_source == "array"]),
    nrow(sample_map[truth_source == "wgs"]),
    nrow(calls),
    mask_sha256,
    nrow(loci),
    ncol(loci),
    ncol(sample_calls),
    metadata$extraction_id,
    length(emitted_info_columns),
    length(emitted_format_columns),
    unique(calls[, .(locus_id, treatment, panel, site_cc_eligible)])[site_cc_eligible == TRUE, .N],
    unique(calls[, .(locus_id, treatment, pooled_site_eligible)])[pooled_site_eligible == TRUE, .N],
    calls[gt_comparable == TRUE, .N],
    "reported_gt",
    analysis_compromise
  )
)
fwrite(summary, opts$summary_out, sep = "\t", quote = FALSE)

cat("Created: ", opts$locus_output, "\n", sep = "")
cat("Created: ", opts$calls_output, "\n", sep = "")
cat("Twelve-run Array/WGS evaluator mask positions: ", nrow(eval_mask), "\n", sep = "")
cat("Final all-samples/all-runs mask loci (150/150 valid comparisons): ", nrow(common_mask), "\n", sep = "")
cat("Call rows: ", nrow(calls), "\n", sep = "")
