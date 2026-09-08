#!/usr/bin/env Rscript
# Analysis compromise: downstream evaluation scores reported GT regardless of GP.
# This starting presence mask never conditions on GP or call correctness.
# The Parquet builder applies all-samples/all-runs GT/truth completeness afterward.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

fail <- function(...) stop(paste0(...), call. = FALSE)

parse_args <- function(args) {
  opts <- list(input = NULL, output = NULL, chromosomes = NULL, validate_only = FALSE)
  key_map <- c(
    "--input" = "input",
    "--output" = "output",
    "--chr" = "chromosomes"
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      cat(paste(
        "Usage: extract_array_evaluation_positions.R",
        "--input concordance.parquet|per_variant_metrics --chr Chr01,...",
        "[--output positions.tsv] [--validate-only]\n"
      ))
      quit(status = 0)
    }
    if (key == "--validate-only") {
      opts$validate_only <- TRUE
      i <- i + 1L
      next
    }
    if (!key %in% names(key_map) || i == length(args)) {
      fail("Unknown or incomplete option: ", key)
    }
    opts[[key_map[[key]]]] <- args[[i + 1L]]
    i <- i + 2L
  }

  required <- c("input", "chromosomes")
  missing <- required[vapply(required, function(name) is.null(opts[[name]]), logical(1))]
  if (length(missing)) fail("Missing required option(s): ", paste(missing, collapse = ", "))
  if (!opts$validate_only && is.null(opts$output)) fail("--output is required unless --validate-only is used")
  opts
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))
if (!file.exists(opts$input)) fail("Parquet input not found: ", opts$input)

selected <- strsplit(opts$chromosomes, ",", fixed = TRUE)[[1]]

# Both evaluators supply standardised CHROM/POS; no ID parsing, allele inference,
# chromosome rewriting or coordinate repair is needed. This helper selects
# coordinates only; the builder checks call validity in the staged VCFs.
# WGS CHROM lives in Hive directory names (CHROM=ChrNN), not in each part file.
# Arrow projects only coordinates and filters chromosomes before collecting.
positions <- open_dataset(opts$input, format = "parquet", partitioning = hive_partition()) |>
  dplyr::select(CHROM, POS) |>
  dplyr::filter(CHROM %in% selected) |>
  dplyr::collect() |>
  as.data.table()
setnames(positions, c("CHROM", "POS"), c("chrom", "pos"))
if (!nrow(positions)) fail("No selected chromosome positions found in ", opts$input)
setorder(positions, chrom, pos)

if (opts$validate_only) {
  cat("Read ", nrow(positions), " selected positions in ", opts$input, "\n", sep = "")
  quit(status = 0)
}
if (file.exists(opts$output)) fail("Refusing to overwrite existing output: ", opts$output)
fwrite(positions, opts$output, sep = "\t", col.names = FALSE, quote = FALSE)
cat("Extracted ", nrow(positions), " selected positions\n", sep = "")
