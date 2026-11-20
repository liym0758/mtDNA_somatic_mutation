#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)

config <- fread(args[1])
output_root <- args[2]
required_config <- c("species", "variant_stats", "numt_variants")
if (!all(required_config %in% names(config))) {
  stop("Configuration requires: ", paste(required_config, collapse = ", "))
}
if (anyDuplicated(config$species)) stop("Duplicate species in configuration")
if (!"allele_prefix" %in% names(config)) config[, allele_prefix := NA_character_]
if (!"ref_allele_file" %in% names(config)) config[, ref_allele_file := NA_character_]
if (!"cell_variant" %in% names(config)) config[, cell_variant := NA_character_]
if (!"cell_strand_counts" %in% names(config)) config[, cell_strand_counts := NA_character_]
if (!"coverage_summary" %in% names(config)) config[, coverage_summary := NA_character_]
if (!"coverage_matrix" %in% names(config)) config[, coverage_matrix := NA_character_]
if (!"cell_metadata" %in% names(config)) config[, cell_metadata := NA_character_]
if (!"coverage_barcode_mode" %in% names(config)) config[, coverage_barcode_mode := "plain"]
if (!"coverage_samples" %in% names(config)) config[, coverage_samples := NA_character_]

read_bias_from_mgatk <- function(prefix, ref_allele_file) {
  if (is.na(prefix) || !nzchar(prefix)) stop("raw_mean_bias missing and allele_prefix not provided")
  if (is.na(ref_allele_file) || !nzchar(ref_allele_file)) stop("ref_allele_file required for mgatk bias calculation")
  ref <- fread(ref_allele_file, col.names = c("pos", "ref"))
  if (anyDuplicated(ref$pos)) stop("Duplicate positions in reference allele table")
  alleles <- rbindlist(lapply(c("A", "C", "G", "T"), function(base) {
    path <- paste0(prefix, ".", base, ".txt.gz")
    if (!file.exists(path)) stop("Missing mgatk allele matrix: ", path)
    x <- fread(path, header = FALSE, col.names = c("pos", "cell", "fwd", "rev"))
    x[, alt := base]
    x
  }))
  alleles <- merge(alleles, ref, by = "pos", all.x = TRUE)
  if (anyNA(alleles$ref)) stop("mgatk allele position absent from reference allele table")
  alleles[, variant := paste0(pos, ref, ">", alt)]
  alleles[, bias := fwd / (fwd + rev)]
  alleles[, .(
    raw_fwd_total = sum(fwd, na.rm = TRUE),
    raw_rev_total = sum(rev, na.rm = TRUE),
    raw_count_total = sum(fwd + rev, na.rm = TRUE),
    n_cells_detected = uniqueN(cell),
    raw_mean_count_PerCell = mean(fwd + rev, na.rm = TRUE),
    raw_mean_bias = mean(bias, na.rm = TRUE)
  ), by = variant]
}

prepare_cell_variant_details <- function(path) {
  x <- fread(path)
  needed <- c("pos", "cell", "variant", "ref", "alt", "fwd", "rev")
  if (!all(needed %in% names(x))) stop("cell_strand_counts table needs: ", paste(needed, collapse = ", "))
  char_cols <- c("cell", "variant", "ref", "alt")
  x[, (char_cols) := lapply(.SD, function(v) gsub('^"|"$', '', as.character(v))), .SDcols = char_cols]
  x[, `:=`(fwd = as.integer(fwd), rev = as.integer(rev))]
  ref <- x[ref == alt, .(ref_fwd = sum(fwd, na.rm = TRUE),
                          ref_rev = sum(rev, na.rm = TRUE)), by = .(pos, cell, ref)]
  ref[, ref_total := ref_fwd + ref_rev]
  alt <- x[ref != alt, .(alt_fwd = sum(fwd, na.rm = TRUE),
                          alt_rev = sum(rev, na.rm = TRUE)), by = .(pos, cell, variant, ref, alt)]
  alt[, alt_total := alt_fwd + alt_rev]
  alt[, alt_bias := fifelse(alt_total > 0, alt_fwd / alt_total, NA_real_)]
  cells <- merge(alt, ref, by = c("pos", "cell", "ref"), all.x = TRUE)
  cells[is.na(ref_fwd), `:=`(ref_fwd = 0L, ref_rev = 0L, ref_total = 0L)]
  cells[, total_depth := ref_total + alt_total]
  cells[, is_conf_detected_cell := alt_fwd >= 2 & alt_rev >= 2]
  cells
}

add_cell_fractions <- function(stats, cells) {
  if ("n_cells_detected" %in% names(stats)) {
    stats[, pct_conf_detected_cell_recomputed :=
            fifelse(n_cells_detected > 0, n_cells_conf_detected / n_cells_detected, NA_real_)]
    if (!"pct_conf_detected_cell" %in% names(stats)) {
      stats[, pct_conf_detected_cell := pct_conf_detected_cell_recomputed]
    } else {
      disagree <- stats[!is.na(pct_conf_detected_cell_recomputed) &
                        abs(pct_conf_detected_cell - pct_conf_detected_cell_recomputed) > 1e-10, .N]
      if (disagree) warning(disagree, " variants differ between supplied pct_conf_detected_cell and n_conf/n_detected; supplied values retained")
    }
  }
  if (!"pct_conf_detected_cell" %in% names(stats)) {
    stop("pct_conf_detected_cell missing: provide it or n_cells_detected from mgatk allele matrices")
  }
  if (is.null(cells)) return(stats)
  needed <- c("variant", "total_depth", "is_conf_detected_cell")
  if (!all(needed %in% names(cells))) stop("cell_variant table needs: ", paste(needed, collapse = ", "))
  pct_depth <- cells[total_depth >= 10,
                     .(pct_conf_detected_cell_DepthFiltered =
                         sum(is_conf_detected_cell == TRUE, na.rm = TRUE) / .N), by = variant]
  if ("pct_conf_detected_cell_DepthFiltered" %in% names(stats)) {
    stats[, pct_conf_detected_cell_DepthFiltered := NULL]
  }
  stats <- merge(stats, pct_depth, by = "variant", all.x = TRUE, sort = FALSE)
  stats
}

read_position_depth <- function(row) {
  sp <- row$species
  if (!is.na(row$coverage_matrix) && nzchar(row$coverage_matrix)) {
    if (is.na(row$cell_metadata) || !nzchar(row$cell_metadata)) stop(sp, ": cell_metadata required with coverage_matrix")
    meta <- fread(row$cell_metadata)
    if (!"barcode" %in% names(meta)) stop(sp, ": metadata needs barcode")
    if (!is.na(row$coverage_samples) && nzchar(row$coverage_samples)) {
      if (!"sample" %in% names(meta)) stop(sp, ": metadata needs sample")
      meta <- meta[sample %in% strsplit(row$coverage_samples, ";", fixed = TRUE)[[1]]]
    }
    mode <- row$coverage_barcode_mode
    if (is.na(mode) || !nzchar(mode)) mode <- "plain"
    if (mode == "sample_prefix") {
      if (!"sample" %in% names(meta)) stop(sp, ": sample_prefix requires metadata sample")
      barcodes <- paste0(meta$sample, "_", meta$barcode)
    } else if (mode == "plain") {
      barcodes <- meta$barcode
    } else stop(sp, ": invalid coverage_barcode_mode")
    cov <- fread(row$coverage_matrix, header = FALSE,
                 col.names = c("position", "barcode", "depth"))
    cov <- cov[barcode %in% barcodes]
    if (!nrow(cov)) stop(sp, ": no coverage barcodes match metadata")
    cov[, position := as.integer(position)]
    cov <- cov[, .(depth = mean(depth, na.rm = TRUE)), by = position]
    cov[, species_code := sp]
    return(cov)
  }
  if (is.na(row$coverage_summary) || !nzchar(row$coverage_summary)) {
    stop(sp, ": provide coverage_matrix + cell_metadata, or coverage_summary")
  }
  fread(row$coverage_summary)
}

filter_species <- function(row) {
  sp <- row$species
  stats <- fread(row$variant_stats)
  need <- c("variant", "position", "vmr", "n_cells_conf_detected")
  if (!all(need %in% names(stats))) stop(sp, ": variant stats need: ", paste(need, collapse = ", "))
  need_mgatk_summary <- !"raw_mean_bias" %in% names(stats) ||
    (!"pct_conf_detected_cell" %in% names(stats) && !"n_cells_detected" %in% names(stats))
  if (need_mgatk_summary) {
    bias <- read_bias_from_mgatk(row$allele_prefix, row$ref_allele_file)
    missing_fields <- setdiff(names(bias), names(stats))
    stats <- merge(stats, bias[, c("variant", missing_fields), with = FALSE],
                   by = "variant", all.x = TRUE)
  }
  cells <- NULL
  if (!is.na(row$cell_strand_counts) && nzchar(row$cell_strand_counts)) {
    cells <- prepare_cell_variant_details(row$cell_strand_counts)
  } else if (!is.na(row$cell_variant) && nzchar(row$cell_variant)) {
    cells <- fread(row$cell_variant)
  }
  stats <- add_cell_fractions(stats, cells)
  numts <- fread(row$numt_variants)
  if (!"Formatted_Variant" %in% names(numts)) stop(sp, ": NUMT table needs Formatted_Variant")
  cov <- read_position_depth(row)
  if (!all(c("species_code", "position", "depth") %in% names(cov))) {
    stop(sp, ": coverage summary needs species_code, position, depth")
  }
  low_coverage <- unique(cov[species_code == sp & depth < 10, position])
  before <- nrow(stats)
  stats <- stats[!(vmr <= 0.01 & n_cells_conf_detected > 3)]
  stats <- stats[!variant %in% numts$Formatted_Variant]
  stats <- stats[!position %in% low_coverage]
  stats <- stats[n_cells_conf_detected >= 1]
  stats <- stats[raw_mean_bias > 0.45 & raw_mean_bias < 0.55]
  stats <- stats[pct_conf_detected_cell >= 0.03]
  if (sp == "human") stats <- stats[variant != "16278C>T"]
  out_dir <- file.path(output_root, sp)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, "variant_stats_final_filtered.tsv")
  if (file.exists(out_file)) stop("Refusing to overwrite: ", out_file)
  fwrite(stats, out_file, sep = "\t", quote = FALSE)
  message(sp, ": ", before, " -> ", nrow(stats), " variants; ", out_file)
  invisible(out_file)
}

for (i in seq_len(nrow(config))) filter_species(config[i])
