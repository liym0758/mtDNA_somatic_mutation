# ================================
# Region UNIQUE-variant density 
# density = (# unique variants in region) / (effective length of region)
# effective length = positions in region after dropping low_coverage_position
# keep_cells rule: mtCov_cell > 10 * mtDNA_length
# ================================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(stringr)
  library(tibble)
})

AnnoVar_Region_UniqueVarDensity_BarPlt <- function(
    specie,
    gene_pred_file,
    cell_variant_file,
    coverage_file,
    mtDNA_length,
    low_coverage_position = NULL,
    valid_variants = NULL,
    out_dir = NULL,
    output_pdf = NULL
) {
  
  annovar_color <- c(
    `protein-coding` = "#625D9E",
    tRNA = "#C7B1D4",
    rRNA = "#B3659F",
    intergenic = "#D3D3D3"
  )
  
  # output folders
  if (is.null(out_dir)) {
    out_dir <- file.path(
      "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result",
      specie, "result"
    )
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  inter_dir <- file.path(out_dir, "intermediate_uniqueVarDensity_region")
  dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("============== ", specie, " ==============")
  message("mtDNA_length = ", mtDNA_length)
  message("out_dir = ", out_dir)
  message("inter_dir = ", inter_dir)
  
  # --------------------------
  # 1) Build pos_type (position -> region) from genePred
  # --------------------------
  message("Reading genePred and building pos_type (position -> region)...")
  
  gene_pred <- readr::read_tsv(
    gene_pred_file,
    col_names = c(
      "bin","transcript","chr","strand","txStart","txEnd",
      "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds",
      "score","gene","cdsStartStat","cdsEndStat","exonFrames"
    ),
    show_col_types = FALSE
  )
  
  gene_pred_processed <- gene_pred %>%
    dplyr::mutate(
      processed_gene = gsub("^(mt-|MT-)", "", gene, ignore.case = TRUE),
      gene_type = dplyr::case_when(
        cdsStart != cdsEnd ~ "protein_coding",
        grepl("^[Tt]", processed_gene) ~ "tRNA",
        grepl("^[Rr]", processed_gene) ~ "rRNA",
        TRUE ~ "protein_coding"
      )
    )
  
  pos_type <- rep("intergenic", mtDNA_length)
  
  mark_ranges <- function(start_vec, end_vec, label) {
    if (length(start_vec) == 0) return(invisible(NULL))
    df <- data.frame(st = as.integer(start_vec), en = as.integer(end_vec))
    df <- df[is.finite(df$st) & is.finite(df$en), , drop = FALSE]
    if (nrow(df) == 0) return(invisible(NULL))
    
    gr <- GenomicRanges::GRanges(
      seqnames = "MT",
      ranges = IRanges::IRanges(df$st, df$en)
    )
    gr <- GenomicRanges::reduce(gr)
    
    ir <- IRanges::ranges(gr)
    for (k in seq_along(ir)) {
      st <- max(1L, IRanges::start(ir[k]))
      en <- min(mtDNA_length, IRanges::end(ir[k]))
      if (st <= en) pos_type[st:en] <<- label
    }
    invisible(NULL)
  }
  
  # protein-coding uses CDS
  pc <- gene_pred_processed %>% dplyr::filter(gene_type == "protein_coding")
  if (nrow(pc) > 0) mark_ranges(pc$cdsStart + 1L, pc$cdsEnd, "protein-coding")
  
  # tRNA/rRNA use tx ranges
  tr <- gene_pred_processed %>% dplyr::filter(gene_type == "tRNA")
  rr <- gene_pred_processed %>% dplyr::filter(gene_type == "rRNA")
  if (nrow(tr) > 0) mark_ranges(tr$txStart + 1L, tr$txEnd, "tRNA")
  if (nrow(rr) > 0) mark_ranges(rr$txStart + 1L, rr$txEnd, "rRNA")
  
  saveRDS(pos_type, file.path(inter_dir, "pos_type.RDS"))
  
  # --------------------------
  # 2) Read coverage long table, define kept cells (mean cov > 10)
  # --------------------------
  message("Reading long-format coverage table, and selecting retained cells (mean cov > 10)...")
  if (!file.exists(coverage_file)) stop("Coverage file not found: ", coverage_file)
  
  cov_dt <- data.table::fread(coverage_file, header = FALSE)
  data.table::setnames(cov_dt, 1:3, c("position","barcode","coverage"))
  cov_dt[, position := as.integer(position)]
  cov_dt[, coverage := as.numeric(coverage)]
  
  # Drop low coverage positions (blacklist) - same as your old code
  if (!is.null(low_coverage_position) && length(low_coverage_position) > 0) {
    low_coverage_position <- unique(as.integer(low_coverage_position))
    low_coverage_position <- low_coverage_position[is.finite(low_coverage_position)]
    
    n_before <- nrow(cov_dt)
    pos_before <- data.table::uniqueN(cov_dt$position)
    
    cov_dt <- cov_dt[!position %in% low_coverage_position]
    
    n_after <- nrow(cov_dt)
    pos_after <- data.table::uniqueN(cov_dt$position)
    
    msg <- paste0(
      "Drop low_coverage positions: n_pos_drop=", length(low_coverage_position),
      " | unique_pos: ", pos_before, " -> ", pos_after,
      " | n_rows: ", n_before, " -> ", n_after
    )
    message(msg)
    writeLines(msg, file.path(inter_dir, "low_coverage_drop_log.txt"))
    saveRDS(low_coverage_position, file.path(inter_dir, "low_coverage_position_used.RDS"))
  } else {
    message("No low_coverage_position provided; skip dropping positions.")
    low_coverage_position <- integer(0)
  }
  
  mtCov_cell_dt <- cov_dt[, .(mtCov_cell = sum(coverage, na.rm = TRUE)), by = barcode]
  keep_cells <- mtCov_cell_dt[mtCov_cell > 10 * mtDNA_length, barcode]
  
  data.table::fwrite(mtCov_cell_dt, file.path(inter_dir, "mtCov_cell.tsv"), sep = "\t")
  writeLines(keep_cells, file.path(inter_dir, "keep_cells.txt"))
  saveRDS(keep_cells, file.path(inter_dir, "keep_cells.RDS"))
  message("kept cells: ", length(keep_cells))
  
  # --------------------------
  # 3) Effective length per region
  #    = count of positions in that region after dropping low_coverage_position
  # --------------------------
  message("Computing effective length for each region (excluding low-coverage positions)...")
  
  kept_positions <- sort(unique(cov_dt$position))
  kept_positions <- kept_positions[kept_positions >= 1L & kept_positions <= mtDNA_length]
  
  effective_len_df <- tibble::tibble(
    position = kept_positions,
    mutation_type = pos_type[kept_positions]
  ) %>%
    dplyr::count(mutation_type, name = "effective_length") %>%
    dplyr::mutate(
      mutation_type = as.character(mutation_type)
    )
  
  # ensure all 4 regions exist
  all_regions <- c("protein-coding","tRNA","rRNA","intergenic")
  effective_len_df <- tibble::tibble(mutation_type = all_regions) %>%
    dplyr::left_join(effective_len_df, by = "mutation_type") %>%
    dplyr::mutate(effective_length = dplyr::if_else(is.na(effective_length), 0L, effective_length))
  
  readr::write_tsv(effective_len_df, file.path(inter_dir, "region_effective_length.tsv"))
  
  # --------------------------
  # 4) Numerator: unique variants per region (in kept cells; optional valid_variants)
  # --------------------------
  message("Reading cell-variant table, and counting unique variants per region...")
  
  if (!file.exists(cell_variant_file)) stop("Cell-variant file not found: ", cell_variant_file)
  
  dat <- readr::read_tsv(cell_variant_file, show_col_types = FALSE) %>%
    dplyr::rename(cell = barcode, variant = mutation) %>%
    dplyr::mutate(pos = as.integer(pos)) %>%
    dplyr::filter(cell %in% keep_cells)
  
  if (!is.null(valid_variants)) {
    dat <- dat %>% dplyr::filter(variant %in% valid_variants)
  }
  
  # IMPORTANT: also drop low_coverage positions here to match effective length
  if (length(low_coverage_position) > 0) {
    dat <- dat %>% dplyr::filter(!pos %in% low_coverage_position)
  }
  
  dat_dt <- data.table::as.data.table(dat)
  dat_dt[, mutation_type := pos_type[pos]]
  
  uniqueVar_dt <- dat_dt[, .(
    n_unique_variants = data.table::uniqueN(variant),
    n_cells_with_any_var = data.table::uniqueN(cell),
    n_rows_cellVariant = .N
  ), by = mutation_type]
  
  uniqueVar_df <- as.data.frame(uniqueVar_dt) %>%
    dplyr::mutate(mutation_type = as.character(mutation_type))
  
  readr::write_tsv(uniqueVar_df, file.path(inter_dir, "region_unique_variants.tsv"))
  
  # --------------------------
  # 5) Final stats: density = n_unique_variants / effective_length
  # --------------------------
  message("Combining data and computing unique-variant density ...")
  
  region_stats <- tibble::tibble(mutation_type = c("protein-coding","tRNA","rRNA","intergenic")) %>%
    dplyr::left_join(effective_len_df, by = "mutation_type") %>%
    dplyr::left_join(uniqueVar_df, by = "mutation_type") %>%
    dplyr::mutate(
      n_unique_variants = dplyr::if_else(is.na(n_unique_variants), 0L, n_unique_variants),
      n_cells_with_any_var = dplyr::if_else(is.na(n_cells_with_any_var), 0L, n_cells_with_any_var),
      n_rows_cellVariant = dplyr::if_else(is.na(n_rows_cellVariant), 0L, n_rows_cellVariant),
      uniqueVar_density = dplyr::if_else(effective_length > 0,
                                         n_unique_variants / effective_length,
                                         NA_real_),
      mutation_type = factor(
        mutation_type,
        levels = rev(c("protein-coding","tRNA","rRNA","intergenic"))
      )
    )
  
  stats_path <- file.path(out_dir, "Region_AnnoVar_uniqueVarDensity_stats.tsv")
  readr::write_tsv(region_stats, stats_path)
  readr::write_tsv(region_stats, file.path(inter_dir, "Region_AnnoVar_uniqueVarDensity_stats.tsv"))
  
  # --------------------------
  # 6) Plot (same style as your old bar plot)
  # --------------------------
  p <- ggplot(region_stats,
              aes(x = mutation_type, y = uniqueVar_density, fill = mutation_type)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = annovar_color, drop = FALSE) +
    labs(
      title = stringr::str_to_title(specie),
      x = "Region",
      y = "Unique variants / effective bp"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  if (!is.null(output_pdf)) {
    grDevices::pdf(output_pdf, width = 6, height = 3)
    print(p)
    grDevices::dev.off()
  }
  
  return(list(
    stats = region_stats,
    plot = p,
    keep_cells = keep_cells
  ))
}


# ================================
# Batch run ------------------------------
# ================================

base_path <- "/PATH/TO/pbmc.mt/SNP/"
species_list <- c("human","mouse","rat","chicken","cobra","sturgeon","stingray")

gene_pred_file_list <- c(
  "human" = "/PATH/TO/biosoft/annovar_2020/humandb/hg38_chrM_ensGene.txt",
  "mouse" = "/PATH/TO/biosoft/annovar_2020/mm39_MT_db/mm39_MT_ensGene.txt",
  "rat" = "/PATH/TO/biosoft/annovar_2020/rat_MT_db/rn7_MT_ensGene.txt",
  "chicken" = "/PATH/TO/biosoft/annovar_2020/chicken_MT_db/gg1_MT_ensGene.txt",
  "cobra" = "/PATH/TO/biosoft/annovar_2020/cobra_MT_db/na_MT_ensGene.txt",
  "sturgeon" = "/PATH/TO/biosoft/annovar_2020/sturgeon_MT_db/as_CNCB_MT_ensGene.txt",
  "stingray" = "/PATH/TO/biosoft/annovar_2020/stingray_MT_db/pl_MT_ensGene.txt"
)

coverage_files <- c(
  human    = "/PATH/TO/pbmc.mt/SNP/human/mgatk_newCleanBam/merge_pbmc_samples/final/merged.coverage.txt.gz",
  mouse    = "/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/mouse.mgatk.coverage.txt.gz",
  rat      = "/PATH/TO/pbmc.mt/SNP/rat/mgatk_newCleanBam/final/rat.mgatk.coverage.txt.gz",
  chicken  = "/PATH/TO/pbmc.mt/SNP/chicken/mgatk_newCleanBam/final/chicken.mgatk.coverage.txt.gz",
  cobra    = "/PATH/TO/pbmc.mt/SNP/cobra/mgatk_newCleanBam/final/cobra.mgatk.coverage.txt.gz",
  sturgeon = "/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk_newCleanBam/final/sturgeon.mgatk.coverage.txt.gz",
  stingray = "/PATH/TO/pbmc.mt/SNP/stingray/mgatk_newCleanBam/final/stingray.mgatk.coverage.txt.gz"
)

mt_len_list <- c(
  human = 16569,
  mouse = 16299,
  rat = 16313,
  chicken = 16784,
  cobra = 17216,
  sturgeon = 16682,
  stingray = 17500
)

variant.list <- readRDS("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Final_Variant_List_7species_260128.RDS")

low.coverage.position <- readRDS("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/low.coverage.position.7SpeciesTable.RDS")

base_out <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result"

result_unique <- vector("list", length(species_list))
names(result_unique) <- species_list

for (sp in species_list) {
  
  cell_variant_file <- paste0(
    "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/",
    sp, "_detail_filtered_df_260128.tsv"
  )
  
  out_dir <- file.path(base_out, sp, "result")
  pdf_path <- file.path(out_dir, "Region_AnnoVar_uniqueVarDensity_BarPlt.pdf")
  
  result_unique[[sp]] <- AnnoVar_Region_UniqueVarDensity_BarPlt(
    specie = sp,
    gene_pred_file = gene_pred_file_list[[sp]],
    cell_variant_file = cell_variant_file,
    coverage_file = coverage_files[[sp]],
    mtDNA_length = mt_len_list[[sp]],
    valid_variants = variant.list[[sp]],
    low_coverage_position = low.coverage.position[[sp]]$position,
    out_dir = out_dir,
    output_pdf = pdf_path
  )
  
  message("[DONE] ", sp, " -> ", pdf_path)
}

# ================================
# Heatmap: species x region (row z-score) -----------------
# ================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(scales)
})

all_regions <- c("protein-coding","rRNA","tRNA","intergenic")


species_order_top_to_bottom <- c("human","mouse","rat","chicken","cobra","sturgeon","stingray")

red.blue <- c("#A42227","#FFFFFF", "#003E7A")  # high -> low

hm_df <- purrr::map_dfr(species_list, function(sp) {
  stats_path <- file.path(base_out, sp, "result", "Region_AnnoVar_uniqueVarDensity_stats.tsv")
  if (!file.exists(stats_path)) {
    warning("[SKIP] stats not found: ", stats_path)
    return(NULL)
  }
  readr::read_tsv(stats_path, show_col_types = FALSE) %>%
    dplyr::mutate(species = sp) %>%
    dplyr::mutate(mutation_type = as.character(mutation_type))
})

hm_df <- hm_df %>%
  tidyr::complete(species, mutation_type = all_regions, fill = list(uniqueVar_density = 0)) %>%
  group_by(species) %>%
  dplyr::mutate(
    row_mean = mean(uniqueVar_density, na.rm = TRUE),
    row_sd   = sd(uniqueVar_density, na.rm = TRUE),
    z = ifelse(is.na(row_sd) | row_sd == 0, 0, (uniqueVar_density - row_mean) / row_sd)
  ) %>%
  ungroup() %>%
  dplyr::mutate(
    mutation_type = factor(mutation_type, levels = rev(all_regions)),
    species = factor(species, levels = species_order_top_to_bottom)
  )%>%
  dplyr::mutate(
    mutation_type = factor(mutation_type, levels = rev(all_regions)),
    species = factor(
      species,
      levels = species_order_top_to_bottom,
      labels = stringr::str_to_title(species_order_top_to_bottom)  
    )
  )

p_hm <- ggplot(hm_df, aes(x = mutation_type, y = species, fill = z)) +
  geom_tile(color = "black", linewidth = 0.2) +
  coord_flip()+
  scale_fill_gradientn(
    colours = rev(red.blue),          # low -> high
    oob = scales::squish
  ) +
  scale_x_discrete(
    labels = c(
      "protein-coding" = "Protein-\ncoding",
      "rRNA" = "rRNA",
      "tRNA" = "tRNA",
      "intergenic" = "Intergenic"
    )
  )+
  labs( fill = "Row\nz-score",
       title = "Regional mutation density") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, color = "black",size = 12),
    axis.text.y = element_text(color = "black",size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5)
  )

pdf("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/plots/Density_heatmap/Heatmap_uniqueVarDensity_rowZ.pdf",
    width = 8, height = 3)
print(p_hm)
dev.off()


write_tsv(hm_df, file = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/plots/Density_heatmap/Heatmap_uniqueVarDensity_HeatmapData.tsv", quote = "none")
