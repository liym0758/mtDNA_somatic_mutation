library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

plot_VAF_hist_fixedDepth <- function(
    sp, 
    meta_path,
    cov_path, #sp/mgatk_newCleanBam/final/sp.mgatk.coverage.txt.gz
    NUMT_path, #sp/mgatk/NUMT.simulate/results/snv_variants.tsv
    low_coverage_position,
    celltype_colors,
    celltype_resolution = "main_celltype2",
    output_file = NULL, 
    return_data = NULL
) {
  message("Processing species: ", sp)
  
  meta_df <- read.csv(meta_path)
  stopifnot("barcode" %in% colnames(meta_df))
  stopifnot(celltype_resolution %in% colnames(meta_df))
  meta_df$celltype_resolution <- meta_df[[celltype_resolution]]
  
  valid_barcodes <- meta_df$barcode[!is.na(meta_df$celltype_resolution)]
  message("Number of valid barcodes in meta: ", length(valid_barcodes))
  

  NUMT_data <- read.table(NUMT_path, sep = "\t", header = T)
  numt <- NUMT_data$Formatted_Variant
  

  coverage_df <- fread(cov_path, col.names = c("pos", "cell", "depth"))
  
  reads_detail_path <- paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/detailed_data/", sp,
                              "/Cell_Var_ReadsCount_PctConf_Data.tsv")
  reads_detail <- read_tsv(reads_detail_path, show_col_types = FALSE)
  
    
  variant_stats_final_filtered_path <- paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/",
                                    sp, "/variant_stats_final_filtered.tsv")
  variant_stats_final_filtered <- read.table(variant_stats_final_filtered_path,sep = "\t", header = T)
  
  reads_detail <- reads_detail %>%
    dplyr::filter(cell %in% valid_barcodes) %>%
    dplyr::filter(variant %in% variant_stats_final_filtered$variant) %>% 
    dplyr::filter(!pos %in% low_coverage_position) %>% 
    dplyr::filter(!variant %in% numt) %>% 
    left_join(coverage_df, by = c("pos", "cell")) %>%
    dplyr::mutate(total_depth = depth) %>%
    dplyr::mutate(VAF = alt_total / total_depth) %>%
    dplyr::mutate(ref_bias = ref_fwd / ref_total) %>%
    dplyr::select(-depth)
  
  
  
  
  variant_filtered <- reads_detail %>%
    dplyr::filter(alt_bias > 0.45 & alt_bias < 0.55)%>%
    dplyr::group_by(variant) %>%
    dplyr::filter(any(is_conf_detected_cell %in% TRUE, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% 
    pull(variant) %>% unique()
  
  reads_filtered <- reads_detail %>%
    # dplyr::filter(alt_bias > 0.45 & alt_bias < 0.55)%>%
    dplyr::filter(variant %in% variant_filtered)
  
  fwrite(reads_filtered, 
         paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/", sp,
                "/Cell_Var_ReadsCount_PctConf_Data_filtered.tsv"),
         sep = "\t", quote = FALSE, na = "NA")
  cat("  Saved:", paste(sp, "Cell_Var_ReadsCount_PctConf_Data_filtered_260128.tsv"), "\n")
  
  
  
  
  # vaf_distribution <- reads_filtered %>%
  #   summarise(
  #     total_variants = n(),
  #     vaf_1_count = sum(VAF == 1),
  #     vaf_1_prop = vaf_1_count / total_variants,
  #     mean_vaf = mean(VAF),
  #     median_vaf = median(VAF)
  #   )
  # 
  
  
  merged_df <- reads_filtered %>%
    dplyr::rename(barcode = cell, mutation = variant, frequency = VAF) %>%
    left_join(meta_df, by = "barcode") %>%
    dplyr::filter(!is.na(celltype_resolution))
  
  filtered_df <- merged_df %>% dplyr::filter(frequency > 0)
  
  bins <- base::seq(0, 1, by = 0.1)
  bin_labels <- base::sprintf("(%.1f-%.1f]", bins[-length(bins)], bins[-1])
  filtered_df <- filtered_df %>%
    dplyr::mutate(
      bin = base::cut(
        frequency,
        breaks = bins,
        include.lowest = FALSE,
        right = TRUE,
        labels = bin_labels
      )
    )
  
  message("Filtered data rows: ", nrow(filtered_df))
  
  fwrite(filtered_df, 
         paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/", sp, "_detail_filtered_df_260128.tsv"),
         sep = "\t", quote = FALSE, na = "NA")
  cat("  Saved:", paste(sp, "detail_filtered_df.tsv"), "\n")
  
  plot_data <- filtered_df %>%
    group_by(celltype_resolution, bin) %>%
    dplyr::summarise(pair_count = dplyr::n(), .groups = "drop") %>%
    complete(celltype_resolution, bin, fill = list(pair_count = 0)) %>%
    dplyr::mutate(bin = factor(bin, levels = bin_labels))
  
  main_type_order <- c("Monocyte","B cell","T cell","NK cell","T/NK cell")
  sub_type_order <- c(
    "HSP", "B cell","Naive B","Memory B","Plasma B",
    "CD4+ T","CD4+ Naive T","CD4+ Memory T",
    "CD8+ T","CD8+ Naive T","CD8+ Memory T",
    "Delta T","MAIT","NKT","T/NK cell","NK cell",
    "Monocyte","Monocyte_1","Monocyte_2",
    "CD14+ Monocyte","CD16+ Monocyte",
    "cDC","pDC","Neutrophil",
    "Platelet","Erythrocyte","Others"
  )
  
  if (celltype_resolution == "main_celltype2" | celltype_resolution == "main_celltype") {
    celltype_levels <- intersect(main_type_order, unique(plot_data$celltype_resolution))
  } else {
    celltype_levels <- intersect(sub_type_order, unique(plot_data$celltype_resolution))
  }
  plot_data$celltype_resolution <- factor(plot_data$celltype_resolution, levels = celltype_levels)
  
  existing_celltypes <- levels(plot_data$celltype_resolution)
  celltype_colors <- celltype_colors[names(celltype_colors) %in% existing_celltypes]
  plot_data <- plot_data %>% dplyr::filter(celltype_resolution %in% existing_celltypes)
  
  p <- ggplot(plot_data, aes(x = bin, y = pair_count, fill = celltype_resolution)) +
    geom_col(position = "dodge", color = "black", width = 0.8) +
    facet_wrap(~celltype_resolution, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = celltype_colors) +
    labs(x = "Variant Allele Frequency (VAF)", 
         y = "Number of Cell-Mutation Pairs",
         title = str_to_title(sp)) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  if (!is.null(output_file)) {
    n_panels <- length(unique(plot_data$celltype_resolution))
    panel_width <- 3
    fixed_height <- 4
    auto_width <- max(6, n_panels * panel_width)
    ggplot2::ggsave(output_file, p, width = auto_width, height = fixed_height, units = "in")
    message("Saved figure: ", output_file)
  }
  
  if (!is.null(return_data)) {
    detailed_data <- list(
      merged_data = merged_df,
      summary_data = plot_data
    )
    saveRDS(detailed_data, file = return_data)
  }
  
  invisible(p)
}


# ==============================
# ==============================

setwd("/PATH/TO/pbmc.mt/SNP/")

mt.celltype.colors <- c('#B3659F','#8E90B9', '#CCE0F5','#3361A5','#7FB0B4')
names(mt.celltype.colors) <- c("B cell","T cell","NK cell","T/NK cell","Monocyte")

species <- c("human","mouse","rat","chicken","cobra","sturgeon","stingray")

NUMTs_path <- c(
  "human"="/PATH/TO/pbmc.mt/SNP/human_mix/mgatk/NUMT.simulate/results/snv_variants.tsv",
  "mouse"="/PATH/TO/pbmc.mt/SNP/mouse/mgatk/false.positive.test/numt.call/freebayes/snv_variants.tsv",
  "rat"="/PATH/TO/pbmc.mt/SNP/rat/mgatk/NUMT.simulate/results/snv_variants.tsv",
  "chicken"="/PATH/TO/pbmc.mt/SNP/chicken/mgatk/NUMT.simulate/results/snv_variants.tsv",
  "cobra"="/PATH/TO/pbmc.mt/SNP/cobra/mgatk/NUMT.simulate/results/snv_variants.tsv",
  "sturgeon"="/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk/NUMT.simulate/results/snv_variants.tsv",
  "stingray"="/PATH/TO/pbmc.mt/SNP/stingray/mgatk/NUMT.simulate/results/snv_variants.tsv"
)
cov_path <- c(
  "human"="/PATH/TO/pbmc.mt/SNP/human/mgatk_newCleanBam/merge_pbmc_samples/final/merged.coverage.txt.gz",
  "mouse"="/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/mouse.mgatk.coverage.txt.gz",
  "rat"="/PATH/TO/pbmc.mt/SNP/rat/mgatk_newCleanBam/final/rat.mgatk.coverage.txt.gz",
  "chicken"="/PATH/TO/pbmc.mt/SNP/chicken/mgatk_newCleanBam/final/chicken.mgatk.coverage.txt.gz",
  "cobra"="/PATH/TO/pbmc.mt/SNP/cobra/mgatk_newCleanBam/final/cobra.mgatk.coverage.txt.gz",
  "sturgeon"="/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk_newCleanBam/final/sturgeon.mgatk.coverage.txt.gz",
  "stingray"="/PATH/TO/pbmc.mt/SNP/stingray/mgatk_newCleanBam/final/stingray.mgatk.coverage.txt.gz"
)

meta_path <- c(
  "human"="/PATH/TO/pbmc.mt/SNP/human/clean_bam/merge_bam/human.PBMConly.sample.barcode.celltype.mtPercent.withArchR.csv",
  "mouse"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/mouse.meta.csv",
  "rat"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/rat.meta.csv",
  "chicken"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/chicken.meta.csv",
  "cobra"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/cobra.meta.csv",
  "sturgeon"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/sturgeon.meta.csv",
  "stingray"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/stingray.meta.csv"
)

low.coverage.position <- readRDS("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/low.coverage.position.7SpeciesTable.RDS")

for (sp in species) {
  plot_VAF_hist_fixedDepth(
    sp = sp,
    meta_path = meta_path[sp],
    cov_path = cov_path[sp],
    NUMT_path = NUMTs_path[sp],
    low_coverage_position=low.coverage.position[[sp]]$position,
    celltype_colors = mt.celltype.colors,
    celltype_resolution = "main_celltype2",
    output_file = paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Reorder/", sp, "_Reorder_celltype_VAF_hist_260128.pdf"),
    return_data = paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Reorder/", sp, "_Reorder_celltype_VAF_data_260128.RDS")
  )
}

variant.list <- list()
for (i in 1:7) {
  data <- read_tsv(paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/", species[i], "_detail_filtered_df_260128.tsv"))
  variant.list[[i]] <- unique(data$mutation)
  names(variant.list)[i] <- species[i]
}

for(i in 1:7){print(length(variant.list[[i]]))}

saveRDS(variant.list, file = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Final_Variant_List_7species_260128.RDS")




variant_stats_final_filtered <- list()
variant_stats_final_filtered_CovBlacklisted <- list()
for (sp in species){
  variant_stats_final_filtered[[sp]] <- read_tsv(file = paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/",sp,"/variant_stats_final_filtered.tsv"))

  variant_stats_final_filtered_CovBlacklisted[[sp]] <- variant_stats_final_filtered[[sp]] %>% 
    dplyr::filter(variant %in% variant.list[[sp]])
  
  fwrite(variant_stats_final_filtered_CovBlacklisted[[sp]], 
         file = paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/",sp,"/variant_stats_final_filtered_CovBlacklisted.tsv"), 
         sep="\t", quote=FALSE)
}
