
setwd(Sys.getenv("MITOEVO_SCMTMPM_ROOT", unset = "/PATH/TO/scmtMPM_output"))

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(ggprism)
  library(forcats)
  library(rstatix)
  library(ggpubr)
  library(patchwork)
  library(stringr)
  library(Biostrings)
  
})


species_list <- c("human", "mouse", "rat", "chicken", "cobra", "sturgeon", "stingray")

low.coverage.position <- readRDS("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/low.coverage.position.7SpeciesTable.RDS")

meta_path <- c(
  "human"   ="/PATH/TO/pbmc.mt/SNP/human/clean_bam/merge_bam/human.PBMConly.sample.barcode.celltype.mtPercent.withArchR.csv",
  "mouse"   ="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/mouse.meta.csv",
  "rat"     ="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/rat.meta.csv",
  "chicken" ="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/chicken.meta.csv",
  "cobra"   ="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/cobra.meta.csv",
  "sturgeon"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/sturgeon.meta.csv",
  "stingray"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/stingray.meta.csv"
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

ref_fasta <- list(
  "human" = "/PATH/TO/pbmc.mt/SNP/human_mix/mgatk/downstream/rm.germline/human.germline.rm.MT.genome.fa",
  "mouse" = "/PATH/TO/pbmc.mt/ref/mouse/new/mm39_mito.fa",
  "rat" = "/PATH/TO/pbmc.mt/SNP/rat/mgatk/downstream/rm.germline/rat.germline.rm.MT.genome.fa",
  "chicken" = "/PATH/TO/pbmc.mt/SNP/chicken/mgatk/downstream/rm.germline/chicken.germline.rm.MT.genome.fa",
  "cobra" = "/PATH/TO/pbmc.mt/SNP/cobra/mgatk/downstream/rm.germline/cobra.germline.rm.MT.genome.fa",
  "sturgeon" = "/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk/downstream/rm.germline/sturgeon.germline.rm.MT.genome.fa",
  "stingray" = "/PATH/TO/pbmc.mt/SNP/stingray/mgatk/downstream/rm.germline/stingray.germline.rm.MT.genome.fa"
)


variant.list <- readRDS("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Final_Variant_List_7species_260128.RDS")

species.colors <- c("#625D9E", "#CCC9E6", "#DCC1DD", "#F2E2B2", "#C5DEBA", "#AAD4F8", "#5184B2")
names(species.colors) <- c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray")

mt.celltype.colors <- c('#B3659F','#8E90B9', '#CCE0F5','#3361A5','#7FB0B4')
names(mt.celltype.colors) <- c("B cell","T cell","NK cell","T/NK cell","Monocyte")

get_mtDNA_length_from_fasta <- function(fasta_file) {
  fasta_seq <- readDNAStringSet(fasta_file)
  
  seq_length <- width(fasta_seq)[1]
  
  message("mtDNA length detected from FASTA: ", seq_length)
  return(seq_length)
}

coerce_meta <- function(meta_raw, sp) {
  meta <- as.data.frame(meta_raw)
  
  if (!("barcode" %in% colnames(meta))) {
    candidate_cols <- c("barcode","cell","Cell","cellID","cell_id","raw_barcode","barcodes","Barcode","X")
    found <- intersect(candidate_cols, colnames(meta))
    if (length(found) > 0) {
      cand <- found[1]
      if (cand == "X" && all(grepl("^[0-9]+$", meta[[cand]]))) {
      } else {
        meta <- meta %>% dplyr::rename(barcode = all_of(cand))
      }
    }
  }
  
  if (!("barcode" %in% colnames(meta))) {
    stop(
      "Meta file for '", sp, "' does not contain per-cell barcodes. ",
      "You likely pointed to an aggregated mt.percent table (Species/Cell_Type/mt.percent). ",
      "Please use a PER-CELL metadata file (e.g., seurat_meta/*.meta.csv) that contains barcodes."
    )
  }
  
  if (!("main_celltype2" %in% colnames(meta))) {
    candidate_ct <- c("main_celltype2","Cell_Type","celltype","cell_type","celltype0","celltype1","main_celltype","MainCelltype2")
    found_ct <- intersect(candidate_ct, colnames(meta))
    if (length(found_ct) == 0) {
      stop("Meta file for '", sp, "' lacks a cell-type column. Need main_celltype2 (or an equivalent).")
    }
    meta <- meta %>% dplyr::rename(main_celltype2 = all_of(found_ct[1]))
  }
  
  if (!("mt.percent" %in% colnames(meta))) {
    candidate_mt <- c("mt.percent","mt_percent","percent.mt","percent_mt","MTpercent","mtPercent")
    found_mt <- intersect(candidate_mt, colnames(meta))
    if (length(found_mt) > 0) {
      meta <- meta %>% dplyr::rename(mt.percent = all_of(found_mt[1]))
    } else {
      meta$mt.percent <- NA_real_
    }
  }
  
  if (!("mt_reads" %in% colnames(meta))) meta$mt_reads <- NA_real_
  if (!("total_reads" %in% colnames(meta))) meta$total_reads <- NA_real_
  
  meta %>% dplyr::select(barcode, main_celltype2, mt_reads, total_reads, mt.percent)
}

read_mtCov_cell_from_long <- function(cov_gz_path, drop_positions = NULL) {
  cov_dt <- data.table::fread(cov_gz_path, header = FALSE)
  data.table::setnames(cov_dt, 1:3, c("position", "barcode", "coverage"))
  
  if (is.null(drop_positions) || length(drop_positions) == 0) {
    cov_dt[, .(mtCov_all = sum(coverage, na.rm = TRUE),
               mtCov_eff = sum(coverage, na.rm = TRUE)), by = barcode]
  } else {
    cov_dt[, .(
      mtCov_all = sum(coverage, na.rm = TRUE),
      mtCov_eff = sum(coverage[!position %in% drop_positions], na.rm = TRUE)
    ), by = barcode]
  }
}


calc_mtMut_counts_cell <- function(dat) {
  cn <- colnames(dat)
  if ("alt_total" %in% cn) {
    return(dat %>% 
             dplyr::group_by(cell) %>% 
             dplyr::summarise(mtMut_counts_cell = sum(alt_total, na.rm = TRUE), .groups = "drop"))
  }
  stop("Cannot compute mtMut_counts_cell: need alt_depth/alt_total or (allele_frequency or vaf) + total_depth.")
}

get_comparisons <- function(species) {
  if (species %in% c("cobra","sturgeon","stingray")) {
    return(list(
      c("T/NK cell", "Monocyte"),
      c("B cell", "Monocyte"),
      c("B cell", "T/NK cell")
    ))
  }
  list(
    c("NK cell", "Monocyte"),
    c("T cell", "Monocyte"),
    c("B cell", "Monocyte"),
    c("B cell", "T cell"),
    c("B cell", "NK cell"),
    c("T cell", "NK cell")
  )
}

compute_species_stats <- function(data, y_col = "scmtMPM") {
  data <- data %>% dplyr::rename(Species = species, Cell_Type = main_celltype2)
  sps <- unique(data$Species)
  
  stat_list <- lapply(sps, function(sp) {
    sp_data <- dplyr::filter(data, Species == sp)
    comps <- get_comparisons(sp)
    
    valid_comps <- Filter(function(pair) all(pair %in% unique(sp_data$Cell_Type)), comps)
    if (length(valid_comps) == 0) return(NULL)
    
    fml <- as.formula(paste0(y_col, " ~ Cell_Type"))
    stat <- sp_data %>% wilcox_test(formula = fml, comparisons = valid_comps) %>%
      add_significance("p") %>% dplyr::mutate(Species = sp)
    
    y_base <- max(sp_data[[y_col]], na.rm = TRUE)
    stat %>% arrange(group1, group2) %>%
      dplyr::mutate(y.position = y_base * (1.05 + 0.07 * row_number()))
  })
  
  bind_rows(stat_list)
}



meta.list <- list()

result.list <- map(species_list, function(sp) {
  message("Processing ", sp, " ...")
  
  drop_pos <- low.coverage.position[[sp]]
  
  meta_raw <- read.csv(meta_path[[sp]], header = TRUE, check.names = FALSE)
  meta <- coerce_meta(meta_raw, sp)
  meta.list[[sp]] <- meta_raw
  
  valid_vars <- variant.list[[sp]]
  
  path <- paste0(
    "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/",
    sp, "_detail_filtered_df_260128.tsv"
  )
  dat <- read_tsv(path, show_col_types = FALSE) %>%
    dplyr::rename(cell = barcode) %>%  
    dplyr::rename(variant = mutation) %>% 
    dplyr::filter(alt_bias > 0.45, alt_bias < 0.55) %>%
    dplyr::filter(variant %in% valid_vars)%>%
    dplyr::filter(!pos %in% drop_pos)
  
  cell_mut_n <- dat %>% 
    dplyr::group_by(cell) %>% 
    dplyr::summarise(n_mut = dplyr::n(), .groups = "drop")
  
  cell_mut_counts <- calc_mtMut_counts_cell(dat) # cell, mtMut_counts_cell
  
  cov_path <- coverage_files[[sp]]
  if (is.null(cov_path) || !file.exists(cov_path)) stop("Coverage file not found for ", sp, ": ", cov_path)
  mtCov_dt <- read_mtCov_cell_from_long(cov_path, drop_positions = drop_pos) # barcode, mtCov_cell
  
  mt_length <- get_mtDNA_length_from_fasta(ref_fasta[[sp]])
  
  df <- meta %>%
    left_join(cell_mut_n, by = c("barcode" = "cell")) %>%
    left_join(cell_mut_counts, by = c("barcode" = "cell")) %>%
    left_join(as.data.frame(mtCov_dt), by = "barcode") %>%
    dplyr::mutate(
      n_mut = replace_na(n_mut, 0),
      mtMut_counts_cell = replace_na(mtMut_counts_cell, 0),
      mtCov_all = replace_na(mtCov_all, 0),
      mtCov_eff = replace_na(mtCov_eff, 0)
    ) %>%
    dplyr::filter(mtCov_all > 10 * mt_length) %>%
    dplyr::mutate(
      mtCov_cell = mtCov_eff,
      scmtMPM = dplyr::if_else(mtCov_cell > 0, (mtMut_counts_cell / mtCov_cell) * 1e6, NA_real_)
    )
  
  
  df$species <- sp
  df
})


mutload_all <- bind_rows(result.list) %>%
  dplyr::filter(main_celltype2 %in% c("Monocyte","B cell","T cell","T/NK cell","NK cell")) %>%
  dplyr::filter(!is.na(scmtMPM), is.finite(scmtMPM), scmtMPM >= 0)

mutload_all$main_celltype2 <- factor(mutload_all$main_celltype2, levels = c("Monocyte","B cell","T cell","T/NK cell","NK cell"))
mutload_all$species <- factor(mutload_all$species, levels = species_list)



mutload_summary <- mutload_all %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    avg_mut_per_cell = mean(n_mut),
    avg_mt_percent = mean(mt.percent, na.rm = TRUE),
    avg_scmtMPM = mean(scmtMPM, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(mutload_all, "PerCell_scmtMPM_260128.tsv")
write_tsv(mutload_summary, "PerSpecies_scmtMPM_summary_260128.tsv")

celltype_levels <- c("Monocyte","B cell","T cell","T/NK cell","NK cell")

plot_data0 <- mutload_all %>%
  dplyr::filter(!is.na(scmtMPM), is.finite(scmtMPM), scmtMPM >= 0) %>%
  dplyr::filter(mt.percent >= 1) %>%                
  dplyr::mutate(
    main_celltype2 = factor(main_celltype2, levels = celltype_levels),
    species = factor(species, levels = species_list)
  ) %>%
  dplyr::filter(!is.na(main_celltype2))

plot_data <- plot_data0 %>%
  dplyr::group_by(species, main_celltype2) %>%
  dplyr::mutate(
    n_in_group = dplyr::n(),
    q05 = ifelse(n_in_group > 100, quantile(scmtMPM, 0.05, na.rm = TRUE), -Inf),
    q95 = ifelse(n_in_group > 100, quantile(scmtMPM, 0.95, na.rm = TRUE),  Inf)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(scmtMPM >= q05, scmtMPM <= q95) %>%
  dplyr::select(-q05, -q95)


stat_df <- compute_species_stats(plot_data, y_col = "scmtMPM") %>%
  dplyr::rename(species = Species) %>%
  dplyr::mutate(
    p.signif = dplyr::case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  )


 
  y_by_sp <- plot_data %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(ymax = max(scmtMPM, na.rm = TRUE), .groups = "drop")
  
  stat_df <- stat_df %>%
    dplyr::left_join(y_by_sp, by = "species") %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(y.position = ymax * (1.05 + 0.07 * dplyr::row_number())) %>%
    dplyr::ungroup() %>%
    dplyr::select(-ymax)


stat_df$species <- factor(stat_df$species, levels = species_list)

p1 <- ggplot(plot_data, aes(x = main_celltype2, y = scmtMPM, fill = main_celltype2)) +
  geom_violin(trim = TRUE, color = "grey40", alpha = 0.9) +
  geom_boxplot(width = 0.1, color = "grey25", fill = "white", outlier.shape = NA) +
  facet_wrap(~ species, ncol = 1, strip.position = "top", scales = "free_y", drop = FALSE) +  # ✅ free_y
  scale_fill_manual(values = mt.celltype.colors) +
  ggpubr::stat_pvalue_manual(
    stat_df,
    label = "p.signif",
    xmin = "group1", xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 3,
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(face = "bold", size = 14, color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 10), color = "black"),
    plot.margin = margin(t = 5, r = 25, b = 5, l = 5)
  ) +
  labs(
    x = "Cell type",
    y = "scmtMPM (mutations per million bases)",
    title = "Per-cell mitochondrial mutational burden"
  )


ggsave("./plots/Reorder_violin_percell_scmtMPM_q05q95_freeY_linear_260128.pdf", p1, width = 5, height = 11)




ct_mean <- plot_data %>%
  dplyr::group_by(species, main_celltype2) %>%
  dplyr::summarise(mean_scmtMPM = mean(scmtMPM, na.rm = TRUE), .groups = "drop")

ct_global <- ct_mean %>%
  dplyr::group_by(main_celltype2) %>%
  dplyr::summarise(mean_mt = mean(mean_scmtMPM, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_mt)

ct_global$main_celltype2 <- factor(ct_global$main_celltype2, levels = c("Monocyte","B cell","T cell","T/NK cell","NK cell"))
ct_mean$main_celltype2 <- factor(ct_mean$main_celltype2, levels = levels(ct_global$main_celltype2))
ct_mean <- ct_mean %>%
  dplyr::mutate(species = tools::toTitleCase(as.vector(species))) %>%
  dplyr::mutate(species = factor(species, levels = names(species.colors)))

p2 <- ggplot(ct_global, aes(x = main_celltype2, y = mean_mt)) +
  geom_col(width = 0.6, color = "black",
           fill = mt.celltype.colors[match(ct_global$main_celltype2, names(mt.celltype.colors))]) +
  geom_point(
    data = ct_mean,
    aes(x = main_celltype2, y = mean_scmtMPM, fill = species),
    shape = 21, color = "black", size = 2, stroke = 0.5,
    position = position_jitter(width = 0.2)
  ) +
  geom_text(aes(label = sprintf("%.2f", mean_mt)), vjust = -0.3, size = 4, color = "black") +
  scale_fill_manual(values = species.colors) +
  theme_prism(base_size = 12) +
  labs(x = "Cell type", y = "Average scmtMPM", title = "Average scmtMPM per cell type")

ggsave("./plots/Reorder_top_bar_celltype_mean_scmtMPM_260128.pdf", p2, width = 6, height = 2.5)


sp_means <- plot_data %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(avg_scmtMPM = mean(scmtMPM, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(species = factor(species, levels = rev(species_list)))

ct_by_sp <- plot_data %>%
  dplyr::group_by(species, main_celltype2) %>%
  dplyr::summarise(subgroup_mean = mean(scmtMPM, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(species = factor(species, levels = rev(species_list)))

p3 <- ggplot(sp_means) +
  geom_col(aes(x = species, y = avg_scmtMPM), width = 0.2, fill = "#8291B0", alpha = 1) +
  geom_point(
    data = ct_by_sp,
    aes(x = species, y = subgroup_mean, fill = main_celltype2),
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.4),
    size = 2, shape = 21, color = "grey40"
  ) +
  scale_fill_manual(values = mt.celltype.colors, name = "Cell Type") +
  geom_text(aes(x = species, y = avg_scmtMPM, label = sprintf("%.2f", avg_scmtMPM)),
            hjust = -0.1, size = 3, color = "black") +
  coord_flip() +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt")) +
  labs(x = "Species", y = "scmtMPM", title = "Across species")

ggsave("./plots/right_bar_species_avg_scmtMPM_260128.pdf", p3, width = 5, height = 15)


celltype_levels <- c("Monocyte","B cell","T cell","T/NK cell","NK cell")

ct_agg <- plot_data %>%
  dplyr::group_by(species, main_celltype2) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    sum_alt = sum(mtMut_counts_cell, na.rm = TRUE),
    sum_cov = sum(mtCov_cell, na.rm = TRUE),
    scmtMPM_agg = 1e6 * (sum_alt / sum_cov),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    main_celltype2 = factor(main_celltype2, levels = celltype_levels),
    species = factor(species, levels = species_list),
    species_title = tools::toTitleCase(as.character(species)),
    species_title = factor(species_title, levels = names(species.colors))
  )

ct_norm <- ct_agg %>%
  group_by(species) %>%
  dplyr::mutate(
    species_mean = mean(scmtMPM_agg),
    relative_burden = scmtMPM_agg / species_mean
  ) %>%
  ungroup()


ct_bar <- ct_norm %>%
  dplyr::group_by(main_celltype2) %>%
  dplyr::summarise(
    mean_relative_burden = mean(relative_burden, na.rm = TRUE),
    .groups = "drop"
  )

p_contrib2 <- ggplot() +
  geom_col(
    data = ct_bar,
    aes(x = main_celltype2, y = mean_relative_burden),
    width = 0.6, color = "black",
    fill = mt.celltype.colors[match(ct_bar$main_celltype2, names(mt.celltype.colors))]
  ) +
  geom_point(
    data = ct_norm,
    aes(x = main_celltype2, y = relative_burden, fill = species_title),
    shape = 21, color = "black", size = 2, stroke = 0.4,
    position = position_jitter(width = 0.18, height = 0)
  ) +
  geom_text(
    data = ct_bar,
    aes(x = main_celltype2, y = mean_relative_burden, label = sprintf("%.2f", mean_relative_burden)),
    vjust = -0.35, size = 3.8, color = "black"
  ) +
  scale_fill_manual(values = species.colors, name = "Species") +
  theme_prism(base_size = 12) +
  labs(
    x = "Cell type",
    y = "Relative mtMPM",
    title = "Within-species normalised cell-type aggregated scmtMPM"
  )

ggsave("./plots/Reorder_bar_celltype_scmtMPMagg_dotBySpecies_260128.pdf",
       p_contrib2, width = 6.5, height = 3.0)




sp_agg <- ct_agg %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    sum_alt_sp = sum(sum_alt, na.rm = TRUE),
    sum_cov_sp = sum(sum_cov, na.rm = TRUE),
    scmtMPM_spAgg = 1e6 * (sum_alt_sp / sum_cov_sp),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    species = factor(species, levels = species_list),
    species_title = tools::toTitleCase(as.character(species)),
    species_title = factor(species_title, levels = names(species.colors))
  )

sp_agg <- sp_agg %>%
  dplyr::mutate(species = factor(species, levels = rev(species_list)))

ct_agg2 <- ct_agg %>%
  dplyr::mutate(
    main_celltype2 = factor(main_celltype2, levels = celltype_levels),
    species = factor(species, levels = rev(species_list))
  )

p3_new <- ggplot() +
  geom_col(
    data = sp_agg,
    aes(y = species, x = scmtMPM_spAgg),
    width = 0.3,
    fill = "#8291B0", alpha = 1
  ) +
  geom_point(
    data = ct_agg2,
    aes(y = species, x = scmtMPM_agg, fill = main_celltype2),
    shape = 21, color = "black", size = 2, stroke = 0.35,
    position = position_jitter(height = 0.18, width = 0)
  ) +
  geom_text(data = sp_agg,
            aes(y = species, x = scmtMPM_spAgg, label = sprintf("%.2f", scmtMPM_spAgg)),
            hjust = -0.1, size = 3, color = "black") +
  scale_fill_manual(values = mt.celltype.colors, name = "Cell type") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt")) +
  labs(
    x = "Aggregated mtMPM",
    y = "Species",
    title = "Species-level aggregated mtMPM"
  )

ggsave("./plots/right_bar_species_aggregated_mtMPM_260128.pdf",
       p3_new, width = 5, height = 15)
