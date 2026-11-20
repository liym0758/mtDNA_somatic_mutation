#QC PLOTS

# 01 average depth per cell violin ----------------
# =========================
# Mean MT depth per cell (from mgatk coverage)
# Filter: mean_depth > 10
# Output: violin plot across species
# =========================

# =========================
# Mean MT depth per cell
# Intersect with meta$barcode
# Filter: mean_depth > 10
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# -------------------------
# Input coverage files
# -------------------------
coverage_files <- c(
  human    = "/PATH/TO/pbmc.mt/SNP/human/mgatk_newCleanBam/merge_pbmc_samples/final/merged.coverage.txt.gz",
  mouse    = "/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/mouse.mgatk.coverage.txt.gz",
  rat      = "/PATH/TO/pbmc.mt/SNP/rat/mgatk_newCleanBam/final/rat.mgatk.coverage.txt.gz",
  chicken  = "/PATH/TO/pbmc.mt/SNP/chicken/mgatk_newCleanBam/final/chicken.mgatk.coverage.txt.gz",
  cobra    = "/PATH/TO/pbmc.mt/SNP/cobra/mgatk_newCleanBam/final/cobra.mgatk.coverage.txt.gz",
  sturgeon = "/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk_newCleanBam/final/sturgeon.mgatk.coverage.txt.gz",
  stingray = "/PATH/TO/pbmc.mt/SNP/stingray/mgatk_newCleanBam/final/stingray.mgatk.coverage.txt.gz"
)

# -------------------------
# Meta files
# -------------------------
meta_path <- c(
  human    = "/PATH/TO/pbmc.mt/SNP/human/clean_bam/merge_bam/human.PBMConly.sample.barcode.celltype.mtPercent.withArchR.csv",
  mouse    = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/mouse.meta.csv",
  rat      = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/rat.meta.csv",
  chicken  = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/chicken.meta.csv",
  cobra    = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/cobra.meta.csv",
  sturgeon = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/sturgeon.meta.csv",
  stingray = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/stingray.meta.csv"
)

species_order <- names(coverage_files)

# -------------------------
# Function
# -------------------------
compute_mean_depth_with_meta <- function(sp) {
  
  message("Processing: ", sp)
  
  # ---- read coverage ----
  dt <- data.table::fread(
    coverage_files[[sp]],
    col.names = c("position","barcode","depth")
  )
  
  # ---- compute mean depth per barcode ----
  mean_dt <- dt[, .(
    mean_depth = mean(depth, na.rm = TRUE)
  ), by = barcode]
  
  # free memory
  rm(dt); gc()
  
  # ---- read meta ----
  meta_df <- readr::read_csv(meta_path[[sp]], show_col_types = FALSE)
  
  if (!"barcode" %in% colnames(meta_df)) {
    stop(paste("No barcode column in meta for", sp))
  }
  
  # ---- intersect with meta barcodes ----
  mean_dt <- mean_dt[barcode %in% meta_df$barcode]
  
  mean_dt[, species := sp]
  
  mean_dt
}

# -------------------------
# Run all species
# -------------------------
all_depth <- rbindlist(
  lapply(species_order, compute_mean_depth_with_meta)
)

# -------------------------
# QC filter
# -------------------------
qc_depth <- all_depth %>%
  dplyr::filter(mean_depth > 10)

qc_depth$species <- factor(qc_depth$species, levels = species_order)

# -------------------------
# Summary
# -------------------------
summary_df <- qc_depth %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_pass = dplyr::n(),
    median_depth = median(mean_depth),
    .groups = "drop"
  )

print(summary_df)

# -------------------------
# Plot violin
# -------------------------
species.colors <- c("#625D9E", "#CCC9E6","#DCC1DD",
                    "#F2E2B2","#C5DEBA","#AAD4F8","#5184B2")
names(species.colors) <- species_order



p_depth <- ggplot(qc_depth, aes(x = species, y = mean_depth, fill = species)) +
  geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.3) +
  stat_summary(fun = median, geom = "point", size = 1.6, color = "black") +
  scale_fill_manual(values = species.colors) +
  labs(
    x = NULL,
    y = "Mean mtDNA depth"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.y = element_text(color = "black", size = 14)
  )


ggsave("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/mean_mtDepth_percell_violin.pdf",
       p_depth, width = 5.5, height = 6, dpi = 300)
write.table(qc_depth, file = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/mean_mtDepth_percell_data.tsv", sep = "\t", quote = F)

p_depth



# 02 n_var stack barplot -------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(scales)
})

# variant.list <- readRDS(file = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Final_Variant_List_7species_260128.RDS")
# 
# for (sp in names(variant.list)){
#   print(sp)
#   print(length(variant.list[[sp]]))
# }

# --- input ---
df <- tibble::tribble(
  ~Species,    ~n_cells, ~n_before, ~n_after,
  "Human",       13889,      1359,     438,
  "Mouse",       11533,      2209,     255,
  "Rat",         14696,       747,     114,
  "Chicken",      9320,      2944,     216,
  "Cobra",       11832,      5391,     198,
  "Sturgeon",    12350,       896,      29,
  "Stingray",     6625,      3162,     219
) %>%
  dplyr::mutate(
    Species = factor(Species, levels = c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray")),
    n_filtered = n_before - n_after
  )

df_long <- df %>%
  dplyr::select(Species, n_after, n_filtered) %>%
  tidyr::pivot_longer(
    cols = c(n_after, n_filtered),
    names_to = "status",
    values_to = "count"
  ) %>%
  dplyr::mutate(
    status = dplyr::recode(status,
                           n_after = "Retained",
                           n_filtered = "Filtered out")
  )

status_colors <- c(
  "Retained" = "#6281C0",
  "Filtered out" = "grey80"
)

p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Species, y = count, fill = status)) +
  ggplot2::geom_col(width = 0.78, color = "black", linewidth = 0.3) +
  ggplot2::scale_fill_manual(values = status_colors) +
  ggplot2::labs(x = NULL, y = "n_variants", fill = NULL) +
  # label retained counts on top of the retained segment
  ggplot2::geom_text(
    data = df,
    mapping = ggplot2::aes(x = Species, y = n_after*1/2, label = n_after),
    inherit.aes = FALSE,
    size = 3,
    color = "black"
  ) +
  ggplot2::expand_limits(y = max(df$n_before) * 1.12) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(color = "black", size = 12),
    axis.text.y = ggplot2::element_text(color = "black", size = 12),
    axis.title.y = ggplot2::element_text(color = "black", size = 14),
    legend.position = "bottom"
  )

p

ggplot2::ggsave("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/kept_variant.pdf",
                p, width = 5.5, height = 6, units = "in", dpi = 300)



