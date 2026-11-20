library(tidyverse)
library(data.table)


coverage_files <- c(
  human    = "/PATH/TO/pbmc.mt/SNP/human/mgatk_newCleanBam/merge_pbmc_samples/final/merged.coverage.txt.gz",
  mouse    = "/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/mouse.mgatk.coverage.txt.gz",
  rat      = "/PATH/TO/pbmc.mt/SNP/rat/mgatk_newCleanBam/final/rat.mgatk.coverage.txt.gz",
  chicken  = "/PATH/TO/pbmc.mt/SNP/chicken/mgatk_newCleanBam/final/chicken.mgatk.coverage.txt.gz",
  cobra    = "/PATH/TO/pbmc.mt/SNP/cobra/mgatk_newCleanBam/final/cobra.mgatk.coverage.txt.gz",
  sturgeon = "/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk_newCleanBam/final/sturgeon.mgatk.coverage.txt.gz",
  stingray = "/PATH/TO/pbmc.mt/SNP/stingray/mgatk_newCleanBam/final/stingray.mgatk.coverage.txt.gz"
)

meta_files <- c(
  human    = "/PATH/TO/pbmc.mt/SNP/human/clean_bam/merge_bam/human.PBMConly.sample.barcode.celltype.mtPercent.withArchR.csv",
  mouse    = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/mouse.meta.csv",
  rat      = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/rat.meta.csv",
  chicken  = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/chicken.meta.csv",
  cobra    = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/cobra.meta.csv",
  sturgeon = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/sturgeon.meta.csv",
  stingray = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/stingray.meta.csv"
)

species_eng <- c(
  human    = "Human",
  mouse    = "Mouse",
  rat      = "Rat",
  chicken  = "Chicken",
  cobra    = "Cobra",
  sturgeon = "Sturgeon",
  stingray = "Stingray"
)

latin_names <- c(
  Human     = "Homo~sapiens",
  Mouse     = "Mus~musculus",
  Rat       = "Rattus~norvegicus",
  Chicken   = "Gallus~gallus",
  Cobra     = "Naja~atra",
  Sturgeon  = "Acipenser~sinensis",
  Stingray  = "Potamotrygon~leopoldi"
)

#all_data <- map_df(names(coverage_files), function(sp) {
all_data <- map_df(names(coverage_files), function(sp) {
  message("Processing species: ", sp)
  
  file_cov <- coverage_files[sp]
  file_meta <- meta_files[sp]
  
  if (sp == "human") {
    seurat.meta <- fread(file_meta, header = TRUE)
    seurat.meta <- seurat.meta %>% dplyr::filter(sample %in% c("PBMC_1", "PBMC_2"))
    valid_bc <- paste0(seurat.meta$sample, "_", seurat.meta$barcode)
    df <- fread(file_cov, sep = ",", header = FALSE, 
                col.names = c("position", "barcode", "depth"))
    df <- df[barcode %in% valid_bc]
    
  } else {
    seurat.meta <- read.csv(file_meta, row.names = 1, header = TRUE)
    valid_bc <- seurat.meta$barcode
    df <- fread(file_cov, header = FALSE, 
                col.names = c("position", "barcode", "depth"))
    df <- df[barcode %in% valid_bc]
  }
  
  df$position <- as.integer(df$position)
  
  df_summary <- df %>%
    dplyr::group_by(position) %>%
    dplyr::summarise(
      depth = mean(depth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      species_code = sp,
      species = species_eng[sp]
    )
  
  message("  Summary rows: ", nrow(df_summary))
  return(df_summary)
})

all_data$whole_species_name_expr <- mapply(function(sp, lat) {
  sprintf("%s~(italic(%s))", sp, lat)
}, as.character(all_data$species), latin_names[as.character(all_data$species)])

expr_list <- parse(text = all_data$whole_species_name_expr)

species_order <- c("Human", "Mouse", "Rat", "Chicken", "Cobra", "Sturgeon", "Stingray")

expr_levels <- sapply(species_order, function(sp) {
  sprintf("%s~(italic(%s))", sp, latin_names[sp])
})

all_data$whole_species_name_expr <- factor(
  all_data$whole_species_name_expr,
  levels = expr_levels
)



species_labels <- c(
  Human = "Human~(italic(Homo~sapiens))",
  Mouse = "Mouse~(italic(Mus~musculus))",
  Rat = "Rat~(italic(Rattus~norvegicus))",
  Chicken = "Chicken~(italic(Gallus~gallus))",
  Cobra = "Cobra~(italic(Naja~atra))",
  Sturgeon = "Sturgeon~(italic(Acipenser~sinensis))",
  Stingray = "Stingray~(italic(Potamotrygon~leopoldi))"
)

all_data$species <- factor(all_data$species, levels = names(species_labels))

all_data$log10.depth <- log10(all_data$depth + 1)

write.table(all_data, 
            "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/all_data.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)


p <- ggplot(all_data, aes(x = position, y = log10.depth)) +
  geom_line(color = "#3A6896", linewidth = 0.5, alpha = 1) +
  facet_wrap(
    ~ whole_species_name_expr,
    ncol = 1,
    scales = "free_x",
    labeller = label_parsed
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2, color = "grey50") +
  labs(
    x = "Position in mitochondrial genome",
    y = "Log10(average depth) per cell",
    title = "Mitochondrial Genome Coverage Across 7 Species"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = NA),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

pdf("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/species.cov.plt.pdf",width = 6, height = 10)
plot(p)
dev.off()


df_mean <- all_data %>%
  group_by(species) %>%
  dplyr::summarise(mean_depth = mean(depth, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    species_label = factor(species, levels = c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray"))
  )

p1 <- ggplot(df_mean, aes(x = species_label, y = mean_depth)) +
  geom_col(width = 0.6, fill = "#578FA3", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f", mean_depth)),
            vjust = -0.1, size = 2, color = "black") +
  labs(x = "", y = "Mean Depth Per Position") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title = element_text(face = "bold", size = 10)
  )

pdf("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/species.mean_depth.plt.pdf",width = 6, height = 2)
plot(p1)
dev.off()

write.table(df_mean, 
            "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/barplot_data.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
