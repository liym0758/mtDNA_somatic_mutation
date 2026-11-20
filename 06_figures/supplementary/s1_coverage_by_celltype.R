
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

mt.celltype.colors.deep <- c(
  "B cell"    = "#8F3E74",
  "T cell"    = "#5F64A5",
  "NK cell"   = "#8BB6E1",
  "T/NK cell" = "#234B8A",
  "Monocyte"  = "#4E8F8B"
)
keep_ct <- names(mt.celltype.colors.deep)

out_dir <- "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_data <- map_df(names(coverage_files), function(sp) {
  message("Processing species: ", sp)
  
  file_cov  <- coverage_files[sp]
  file_meta <- meta_files[sp]
  
  if (sp == "human") {
    
    seurat.meta <- fread(file_meta, header = TRUE) %>%
      dplyr::filter(sample %in% c("PBMC_1", "PBMC_2")) %>%
      dplyr::select(barcode, main_celltype2) %>%
      dplyr::filter(!is.na(main_celltype2)) %>%
      dplyr::filter(main_celltype2 %in% keep_ct)
    
    df <- fread(file_cov, sep = ",", header = FALSE,
                col.names = c("position", "barcode", "depth"))
    
    df <- df %>%
      dplyr::left_join(seurat.meta, by = "barcode") %>%
      dplyr::filter(!is.na(main_celltype2)) %>%
      dplyr::filter(main_celltype2 %in% keep_ct)
    
  } else {
    seurat.meta <- read.csv(file_meta, row.names = 1, header = TRUE) %>%
      dplyr::as_tibble() %>%
      dplyr::select(barcode, main_celltype2) %>%
      dplyr::filter(!is.na(barcode)) %>%
      dplyr::filter(!is.na(main_celltype2)) %>%
      dplyr::filter(main_celltype2 %in% keep_ct)
    
    df <- fread(file_cov, header = FALSE,
                col.names = c("position", "barcode", "depth"))
    
    df <- df %>%
      dplyr::left_join(seurat.meta, by = "barcode") %>%
      dplyr::filter(!is.na(main_celltype2)) %>%
      dplyr::filter(main_celltype2 %in% keep_ct)
  }
  
  df$position <- as.integer(df$position)
  
  df_summary <- df %>%
    dplyr::group_by(position, main_celltype2) %>%
    dplyr::summarise(
      depth = mean(depth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      species_code = sp,
      species = species_eng[sp],
      celltype = main_celltype2
    ) %>%
    dplyr::select(-main_celltype2)
  
  message("  Summary rows: ", nrow(df_summary))
  return(df_summary)
})

all_data$whole_species_name_expr <- mapply(function(sp, lat) {
  sprintf("%s~(italic(%s))", sp, lat)
}, as.character(all_data$species), latin_names[as.character(all_data$species)])

species_order <- c("Human", "Mouse", "Rat", "Chicken", "Cobra", "Sturgeon", "Stingray")

expr_levels <- sapply(species_order, function(sp) {
  sprintf("%s~(italic(%s))", sp, latin_names[sp])
})

all_data$whole_species_name_expr <- factor(
  all_data$whole_species_name_expr,
  levels = expr_levels
)

species_labels <- c(
  Human    = "Human~(italic(Homo~sapiens))",
  Mouse    = "Mouse~(italic(Mus~musculus))",
  Rat      = "Rat~(italic(Rattus~norvegicus))",
  Chicken  = "Chicken~(italic(Gallus~gallus))",
  Cobra    = "Cobra~(italic(Naja~atra))",
  Sturgeon = "Sturgeon~(italic(Acipenser~sinensis))",
  Stingray = "Stingray~(italic(Potamotrygon~leopoldi))"
)
all_data$species <- factor(all_data$species, levels = names(species_labels))

all_data$celltype <- factor(all_data$celltype, levels = keep_ct)

all_data$log10.depth <- log10(all_data$depth + 1)

all_data$log2.depth <- log2(all_data$depth + 1)

write.table(
  all_data,
  file.path(out_dir, "all_data_raw.by_celltype.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

p <- ggplot(all_data, aes(x = position, y = depth, color = celltype)) +
  geom_line(linewidth = 0.2, alpha = 0.8) +
  scale_color_manual(values = mt.celltype.colors.deep, drop = TRUE) +
  facet_wrap(
    ~ whole_species_name_expr,
    ncol = 1,
    scales = "free_x",
    labeller = label_parsed
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2, color = "grey50") +
  labs(
    x = "Position in mitochondrial genome",
    y = "Average depth per cell",
    color = "Cell type",
    title = "Mitochondrial Genome Coverage Across 7 Species (by Cell Type)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = NA),
    strip.text = element_text(face = "bold", size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.position = "top"
  )

pdf(file.path(out_dir, "species.raw_cov.byCellType.plt.pdf"), width = 7, height = 11)
plot(p)
dev.off()

