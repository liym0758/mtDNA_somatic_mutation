
# get the low coverage positions ----------------------
coverage_data <- read_tsv("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/all_data.tsv")

pdf(file ="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/histogram_depth.pdf", width = 5, height = 3)
for (sp in species_list){
  hist(
    coverage_data %>% 
      dplyr::filter(species_code == sp) %>% 
      dplyr::pull(depth),
    main = sp,
    xlab = "Depth",
    ylab = "Frequency"
  )
}
dev.off()


#check filtered variant position: 

variant.pos.list <- list()
low.coverage.position <- list()

for (sp in species_list){
  print(paste0("-------------------------",sp, "-------------------------"))
  variant.pos.list[[sp]] <- as.integer(
    sub("^([0-9]+).*", "\\1", variant.list[[sp]])
  )
  low.coverage.position[[sp]] <- coverage_data %>% 
    dplyr::filter(species_code == sp & depth < 10)
  print("Number of low coverage position: ")
  print(nrow(low.coverage.position[[sp]]))
  print("Position need to be removed: ")
  print(intersect(variant.pos.list[[sp]], low.coverage.position[[sp]]$position))
}


saveRDS(low.coverage.position, file = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/low.coverage.position.7SpeciesTable.RDS")


# strand balance filter-----------------
# strand_coverage_files <- c(human = "/PATH/TO/pbmc.mt/SNP/human/mgatk_newCleanBam/merge_pbmc_samples/final/cell_strand_position_count_data.tsv", 
#                      mouse = "/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/cell_strand_position_count_data.tsv", 
#                      rat = "/PATH/TO/pbmc.mt/SNP/rat/mgatk_newCleanBam/final/cell_strand_position_count_data.tsv", 
#                      chicken = "/PATH/TO/pbmc.mt/SNP/chicken/mgatk_newCleanBam/final/cell_strand_position_count_data.tsv", 
#                      cobra = "/PATH/TO/pbmc.mt/SNP/cobra/mgatk_newCleanBam/final/cell_strand_position_count_data.tsv", 
#                      sturgeon = "/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk_newCleanBam/final/cell_strand_position_count_data.tsv", 
#                      stingray = "/PATH/TO/pbmc.mt/SNP/stingray/mgatk_newCleanBam/final/cell_strand_position_count_data.tsv" )
# 
# 
# strand.bias.position <- list()
# 
# pdf(
#   file = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/strand_bias_fwdFrac_histogram.pdf",
#   width = 5,
#   height = 4
# )
# 
# for (sp in names(strand_coverage_files)) {
#   
#   cat("Plotting:", sp, "\n")
#   
#   df <- read_tsv(
#     strand_coverage_files[[sp]],
#     col_types = cols()
#   )
#   
#   pos_strand <- df %>%
#     group_by(pos) %>%
#     summarise(
#       fwd_sum = sum(fwd, na.rm = TRUE),
#       rev_sum = sum(rev, na.rm = TRUE),
#       total_depth = fwd_sum + rev_sum,
#       fwd_frac = fwd_sum / pmax(total_depth, 1),
#       .groups = "drop"
#     )
#   
#   p <- ggplot(pos_strand, aes(x = fwd_frac)) +
#     geom_histogram(
#       bins = 50,
#       fill = "grey70",
#       color = "black"
#     ) +
#     geom_vline(xintercept = c(0.45, 0.55),
#                linetype = "dashed",
#                color = "red",
#                linewidth = 0.4) +
#     labs(
#       title = paste0(sp, " position-level strand balance"),
# 
#       x = "Forward strand fraction",
#       y = "Number of positions"
#     ) +
#     theme_classic()
#   
#   print(p)
# }
# 
# dev.off()
# 
# 
# # ---- parameters ----
# bias_lower <- 0.4
# bias_upper <- 0.6
# 
# imbalanced.position <- list()
# 
# for (sp in names(strand_coverage_files)) {
#   
#   cat("==============", sp, "==============\n")
#   
#   df <- read_tsv(
#     strand_coverage_files[[sp]],
#     col_types = cols()
#   )
#   
#   # ---- position-level strand summary ----
#   pos_strand <- df %>%
#     group_by(pos) %>%
#     dplyr::summarise(
#       fwd_sum = sum(fwd, na.rm = TRUE),
#       rev_sum = sum(rev, na.rm = TRUE),
#       total_depth = fwd_sum + rev_sum,
#       bias = fwd_sum / pmax(total_depth, 1),
#       .groups = "drop"
#     )
#   
#   # ---- bias-based imbalanced positions ----
#   imbalanced_pos <- pos_strand %>%
#     dplyr::filter(
#       bias < bias_lower | bias > bias_upper
#     )
#   
#   imbalanced.position[[sp]] <- imbalanced_pos
#   
#   cat("Total positions:", nrow(pos_strand), "\n")
#   cat("Imbalanced positions (bias not in [0.4, 0.6]):",
#       nrow(imbalanced_pos), "\n\n")
# }
# 
# saveRDS(
#   imbalanced.position,
#   file = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/imbalanced.position.7SpeciesTable.RDS"
# )
