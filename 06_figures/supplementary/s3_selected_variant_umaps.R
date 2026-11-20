setwd("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# -------------------- Config --------------------
fixedDepth_rds_root <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Reorder"
umap_root <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/LIS_Score/umap_260128"
mutse_path_template <- "/PATH/TO/pbmc.mt/SNP/%s/mgatk_newCleanBam/downstream/mut_se.RDS"

species.colors <- c("#625D9E", "#AF95C4","#B3659F","#AF651E","#274B2F",  "#5184B2", "#2C4A7F")
names(species.colors) <- c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray")

species.display <- c(
  "human"    = "Human",
  "mouse"    = "Mouse",
  "rat"      = "Rat",
  "chicken"  = "Chicken",
  "cobra"    = "Cobra",
  "sturgeon" = "Sturgeon",
  "stingray" = "Stingray"
)

common_mutation <- list(
  "mouse"    = c("6673C>T","3611T>C"),
  "rat"      = c("4850A>C"),
  "chicken"  = c("7306C>T"),
  "cobra"    = c("8148G>A","326C>T","15082C>T"),
  "stingray" = c("4551C>T","16077C>T")
)

# -------------------- Helper: fixedDepth pass cells (STRICT to your real RDS schema) --------------------
# test$merged_data: tibble with columns:
#   mutation (variant id), barcode (cell id), is_conf_detected_cell (TRUE/FALSE), ...
get_fixedDepth_pass_cells <- function(species, variant) {
  
  rds_file <- file.path(fixedDepth_rds_root, paste0(species, "_Reorder_celltype_VAF_data.RDS"))
  if (!file.exists(rds_file)) stop("Cannot find fixedDepth RDS: ", rds_file)
  
  obj <- readRDS(rds_file)
  
  if (!is.list(obj) || !"merged_data" %in% names(obj)) {
    stop("fixedDepth RDS does not contain merged_data. please check file: ", rds_file)
  }
  
  md <- obj$merged_data
  if (!is.data.frame(md)) stop("merged_data is not data.frame/tibble. file: ", rds_file)
  
  need_cols <- c("mutation", "barcode", "is_conf_detected_cell")
  miss <- setdiff(need_cols, colnames(md))
  if (length(miss) > 0) {
    stop("merged_data missing columns: ", paste(miss, collapse = ", "), " | file: ", rds_file)
  }
  
  pass_cells <- md %>%
    dplyr::filter(.data$mutation == variant, .data$is_conf_detected_cell) %>%
    dplyr::pull(.data$barcode) %>%
    as.character() %>%
    unique()
  
  if (length(pass_cells) == 0) {
    stop("merged_data does not contain pass cells for this variant (mutation==variant & is_conf_detected_cell==TRUE): ",
         variant, " | file: ", rds_file)
  }
  
  pass_cells
}

# -------------------- Plotting --------------------
plot_mito_umap <- function(df, sp_disp, variant, sp_col, pass_cells) {
  
  safe_af <- make.names(variant)
  if (!safe_af %in% colnames(df)) {
    stop("No AF column found in mito_df: ", variant, " | make.names -> ", safe_af)
  }
  
  v <- df[[safe_af]] * 100
  v[!(rownames(df) %in% pass_cells)] <- 0
  df$variant_value <- v
  
  df <- df %>% dplyr::arrange(variant_value)
  
  ggplot(df, aes(x = umap_1, y = umap_2, color = variant_value)) +
    geom_point(size = 0.5) +
    scale_color_gradientn(colors = c("lightgrey", sp_col), limits = c(0, 100)) +
    ggtitle(paste0(sp_disp, " | ", variant, " (fixedDepth: only pass cells colored)")) +
    labs(color = "Frequency (%)") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, color = "black")
    )
}

run_umap_fixedDepth_overlay0 <- function(species) {
  
  if (!species %in% names(species.display)) stop("species.display missing species: ", species)
  if (!species %in% names(common_mutation)) stop("common_mutation does not define this species: ", species)
  
  sp_disp <- species.display[[species]]
  if (!sp_disp %in% names(species.colors)) stop("species.colors missing: ", sp_disp)
  sp_col <- unname(species.colors[[sp_disp]])
  
  outdir <- file.path(umap_root, species)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  message("Output directory: ", outdir)
  
  # Seurat (UMAP embedding + meta)
  objlist <- readRDS("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/RDS/objlist.RDS")
  if (!species %in% names(objlist)) stop("objlist.RDS does not contain this species Seurat: ", species)
  species.seurat <- objlist[[species]]
  
  mut_se_path <- sprintf(mutse_path_template, species)
  if (!file.exists(mut_se_path)) stop("Cannot find mut_se.RDS: ", mut_se_path)
  mut_se <- readRDS(mut_se_path)
  if (!inherits(mut_se, "SummarizedExperiment")) stop("mut_se.RDS is not SummarizedExperiment: ", mut_se_path)
  
  keep_cells <- intersect(colnames(mut_se), colnames(species.seurat))
  if (length(keep_cells) == 0) stop("mut_se and seurat have no shared cells: ", species)
  mut_se <- mut_se[, keep_cells]
  species.seurat <- subset(species.seurat, cells = keep_cells)
  
  # variants
  variants_to_plot <- common_mutation[[species]]
  
  if (!"allele_frequency" %in% names(assays(mut_se))) stop("mut_se missing allele_frequency assay")
  all_vars <- rownames(assays(mut_se)[["allele_frequency"]])
  variants_to_plot <- intersect(variants_to_plot, all_vars)
  if (length(variants_to_plot) == 0) {
    message("mut_se contains no common mutations, skipping: ", species)
    return(invisible(NULL))
  }
  
  # AF: variant x cell
  af_mat <- assays(mut_se)[["allele_frequency"]][variants_to_plot, , drop = FALSE]
  
  emb <- as.data.frame(species.seurat@reductions$umap@cell.embeddings)
  colnames(emb) <- c("umap_1", "umap_2")
  meta <- as.data.frame(species.seurat@meta.data)
  
  af_df <- as.data.frame(t(af_mat))
  colnames(af_df) <- make.names(rownames(af_mat), unique = TRUE)
  
  mito_df <- cbind(meta, emb, af_df)
  rownames(mito_df) <- rownames(meta)
  
  for (variant in variants_to_plot) {
    
    pass_cells <- get_fixedDepth_pass_cells(species = species, variant = variant)
    pass_cells <- intersect(pass_cells, rownames(mito_df))
    
    if (length(pass_cells) == 0) {
      message("pass_cells has no overlap with current object, skipping: ", species, " | ", variant)
      next
    }
    
    safe_variant <- gsub(">", "_", variant)
    safe_variant <- gsub(":", "-", safe_variant)
    
    pdf_file <- file.path(outdir, paste0(species, "_", safe_variant, ".umap.fixedDepthOverlay0.pdf"))
    pdf(pdf_file, width = 6, height = 5)
    print(plot_mito_umap(mito_df, sp_disp, variant, sp_col, pass_cells))
    dev.off()
    
    message("Saved: ", basename(pdf_file),
            " | total_cells = ", nrow(mito_df),
            " | pass_cells = ", length(pass_cells))
  }
  
  invisible(outdir)
}

# -------------------- Batch run --------------------
species_to_analyze <- c("mouse", "rat", "chicken", "cobra", "stingray")
for (sp in species_to_analyze) {
  run_umap_fixedDepth_overlay0(sp)
}
