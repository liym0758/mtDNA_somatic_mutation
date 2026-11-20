suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
})

rna_file <- Sys.getenv("MITOEVO_EXAMPLE_RNA_RDS")
atac_file <- Sys.getenv("MITOEVO_EXAMPLE_ATAC_RDS")
output_file <- Sys.getenv("MITOEVO_EXAMPLE_ANNOTATED_ATAC_RDS")
stopifnot(nzchar(rna_file), nzchar(atac_file), nzchar(output_file))

pbmc.rna <- readRDS(rna_file)
pbmc.atac <- readRDS(atac_file)
pbmc.atac <- FindClusters(pbmc.atac, verbose = FALSE, algorithm = 3,
                          resolution = c(0.08, 1, 1.5))
Idents(pbmc.rna) <- "cell_type"
gene.activities <- GeneActivity(pbmc.atac, features = VariableFeatures(pbmc.rna))
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(pbmc.atac) <- "ACTIVITY"
pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))
transfer.anchors <- FindTransferAnchors(
  reference = pbmc.rna, query = pbmc.atac,
  features = VariableFeatures(pbmc.rna),
  reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca"
)
celltype_predictions <- TransferData(
  anchorset = transfer.anchors, refdata = Idents(pbmc.rna),
  weight.reduction = pbmc.atac[["lsi"]], dims = 2:50
)
ATAC_seurat <- AddMetaData(pbmc.atac, metadata = celltype_predictions)
saveRDS(ATAC_seurat, output_file)

