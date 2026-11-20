suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(rtracklayer)
  library(GenomicRanges)
})

input_root <- Sys.getenv("MITOEVO_EXAMPLE_CELLRANGER_OUTS")
gtf_file <- Sys.getenv("MITOEVO_EXAMPLE_GTF")
output_file <- Sys.getenv("MITOEVO_EXAMPLE_ATAC_RDS")
genome_build <- Sys.getenv("MITOEVO_EXAMPLE_GENOME_BUILD")
stopifnot(nzchar(input_root), nzchar(gtf_file), nzchar(output_file), nzchar(genome_build))

counts <- Read10X_h5(file.path(input_root, "filtered_peak_bc_matrix.h5"))
metadata <- read.csv(file.path(input_root, "singlecell.csv"), header = TRUE, row.names = 1)
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = file.path(input_root, "fragments.tsv.gz"),
  min.cells = 10,
  min.features = 200
)
pbmc <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", meta.data = metadata)
annotations <- as(rtracklayer::import(gtf_file), "GRanges")
genome(annotations) <- genome_build
type_to_biotype <- list(
  transcript = "protein_coding", exon = "protein_coding", CDS = "protein_coding",
  start_codon = "protein_coding", stop_codon = "protein_coding",
  `5UTR` = "protein_coding", `3UTR` = "protein_coding"
)
annotations$gene_biotype <- unlist(lapply(annotations$type, function(x) type_to_biotype[[x]]))
Annotation(pbmc) <- annotations

pbmc <- NucleosomeSignal(pbmc)
pbmc <- TSSEnrichment(pbmc, fast = FALSE)
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
ATAC_seurat <- subset(
  pbmc,
  subset = nCount_peaks > 1500 & nCount_peaks < 15000 &
    nucleosome_signal < 1 & TSS.enrichment > 3 & pct_reads_in_peaks > 40
)
ATAC_seurat <- RunTFIDF(ATAC_seurat)
ATAC_seurat <- FindTopFeatures(ATAC_seurat, min.cutoff = "q0")
ATAC_seurat <- RunSVD(ATAC_seurat)
ATAC_seurat <- RunUMAP(ATAC_seurat, reduction = "lsi", dims = 2:20)
ATAC_seurat <- FindNeighbors(ATAC_seurat, reduction = "lsi", dims = 2:20)
ATAC_seurat <- FindClusters(ATAC_seurat, verbose = FALSE, algorithm = 3, resolution = 0.05)
saveRDS(ATAC_seurat, output_file)

