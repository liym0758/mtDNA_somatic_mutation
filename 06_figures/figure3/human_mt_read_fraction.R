setwd("/PATH/TO/pbmc.mt/SNP/human/clean_bam/merge_bam/")
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# ================================
# ================================
anno_file <- "/PATH/TO/pbmc.mt/data/human_bottleneck/output_pbmc/proj_PBMC_only/proj_pbmc_annotation.tsv"

anno <- read_tsv(
  anno_file
)

# ================================
# ================================
anno <- anno %>%
  mutate(
    sample_raw = str_extract(raw_barcode, "^[^#]+"),

    sample = sample_raw %>%
      str_replace("_rep", "_"),

    barcode = str_extract(raw_barcode, "(?<=#)[ACGT]+-1")
  )

table_C <- anno %>% select(sample, barcode, celltype0, sub_celltype, main_celltype)

# ================================
# ================================
sample_paths <- list(
  PBMC_1 = "/PATH/TO/pbmc.mt/SNP/human/cellranger/PBMC_1/PBMC_1/outs/singlecell.csv",
  PBMC_2 = "/PATH/TO/pbmc.mt/SNP/human/cellranger/PBMC_2/PBMC_2/outs/singlecell.csv"#,
  #CD34_1 = "/PATH/TO/pbmc.mt/SNP/human/cellranger/CD34_1/CD34_1/outs/singlecell.csv",
  #CD34_2 = "/PATH/TO/pbmc.mt/SNP/human/cellranger/CD34_2/CD34_2/outs/singlecell.csv"
)

metrics_list <- imap(sample_paths, function(path, sample) {
  read_csv(path) %>%
    select(barcode = barcode, total_reads = total, mt_reads = mitochondrial) %>%
    mutate(sample = sample)
})

metrics_df <- bind_rows(metrics_list)

# ================================
# ================================
final_df <- table_C %>%
  left_join(metrics_df, by = c("sample", "barcode")) %>%
  mutate(mt.percent = mt_reads / total_reads * 100)

write_tsv(final_df, "human.PBMConly.sample.barcode.celltype.mtPercent.tsv")



library(dplyr)
library(stringr)

mt <- read.table("/PATH/TO/pbmc.mt/SNP/human/clean_bam/merge_bam/human.PBMConly.sample.barcode.celltype.mtPercent.tsv",
                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

mt <- mt %>% 
  mutate(sample_barcode = paste(sample, barcode, sep = "_"))

# ──────────────────────────────────────
# ──────────────────────────────────────

#proj_res1 <- readRDS("/PATH/TO/pbmc.mt/data/human_bottleneck/output_pbmc/proj_res1_annotated/lym_anno-ArchR-Project.rds")

archr_df <- as.data.frame(proj_pbmc@cellColData)

archr_df$archr_barcode_original <- rownames(archr_df)

archr_df$sample_barcode <- gsub("#", "_", archr_df$archr_barcode_original)

 head(archr_df$sample_barcode)

# ──────────────────────────────────────
# ──────────────────────────────────────
archr_df2 <- archr_df %>% 
  select(sample_barcode, sub_celltype, main_celltype)

#archr_df2$sample_barcode <- gsub("^CD34_rep1", "CD34_1", archr_df2$sample_barcode)
archr_df2$sample_barcode <- gsub("^PBMC_rep1", "PBMC_1", archr_df2$sample_barcode)
archr_df2$sample_barcode <- gsub("^PBMC_rep2", "PBMC_2", archr_df2$sample_barcode)
#archr_df2$sample_barcode <- gsub("^CD34_rep2", "CD34_2", archr_df2$sample_barcode)

# ──────────────────────────────────────
# 5. merge
# ──────────────────────────────────────
mt_merged <- merge(mt, archr_df2, by = c("sample_barcode"))

head(mt_merged)

out_file <- "/PATH/TO/pbmc.mt/SNP/human/clean_bam/merge_bam/human.PBMConly.sample.barcode.celltype.mtPercent.withArchR.tsv"

write.table(mt_merged, file = out_file,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Saved merged file to:\n", out_file, "\n")
