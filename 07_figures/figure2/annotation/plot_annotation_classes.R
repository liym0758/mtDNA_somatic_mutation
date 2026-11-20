

celltype_AnnoVar_BarPlt <- function(cell_frequency_matrix_file, 
                                    cell_annotation_file, #VAF RDS[[1]]
                                    annovar_annotation_file, 
                                    gene_pred_file, 
                                    mtDNA_length = 16299,
                                    filtered_mutations = NULL,
                                    output_file = NULL,
                                    color_palette = NULL,
                                    width = 10,
                                    height = 7,
                                    min_frequency = 0) {
  library(tidyverse)
  
  message("Reading mutation matrix data...")
  mutation_matrix <- read_tsv(cell_frequency_matrix_file) 
  
  message("Reading cell-type annotations...")
  cell_annotations <- cell_annotation_file %>%
    dplyr::select(barcode, main_celltype2) %>%
    distinct(barcode, .keep_all = TRUE)
  
  colnames(cell_annotations) <- c("barcode","celltype")
  
  valid_cells <- intersect(unique(mutation_matrix$barcode), cell_annotations$barcode)
  mutation_matrix <- mutation_matrix %>% 
    dplyr::filter(barcode %in% valid_cells)
  cell_annotations <- cell_annotations[cell_annotations$barcode %in% valid_cells, ]
  
  if (!is.null(filtered_mutations)) {
    if (!is.null(filtered_mutations)) {
      mutation_matrix <- mutation_matrix %>% 
        dplyr::filter(mutation %in% filtered_mutations)
    }
  }
  
  message("Processing genePred and computing region lengths...")
  gene_pred <- read_tsv(gene_pred_file, 
                        col_names = c(
    "bin", "transcript", "chr", "strand", "txStart", "txEnd", 
    "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", 
    "score", "gene", "cdsStartStat", "cdsEndStat", "exonFrames"
  ),
  col_types = cols(
    bin = col_double(),
    transcript = col_character(),
    chr = col_character(),
    strand = col_character(),
    txStart = col_double(),
    txEnd = col_double(),
    cdsStart = col_double(),
    cdsEnd = col_double(),
    exonCount = col_double(),
    exonStarts = col_character(),
    exonEnds = col_character(),
    score = col_double(),
    gene = col_character(),
    cdsStartStat = col_character(),
    cdsEndStat = col_character(),
    exonFrames = col_character()
  ))
  
  last_gene_end <- max(gene_pred$txEnd)
  print(paste("last gene end:",last_gene_end))
  
  first_gene_start <- min(gene_pred$txStart)
  print(paste("first gene start:",first_gene_start))
  
  gene_ranges <- GRanges(
    seqnames = "MT",
    ranges = IRanges(start = gene_pred$txStart + 1, end = gene_pred$txEnd),
    strand = gene_pred$strand,
    gene = gene_pred$gene,
    cdsStart = gene_pred$cdsStart,
    cdsEnd = gene_pred$cdsEnd
  )
  

  calculate_region_length <- function() {
    gene_pred_processed <- gene_pred %>%
      dplyr::mutate(
        processed_gene = gsub("^(mt-|MT-)", "", gene, ignore.case = TRUE),
        gene_type = case_when(
          cdsStart != cdsEnd ~ "protein_coding",
          grepl("^[Tt]", processed_gene) ~ "tRNA",
          grepl("^[Rr]", processed_gene) ~ "rRNA",
          TRUE ~ "protein_coding"
        )
      )
    
    safe_calc_length <- function(genes_df, seqname = "MT") {
      if (nrow(genes_df) == 0) {
        return(0)
      }
      ranges <- GRanges(
        seqnames = rep(seqname, nrow(genes_df)),
        ranges = IRanges(start = genes_df$txStart + 1, end = genes_df$txEnd),
        strand = genes_df$strand
      )
      return(sum(width(reduce(ranges))))
    }
    
    protein_genes <- gene_pred_processed %>% dplyr::filter(gene_type == "protein_coding")
    cds_ranges <- GRanges(
      seqnames = rep("MT", nrow(protein_genes)),
      ranges = IRanges(start = protein_genes$cdsStart + 1, end = protein_genes$cdsEnd),
      strand = protein_genes$strand
    )
    cds_length <- ifelse(nrow(protein_genes) > 0, sum(width(reduce(cds_ranges))), 0)
    
    utr5_length <- 0
    if (nrow(protein_genes) > 0) {
      utr5_ranges <- GRanges(
        seqnames = rep("MT", nrow(protein_genes)),
        ranges = IRanges(
          start = ifelse(protein_genes$strand == "+", protein_genes$txStart + 1, protein_genes$cdsEnd + 1),
          end = ifelse(protein_genes$strand == "+", protein_genes$cdsStart, protein_genes$txEnd)
        ),
        strand = protein_genes$strand
      )
      utr5_length <- sum(width(reduce(utr5_ranges)))
    }
    
    utr3_length <- 0
    if (nrow(protein_genes) > 0) {
      utr3_ranges <- GRanges(
        seqnames = rep("MT", nrow(protein_genes)),
        ranges = IRanges(
          start = ifelse(protein_genes$strand == "+", protein_genes$cdsEnd + 1, protein_genes$txStart + 1),
          end = ifelse(protein_genes$strand == "+", protein_genes$txEnd, protein_genes$cdsStart)
        ),
        strand = protein_genes$strand
      )
      utr3_length <- sum(width(reduce(utr3_ranges)))
    }
    
    trna_genes <- gene_pred_processed %>% dplyr::filter(gene_type == "tRNA")
    trna_length <- safe_calc_length(trna_genes)
    
    rrna_genes <- gene_pred_processed %>% dplyr::filter(gene_type == "rRNA")
    rrna_length <- safe_calc_length(rrna_genes)
    
    other_ncrna_genes <- gene_pred_processed %>% dplyr::filter(gene_type == "other_ncRNA")
    other_ncrna_length <- safe_calc_length(other_ncrna_genes)
    
    dloop_length <- mtDNA_length - last_gene_end + first_gene_start
    
    all_genes_ranges <- GRanges(
      seqnames = rep("MT", nrow(gene_pred)),
      ranges = IRanges(start = gene_pred$txStart + 1, end = gene_pred$txEnd),
      strand = gene_pred$strand
    )
    gene_coverage <- reduce(all_genes_ranges)
    intergenic_length <- mtDNA_length - sum(width(gene_coverage))# - dloop_length
    
    result_df <- data.frame(
      mutation_type = c("protein-coding", "5'UTR", "3'UTR", 
                        "tRNA", "rRNA", "ncRNA", "intergenic"),#, "D-loop"),
      length = c(cds_length, utr5_length, utr3_length, 
                 trna_length, rrna_length, other_ncrna_length, 
                 intergenic_length)#, dloop_length)
    )
    
    result_df$length <- pmax(result_df$length, 0)
    
    return(result_df)
  }
  
  region_lengths <- calculate_region_length()
  message("Region-length calculation complete:")
  print(region_lengths)
  

  
  message("Processing ANNOVAR annotations...")
  
  annovar_annotations <- read_tsv(annovar_annotation_file) %>%
    dplyr::mutate(
      position = as.numeric(Start),
      mutation = paste0(position, Ref, ">", Alt),
      is_dloop = (position > last_gene_end)|(position < first_gene_start),
      
      processed_gene = gsub("^(mt-|MT-)", "", Gene.ensGene, ignore.case = TRUE),
      gene_first_char = ifelse(!is.na(processed_gene), substr(processed_gene, 1, 1), NA),
      
      mutation_type = case_when(
        Func.ensGene %in% c("upstream;downstream","upstream","downstream") & (position > last_gene_end | position < first_gene_start) ~ "intergenic",#"D-loop" ,
        Func.ensGene %in% c("upstream;downstream","upstream","downstream") ~ "intergenic",
        
        position > last_gene_end ~ "intergenic",#"D-loop",
        
        Func.ensGene == "exonic" ~ "protein-coding",
        Func.ensGene == "ncRNA_exonic" & grepl("^[Tt]", processed_gene) ~ "tRNA",
        Func.ensGene == "ncRNA_exonic" & grepl("^[Rr]", processed_gene) ~ "rRNA",
        Func.ensGene == "ncRNA_exonic" ~ "ncRNA",
        Func.ensGene %in% c("intronic", "ncRNA_intronic") ~ "intronic",
        Func.ensGene %in% c("upstream", "downstream", "intergenic") ~ "intergenic",
        Func.ensGene == "splicing" ~ "splicing",
        Func.ensGene == "UTR5" ~ "5'UTR",
        Func.ensGene == "UTR3" ~ "3'UTR",
        TRUE ~ "intergenic"
      )
    ) %>%
    dplyr::select(mutation, mutation_type) %>%
    distinct()
  
  message("Mutation-type summary:")

  classification_stats_simple <- annovar_annotations %>%
    group_by(mutation_type) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    ungroup()
  
  print("Simple category summary:")
  print(classification_stats_simple)

  
  if (any(is.na(annovar_annotations$mutation_type))) {
    unclassified <- annovar_raw %>%
      dplyr::filter(mutation %in% annovar_annotations$mutation[is.na(annovar_annotations$mutation_type)]) %>%
      distinct(Func.ensGene)
    warning("Unclassified Func.ensGene values: ", paste(unclassified$Func.ensGene, collapse = ", "))
  }
  
  message("Reshaping mutation matrix to long format...")
  mutation_data <- mutation_matrix %>%
    dplyr::filter(frequency > min_frequency) %>%
    left_join(cell_annotations, by = "barcode") %>%
    dplyr::select(celltype, mutation) %>%
    distinct()
  
  mutation_with_type <- mutation_data %>%
    left_join(annovar_annotations, by = "mutation") %>%
    dplyr::filter(!is.na(mutation_type)) %>%
    distinct(celltype, mutation, mutation_type)
  
  message("Summarizing mutation types by cell type (with length correction)...")
  celltype_mutation_stats <- mutation_with_type %>%
    group_by(celltype, mutation_type) %>%
    dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
    left_join(region_lengths, by = "mutation_type") %>%
    dplyr::mutate(
      density = count / (length/mtDNA_length),
      density = replace_na(density, 0)
    ) %>%
    group_by(celltype) %>%
    dplyr::mutate(
      total_density = sum(density),
      proportion = density / total_density
    ) %>%
    ungroup()
  
  
  if (is.null(color_palette)) {
    color_palette <- c(
      `protein-coding` = "#d62728",
      splicing = "#e377c2",
      ncRNA = "#17becf",
      intronic = "#bcbd22",
      intergenic = "#D4D4D4",#"#7f7f7f",
      `5'UTR` =  "#aec7e8",
      `3'UTR` = "#ff7f0e"#,
    )
  }
  
  mutation_order <- c(
    "protein-coding",
    "tRNA", "rRNA",
    "splicing", "5'UTR", "3'UTR",
    "intronic", "intergenic"#,"D-loop"
  )
  
  celltype_mutation_stats <- celltype_mutation_stats %>%
    dplyr::mutate(mutation_type = factor(mutation_type, levels = intersect(mutation_order, unique(mutation_type))))
  
  celltype_order <-  c("B cell","T cell","T/NK cell","NK cell","Monocyte")
  celltype_levels <- intersect(celltype_order, unique(celltype_mutation_stats$celltype))
  celltype_mutation_stats <- celltype_mutation_stats %>% 
    dplyr::filter(celltype %in% celltype_levels)
  celltype_mutation_stats$celltype <- factor(celltype_mutation_stats$celltype, levels = rev(celltype_levels))
  
  
  message("Generating plot...")
  p <- ggplot(celltype_mutation_stats, 
              aes(x = mutation_type, 
                  y = density,
                  fill = celltype)) +
    geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
    scale_fill_manual(values = color_palette) +
    scale_y_continuous() +
    labs(
      title = paste(str_to_title(species_list[i]),"- Mutation Annotation"),
      x = "Region",
      y = "Mutation enrichment (normalized by region length)",
      fill = "Cell Type"
    ) +
    coord_flip()+
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 9),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = width, height = height, units = "in", dpi = 300)
    message("Plot saved to: ", output_file)
  }
  
  return(list(
    stats = celltype_mutation_stats,
    plot = p
  ))
}






annovar_color <- c(
  `protein-coding` = "#625D9E",
  tRNA = "#C7B1D4",
  rRNA = "#B3659F",
  intergenic = "#D3D3D3"#"#A63376",
)


base_path <- "/PATH/TO/pbmc.mt/SNP/"

species_list <- c("human","mouse","rat","chicken","cobra","sturgeon","stingray")

annovar_annotation_file_list <- c(
  "human" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/human/result/hg38_anno.hg38_chrM_multianno.txt",
  "mouse" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/mouse/result/mm39_anno.mm39_MT_multianno.txt",
  "rat" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/rat/result/rn7_MT_anno.rn7_MT_multianno.txt",
  "chicken" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/chicken/result/gg1_MT_anno.gg1_MT_multianno.txt",
  "cobra" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/cobra/result/na_MT_anno.na_MT_multianno.txt",
  "sturgeon" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/sturgeon/result/as_CNCB_MT_anno.as_CNCB_MT_multianno.txt",
  "stingray" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/stingray/result/pl_MT_anno.pl_MT_multianno.txt"
)
gene_pred_file_list <- c(
  "human" = "/PATH/TO/annovar/humandb/hg38_chrM_ensGene.txt",
  "mouse" = "/PATH/TO/annovar/mm39_MT_db/mm39_MT_ensGene.txt",
  "rat" = "/PATH/TO/annovar/rat_MT_db/rn7_MT_ensGene.txt",
  "chicken" = "/PATH/TO/annovar/chicken_MT_db/gg1_MT_ensGene.txt",
  "cobra" = "/PATH/TO/annovar/cobra_MT_db/na_MT_ensGene.txt",
  "sturgeon" = "/PATH/TO/annovar/sturgeon_MT_db/as_CNCB_MT_ensGene.txt",
  "stingray" = "/PATH/TO/annovar/stingray_MT_db/pl_MT_ensGene.txt"
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

meta_path <- c(
  "human"="/PATH/TO/pbmc.mt/SNP/human/clean_bam/merge_bam/human.PBMConly.sample.barcode.celltype.mtPercent.withArchR.csv",
  "mouse"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/mouse.meta.csv",
  "rat"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/rat.meta.csv",
  "chicken"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/chicken.meta.csv",
  "cobra"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/cobra.meta.csv",
  "sturgeon"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/sturgeon.meta.csv",
  "stingray"="/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/stingray.meta.csv"
)


get_mtDNA_length_from_fasta <- function(fasta_file) {
  fasta_seq <- readDNAStringSet(fasta_file)
  
  seq_length <- width(fasta_seq)[1]
  
  message("mtDNA length detected from FASTA: ", seq_length)
  return(seq_length)
}

filtered_mutation_list <- readRDS("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Final_Variant_List_7species.RDS")
result <- list()
for (i in 1:7){
  print(paste("==============",species_list[i],"=============="))
  
  filtered_mutataions <- unique(filtered_mutation_list[[i]])
  
  seq_length <- get_mtDNA_length_from_fasta(ref_fasta[[i]])
  
  result[[i]] <- celltype_AnnoVar_BarPlt(
    cell_frequency_matrix_file = paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/",species_list[i],"_detail_filtered_df.tsv"),
    cell_annotation_file = read.csv(meta_path[i],header = T),
    gene_pred_file = gene_pred_file_list[i],
    annovar_annotation_file = annovar_annotation_file_list[i],
    mtDNA_length = seq_length,
    filtered_mutations = filtered_mutataions,
    color_palette = annovar_color
  )
  
  pdf(file = paste0(base_path,"mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/",species_list[i],"/result/CellType_AnnoVar_density_CombineDLoop_BarPlt.pdf"), width = 6, height = 3)
  print(result[[i]]$plot)
  dev.off()
  
   write.table(result[[i]]$stats, file = paste0(base_path,"mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/",species_list[i],"/result/CellType_AnnoVar_CombineDLoop_BarPlt.tsv"), sep = "\t", quote = F)
}



multianno <- list()
for (i in 1:7){
  print(species_list[i])
  multianno[[i]] <- read.table(annovar_annotation_file_list[i], sep = "\t", header = T)
  print(unique(multianno[[i]]$Func.ensGene))
}




