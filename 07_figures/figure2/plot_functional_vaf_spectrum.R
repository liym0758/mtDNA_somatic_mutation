# =========================
# 3) plot_vaf_spectrum_by_type ------------
# =========================

read_genepred <- function(path) {
  message("Reading GenePred file: ", path)
  
  gp <- fread(path, header = FALSE, sep = "\t")
  
  message("Detected column count: ", ncol(gp))
  
  if (ncol(gp) >= 12) {
    colnames(gp)[1:12] <- c("name", "chrom", "strand", "txStart", "txEnd", 
                            "cdsStart", "cdsEnd", "exonCount", "exonStarts", 
                            "exonEnds", "score", "name2")
    
    if (ncol(gp) >= 15) {
      colnames(gp)[13:15] <- c("cdsStartStat", "cdsEndStat", "exonFrames")
    }
    if (ncol(gp) >= 16) {
      colnames(gp)[1] <- "bin"
      colnames(gp)[2:16] <- c("name", "chrom", "strand", "txStart", "txEnd", 
                              "cdsStart", "cdsEnd", "exonCount", "exonStarts", 
                              "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
    }
  }
  
  message("Column names used: ", paste(colnames(gp), collapse = ", "))
  
  gp <- gp %>%
    filter(toupper(chrom) %in% c("MT", "CHRM", "CHRMT", "M") | 
             grepl("MT|MITO", chrom, ignore.case = TRUE))
  
  numeric_cols <- c("txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "score")
  for (col in numeric_cols) {
    if (col %in% names(gp)) {
      gp[[col]] <- as.numeric(gp[[col]])
    }
  }
  
  gp <- gp %>%
    mutate(
      gene_name = ifelse(!is.na(name2) & name2 != "" & name2 != "none", name2, name),
      gene_length = txEnd - txStart + 1
    )
  
  gp$gene_type <- character(nrow(gp))
  for (i in 1:nrow(gp)) {
    gp$gene_type[i] <- classify_gene_type_genepred(gp$gene_name[i], gp$gene_length[i])
  }
  
  gp_genes <- gp
  unique_genes <- unique(gp_genes$gene_name)
  result_list <- list()
  
  for (gene in unique_genes) {
    gene_subset <- gp_genes[gp_genes$gene_name == gene, ]
    max_length_idx <- which.max(gene_subset$gene_length)
    if (length(max_length_idx) > 0) {
      result_list[[gene]] <- gene_subset[max_length_idx[1], ]
    }
  }
  
  gp_genes <- do.call(rbind, result_list)
  rownames(gp_genes) <- NULL
  
  gp_genes <- gp_genes[, c("gene_name", "chrom", "strand", "txStart", "txEnd", 
                           "exonStarts", "exonEnds", "exonCount", "gene_length", "gene_type")]
  
  message("Identified in GenePred file ", nrow(gp_genes), " mitochondrial genes")
  message("Gene-type summary:")
  print(table(gp_genes$gene_type))
  
  if (nrow(gp_genes) > 0) {
    message("First five gene examples:")
    print(head(gp_genes[, c("gene_name", "gene_type", "gene_length")], 5))
  }
  
  return(gp_genes)
}

classify_gene_type_genepred <- function(gene_name, gene_length = NA) {
  gene_name <- as.character(gene_name)
  if (is.na(gene_name)) return("ncRNA")
  
  gene_upper <- toupper(gene_name)
  
  if (grepl("^MT-T|^T", gene_upper)) return("tRNA")
  
  if (grepl("^MT-R|^R|16S|12S", gene_upper)) return("rRNA")
  
  if (grepl("^MT-", gene_upper)) {
    if (grepl("^MT-ND|^MT-CO|^MT-ATP|^MT-CY", gene_upper)) {
      return("CDS")
    }
  }
  
  cds_pattern <- "ND[1-6]|CO[I-III1-3]|ATP[6-8]|CYB|CYTB|NAD[1-6]|COX[1-3]"
  if (grepl(cds_pattern, gene_upper)) return("CDS")
  
  if (!is.na(gene_length)) {
    if (gene_length >= 50 && gene_length <= 100) return("tRNA")
    if (gene_length >= 500 && gene_length <= 3000) return("rRNA")
    if (gene_length >= 100 && gene_length <= 2000) return("CDS")
  }
  
  "ncRNA"
}


read_and_classify_annovar_with_gp <- function(annovar_file, genepred_file){
  gp <- read_genepred(genepred_file)
  if(nrow(gp)==0) stop("GenePred file is empty")
  first_start <- min(gp$txStart + 1)
  last_end <- max(gp$txEnd)
  mt_len <- last_end
  
  ann <- readr::read_tsv(annovar_file, show_col_types = FALSE)
  names(ann) <- gsub("\\.", "_", names(ann))
  
  Func_col <- if("Func_ensGene" %in% names(ann)) "Func_ensGene" else grep("Func", names(ann), value = TRUE)[1]
  Gene_col <- if("Gene_ensGene" %in% names(ann)) "Gene_ensGene" else grep("Gene", names(ann), value = TRUE)[1]
  Exonic_col <- if("ExonicFunc_ensGene" %in% names(ann)) "ExonicFunc_ensGene" else grep("ExonicFunc", names(ann), value = TRUE)[1]
  
  ann2 <- ann %>% mutate(
    position = as.integer(Start),
    mutation_id = paste0(Start,Ref,">",Alt),
    processed_gene = ifelse(!is.na(.data[[Gene_col]]), gsub("^(mt-|MT-)","", .data[[Gene_col]],ignore.case = TRUE),NA),
    Func_col2 = .data[[Func_col]],
    Exonic_col2 = .data[[Exonic_col]]
  )
  
  ann_class <- ann2 %>%
    rowwise() %>%
    mutate(
      mutation_type = {
        pos <- position
        if (!is.na(pos) && Func_col2 %in% c("upstream;downstream", "upstream", "downstream") &&
            (pos > last_end || pos < first_start)) {
          "D-loop"
        } else if (!is.na(Func_col2) && Func_col2 %in% c("upstream;downstream", "upstream", "downstream", "intergenic")) {
          "intergenic"
        } else if (!is.na(Exonic_col2) && grepl("^(nonsynonym|stopgain|stoploss|startloss)", Exonic_col2, ignore.case = TRUE)) {
          "nonsynonymous"
        } else if (!is.na(Exonic_col2) && grepl("^synonym", Exonic_col2, ignore.case = TRUE)) {
          "synonymous"
        } else if (!is.na(Func_col2) && Func_col2 == "ncRNA_exonic" &&
                   !is.na(processed_gene) && grepl("^[Tt]", processed_gene)) {
          "tRNA"
        } else if (!is.na(Func_col2) && Func_col2 == "ncRNA_exonic" &&
                   !is.na(processed_gene) && grepl("^[Rr]", processed_gene)) {
          "rRNA"
        } else {
          "other"
        }
      }
    ) %>%
    ungroup() %>%
    select(mutation_id, mutation_type) %>%
    distinct()
  
  
  return(ann_class)
}

plot_vaf_spectrum_by_type <- function(species, annovar_file, genepred_file, vaf_file, output_dir="."){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ann <- read_and_classify_annovar_with_gp(annovar_file, genepred_file)
  ann <- ann %>%
    dplyr::mutate(
      mutation_type = dplyr::if_else(
        mutation_type == "D-loop",
        "intergenic",
        mutation_type
      )
    )
  
  vaf <- readr::read_tsv(vaf_file, show_col_types = FALSE) %>% dplyr::mutate(mutation_id = mutation)
  
  df <- vaf %>%
    left_join(ann, by = "mutation_id") %>%
    dplyr::mutate(
      mutation_type = if_else(is.na(mutation_type), "unknown", mutation_type),
      vaf_bin = cut(frequency, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
    ) %>%
    group_by(mutation_type, vaf_bin) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    group_by(mutation_type) %>%
    dplyr::mutate(freq = n / sum(n)) %>% 
    ungroup()
  
  bins <- levels(cut(0.05, breaks=seq(0,1,0.1), include.lowest=TRUE))
  # bin_labels <- bins
  
  
  mutation_types <- unique(df$mutation_type)
  full_df <- expand_grid(
    vaf_bin=factor(bins,levels=bins), 
    mutation_type=mutation_types) %>%
    left_join(df, by=c("vaf_bin","mutation_type")) %>%
    dplyr::mutate(n=replace_na(n,0), freq=replace_na(freq,0))
  

  mutation_levels <- c(
    "synonymous",
    "nonsynonymous",
    "rRNA",
    "tRNA",
    "intergenic"
  )
  
  full_df <- full_df %>%
    dplyr::mutate(
      mutation_type = factor(mutation_type, levels = mutation_levels)
    )
  
  p <- ggplot(full_df, aes(x=vaf_bin, y=freq, color=mutation_type, group=mutation_type)) +
    geom_line(size=1) + 
    geom_point(size=1.5) + 
    theme_bw(base_size=12) +
    theme(panel.grid = element_blank(),
          text = element_text(color = "black"),
          axis.text.x = element_text(color = "black",angle = 45, hjust = 1),
          axis.text.y = element_text(color = "black"),
          axis.title.x = element_text(color = "black"),
          axis.title.y = element_text(color = "black"),
          
          plot.title = element_text(color = "black",hjust = 0.5),
          
          legend.title = element_text(color = "black"),
          legend.text = element_text(color = "black")
           )+
    labs(title=paste0(str_to_title(species)),
         x="Variant Allele Frequency (VAF)",
         y="Proportion of mutations",
         color="Type") + 
    scale_color_manual(
      values = c(
        "synonymous" = "#E78AC3",
        "nonsynonymous" = "#FC8D62",
        "rRNA" = "#8DA0CB",
        "tRNA" = "#A6D854",
        "intergenic" = "#66C2A5"
      ))
  
  ggsave(file.path(output_dir,paste0(species,"_VAF_spectrum_by_type.pdf")), p, width=6, height=3.5)
  
  # NS / SY count per bin
  ns_vec <- full_df %>% dplyr::filter(mutation_type=="nonsynonymous") %>% arrange(vaf_bin) %>% pull(n)
  sy_vec <- full_df %>% dplyr::filter(mutation_type=="synonymous") %>% arrange(vaf_bin) %>% pull(n)
  if(length(ns_vec)==0) ns_vec <- rep(0,10)
  if(length(sy_vec)==0) sy_vec <- rep(0,10)
  bins_num <- seq(0.1,1.0,0.1)
  
  return(list(df=full_df, NS_count=ns_vec, SY_count=sy_vec, bins=bins_num))
}

outdir <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/plots/VAF_vs_MutationType/Single_plot/new_human/" 

species_list <- c('human','mouse','rat','chicken','cobra','sturgeon','stingray')

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
  "human" = "/PATH/TO/biosoft/annovar_2020/humandb/hg38_chrM_ensGene.txt",
  "mouse" = "/PATH/TO/biosoft/annovar_2020/mm39_MT_db/mm39_MT_ensGene.txt",
  "rat" = "/PATH/TO/biosoft/annovar_2020/rat_MT_db/rn7_MT_ensGene.txt",
  "chicken" = "/PATH/TO/biosoft/annovar_2020/chicken_MT_db/gg1_MT_ensGene.txt",
  "cobra" = "/PATH/TO/biosoft/annovar_2020/cobra_MT_db/na_MT_ensGene.txt",
  "sturgeon" = "/PATH/TO/biosoft/annovar_2020/sturgeon_MT_db/as_CNCB_MT_ensGene.txt",
  "stingray" = "/PATH/TO/biosoft/annovar_2020/stingray_MT_db/pl_MT_ensGene.txt"
)



vaf_results <- list()
for(sp in species_list){
  message("=== START ", sp," ===")
  annovar_file <- annovar_annotation_file_list[[sp]]
  genepred_file <- gene_pred_file_list[[sp]]
  vaf_file <- paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/",sp,"_detail_filtered_df_260128.tsv")
  
  res <- plot_vaf_spectrum_by_type(sp,
                                   annovar_file,
                                   genepred_file,
                                   vaf_file,
                                   "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/plots/VAF_vs_MutationType_260128")
}

# cd /PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/plots/VAF_vs_MutationType_260128/
# pdfunite *.pdf AllSpecies_VAF_spectrum_by_type.pdf
