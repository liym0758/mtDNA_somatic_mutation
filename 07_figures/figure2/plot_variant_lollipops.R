setwd("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/lolipop_chart/mutation_cells/new_human/")


plot_mt_lollipop <- function(
    species_name = NULL,
    chr_name,
    vaf_rds_file = NULL,#data from histogram
    cells_include = NULL,
    variants_include = NULL,
    ref_fasta,
    variant_stat_dir = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter",
    gene_pred,
    p_cutoff = 0,
    lollipop_colors = c("#E9F3E9","#C8E0C9","#A9CFAB","#7AAE7C","#527A54"),
    cds_color   = "#217849",
    trna_color  = "#D9C4DD",
    rrna_color  = "#5386A3",
    dloop_color = "#BF5D58"
    ) {
      
      library(data.table)
      library(dplyr)
      library(stringr)
      library(ggplot2)
      library(rtracklayer)
      library(GenomicRanges)
      library(Biostrings)
      library(patchwork)
      library(ggrepel)
      
      # =========================
      # Part 1: Mean heteroplasmy
      # =========================
      
      get_pos <- function(x) as.integer(str_extract(x, "^\\d+"))
      
      if (!is.null(vaf_rds_file)) {
        
        vaf <- readRDS(vaf_rds_file)[[1]]
        
        if (!is.null(cells_include))
          vaf <- vaf[vaf$barcode %in% cells_include, ]
        
        if (!is.null(variants_include))
          vaf <- vaf[vaf$mutation %in% variants_include, ]
        
        df <- vaf %>%
          group_by(mutation) %>%
          dplyr::summarise(
            mean_freq = mean(frequency, na.rm = TRUE),
            n_cells = sum(frequency > 0, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          dplyr::mutate(
            position = get_pos(mutation),
            mut_label = mutation
          ) %>%
          dplyr::filter(mean_freq > p_cutoff)
        
      } else {
        
        dt <- fread(mgatk_file)
        mat <- as.matrix(dt[, -1])
        rownames(mat) <- dt[[1]]
        
        if (!is.null(cells_include))
          mat <- mat[intersect(cells_include, rownames(mat)), , drop = FALSE]
        
        if (!is.null(variants_include))
          mat <- mat[, intersect(variants_include, colnames(mat)), drop = FALSE]
        
        avg_freq <- apply(mat, 2, function(x) mean(x[x > 0], na.rm = TRUE))
        n_cells  <- apply(mat, 2, function(x) sum(x > 0, na.rm = TRUE))
        
        df <- data.frame(
          position  = get_pos(colnames(mat)),
          mean_freq = avg_freq,
          n_cells   = n_cells,
          mut_label = colnames(mat)
        ) %>% dplyr::filter(mean_freq > p_cutoff)
      }
      df <- df %>%
        dplyr::mutate(
          n_cells_safe = pmax(n_cells, 1L),
          size_bin = factor(
            pmin(3L, floor(log10(n_cells_safe))),
            levels = 0:3
          )
        )
      
      # =========================
      # Add pct_conf_detected_cell
      # =========================
      
      pct_file <- file.path(
        variant_stat_dir,
        species_name,
        "variant_stats_final_filtered.tsv"
      )
      
      pct_df <- fread(pct_file) %>%
        dplyr::select(
          variant,
          pct_conf_detected_cell
        )
      
      df <- df %>%
        dplyr::left_join(
          pct_df,
          by = c("mut_label" = "variant")
        ) %>%
        arrange(pct_conf_detected_cell)
      
      
      # ---------- main lollipop ----------
      p_main <- ggplot(df, aes(position, mean_freq)) +
        geom_segment(aes(xend = position, yend = 0),
                     color = "grey60", linewidth = 0.4) +
        geom_point(aes(size = size_bin, color = pct_conf_detected_cell), alpha = 0.9) +
        scale_color_gradientn(
          colors = lollipop_colors,
          limits = c(0.03, 0.12),
          breaks = c(0.03, 0.06, 0.09, 0.12),
          oob = scales::squish,
          labels = scales::percent_format(accuracy = 1)
        )+
        scale_size_manual(
          values = c("0"=0.5, "1"=1, "2"=2.0,"3"=4.0, "4"=5.5),
          breaks = as.character(0:4),
          drop = FALSE
        )+
      
        labs(
          y = "Mean heteroplasmy",
          x = "Position",
          color = "Confidently\ndetected cells%",
          size = "log10(n cells)",
          title = str_to_title(paste0(species_name))
        ) +
        theme_classic(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))+
        coord_cartesian(ylim = c(0, 0.7))
      
      # =========================
      # Part 2: Gene annotation
      # =========================
      mt_seq <- readDNAStringSet(ref_fasta)
      mt_len <- width(mt_seq)[1]
      
      
      read_genepred_mt <- function(genepred_file, chr_name, mt_len) {
        
        library(data.table)
        library(dplyr)
        
        gp <- fread(genepred_file, header = FALSE)
        
        setnames(
          gp,
          c(
            "bin",
            "tx_id",
            "chrom",
            "strand",
            "txStart",
            "txEnd",
            "cdsStart",
            "cdsEnd",
            "exonCount",
            "exonStarts",
            "exonEnds",
            "score",
            "gene_name",
            "cdsStartStat",
            "cdsEndStat",
            "exonFrames"
          )
        )
        
        gp <- gp[chrom == chr_name]
        
        df_gene <- gp %>%
          transmute(
            start = txStart + 1,
            end   = txEnd,
            gene_name = gene_name,
            strand = strand,
            gene_type = case_when(
              grepl("^MT-T|^T[A-Z]|^trn", gene_name, ignore.case = TRUE) ~ "tRNA",
              grepl("^MT-RNR|16S|12S|^Rnr|^rrn", gene_name, ignore.case = TRUE) ~ "rRNA",
              TRUE ~ "CDS"
            )
          ) %>%
          arrange(start)
        
        gene_start_min <- min(df_gene$start)
        gene_end_max   <- max(df_gene$end)
        
        dloop_df <- bind_rows(
          data.frame(
            start = 1,
            end = gene_start_min - 1,
            gene_name = "D-loop",
            strand = "*",
            gene_type = "Dloop"
          ),
          data.frame(
            start = gene_end_max + 1,
            end = mt_len,
            gene_name = "D-loop",
            strand = "*",
            gene_type = "Dloop"
          )
        ) %>% filter(start <= end)
        
        # ---------- Intergenic: gaps BETWEEN genes ----------
        intergenic_df <- df_gene %>%
          arrange(start) %>%
          mutate(
            next_start = lead(start),
            gap_start  = end + 1,
            gap_end    = next_start - 1
          ) %>%
          filter(!is.na(next_start), gap_start <= gap_end) %>%
          transmute(
            start = gap_start,
            end   = gap_end,
            gene_name = "Intergenic",
            strand = "*",
            gene_type = "Intergenic"
          )
        
        bind_rows(
          df_gene,
          intergenic_df,
          dloop_df
        ) %>% arrange(start)
        

        
        
      }
      
      
      
      
      df_gene <- read_genepred_mt(
        genepred_file = gene_pred,
        chr_name = chr_name,
        mt_len = mt_len
      )
      
      df_gene$gene_type <- factor(
        df_gene$gene_type,
        levels = c("CDS", "tRNA", "rRNA", "Dloop", "Intergenic")
      )
      
      fill_colors <- c(
        CDS = cds_color,
        tRNA = trna_color,
        rRNA = rrna_color,
        Dloop = dloop_color,
        Intergenic = "grey80"
      )
      
      p_gene <- ggplot(df_gene) +
        geom_rect(
          aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = gene_type),
          color = "black", linewidth = 0.25
        ) +
        labs(fill = "Region",
             x = "Position")+
        scale_fill_manual(values = fill_colors) +
        scale_x_continuous(limits = c(1, mt_len), expand = c(0,0)) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom"
        ) +
        ylim(0,1.4)
      
      # ---------- combine ----------
      p_main <- p_main +
        scale_x_continuous(limits = c(1, mt_len), expand = c(0,0))
      
      p_main / p_gene + plot_layout(heights = c(4,1))
    }

chr_name <- c("human" = "chrM",
              "rat" = "MT",
              "mouse" = "MT",
              "chicken" = "MT",
              "cobra" = "MT",
              "sturgeon" = "MT",
              "stingray" = "MT")

base_path <- "/PATH/TO/pbmc.mt/SNP/"
species_list <- c("human","mouse","rat","chicken","cobra","sturgeon","stingray")
ref_fasta <- list(
  "human" = "/PATH/TO/pbmc.mt/SNP/human_mix/mgatk/downstream/rm.germline/human.germline.rm.MT.genome.fa",
  "mouse" = "/PATH/TO/pbmc.mt/ref/mouse/new/mm39_mito.fa",
  "rat" = "/PATH/TO/pbmc.mt/SNP/rat/mgatk/downstream/rm.germline/rat.germline.rm.MT.genome.fa",
  "chicken" = "/PATH/TO/pbmc.mt/SNP/chicken/mgatk/downstream/rm.germline/chicken.germline.rm.MT.genome.fa",
  "cobra" = "/PATH/TO/pbmc.mt/SNP/cobra/mgatk/downstream/rm.germline/cobra.germline.rm.MT.genome.fa",
  "sturgeon" = "/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk/downstream/rm.germline/sturgeon.germline.rm.MT.genome.fa",
  "stingray" = "/PATH/TO/pbmc.mt/SNP/stingray/mgatk/downstream/rm.germline/stingray.germline.rm.MT.genome.fa"
)

# ref_gtf <- list(
#   "human" = "/PATH/TO/pbmc.mt/ref/human/mt.mask/Homo_Shapien_GRCh38/genes/MT_genes.gtf",
#   "mouse" = "/PATH/TO/pbmc.mt/ref/mouse/new/genome.mask/Mus_musculus_toplevel/genes/MT_genes.gtf.gz",
#   "rat" = "/PATH/TO/biosoft/annovar_2020/rat_MT_db/MT_genes.gtf",
#   "chicken" = "/PATH/TO/biosoft/annovar_2020/chicken_MT_db/MT_genes.gtf",
#   "cobra" = "/PATH/TO/biosoft/annovar_2020/cobra_MT_db/Naja_atra_mt.gtf",
#   "sturgeon" = "/PATH/TO/biosoft/annovar_2020/sturgeon_MT_db/CNCBdb_fixed_Acipenser_sinensis_MT.gtf",
#   "stingray" = "/PATH/TO/biosoft/annovar_2020/stingray_MT_db/MT_genes.gtf"
# )
gene_pred_file_list <- c(
  "human" = "/PATH/TO/biosoft/annovar_2020/humandb/hg38_chrM_ensGene.txt",
  "mouse" = "/PATH/TO/biosoft/annovar_2020/mm39_MT_db/mm39_MT_ensGene.txt",
  "rat" = "/PATH/TO/biosoft/annovar_2020/rat_MT_db/rn7_MT_ensGene.txt",
  "chicken" = "/PATH/TO/biosoft/annovar_2020/chicken_MT_db/gg1_MT_ensGene.txt",
  "cobra" = "/PATH/TO/biosoft/annovar_2020/cobra_MT_db/na_MT_ensGene.txt",
  "sturgeon" = "/PATH/TO/biosoft/annovar_2020/sturgeon_MT_db/as_CNCB_MT_ensGene.txt",
  "stingray" = "/PATH/TO/biosoft/annovar_2020/stingray_MT_db/pl_MT_ensGene.txt"
)

variant.list <- readRDS("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Final_Variant_List_7species_260128.RDS")


p <- list()

for (sp in species_list){
  p[[sp]] <- plot_mt_lollipop(
    species_name = sp,
    chr_name = chr_name[sp],
    #mgatk_file = heteroplasmic_df_path[spcies], 
    vaf_rds_file = paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Reorder/",sp,"_Reorder_celltype_VAF_data_260128.RDS"),
    variant_stat_dir = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter",
    variants_include = variant.list[[sp]], 
    ref_fasta = ref_fasta[[sp]],    
    gene_pred  = gene_pred_file_list[[sp]],   
    cds_color   = "#CCE0F5",
    trna_color  = "#3361A5",
    rrna_color  = "#7FB0B4",
    dloop_color = "#E5DBB9",
    p_cutoff = 0,       # min mean frequency threshold
    lollipop_colors = c('#A8DDE1','#75B5DC','#478ECC','#326DB6','#313772')
  )
  
  
  
  write.table(p[[sp]][[1]]$data,
              file = paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/lolipop_chart/mutation_cells/new_human/",sp,"_lolipop_SameYscale_DotScale_final_filter_data_260128.tsv"),
              sep = "\t",quote = F
  )
}


pdf(file = paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/lolipop_chart/mutation_cells/new_human/7species_lolipop_SameYscale_DotScale_final_filter_260128.pdf"), width = 10, height = 4)
for (sp in species_list){
  print(p[[sp]])
}
dev.off()

