setwd("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/")
library(SummarizedExperiment)
library(Seurat)
library(SeuratObject)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(magrittr)
library(purrr)
library(tibble)
library(stringr)

mgatk.out.path <- c(
  "human" = "/PATH/TO/pbmc.mt/SNP/human/mgatk_newCleanBam/final/",
  "mouse" = "/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/",
  "rat" = "/PATH/TO/pbmc.mt/SNP/rat/mgatk_newCleanBam/final/",
  "chicken" = "/PATH/TO/pbmc.mt/SNP/chicken/mgatk_newCleanBam/final/",
  "cobra" = "/PATH/TO/pbmc.mt/SNP/cobra/mgatk_newCleanBam/final/",
  "sturgeon" = "/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk_newCleanBam/final/",
  "stingray" = "/PATH/TO/pbmc.mt/SNP/stingray/mgatk_newCleanBam/final/"
)

run_mgatk_analysis <- function(species) {
  if (!species %in% names(mgatk.out.path)) {
    stop("Species is not in the predefined list, pleaseCheckingspeciesname: ", species)
  }
  
  prepare_data <- function(species) {
    message("Step 1: prepare data...")
    
    base_path <- mgatk.out.path[species]
    rds_path <- file.path(base_path, paste0(species, ".mgatk.rds"))
    downstream_dir <- file.path(dirname(base_path), "downstream")
    
    if (!dir.exists(downstream_dir)) {
      dir.create(downstream_dir, recursive = TRUE)
    }
    
    mgatk_se <- readRDS(rds_path)
    
    objlist <- readRDS("/PATH/TO/pbmc.mt/celltype.anno/seurat.meta/objlist.RDS")
    species.seurat <- objlist[[species]]
    
    mgatk_se <- mgatk_se[, colnames(species.seurat)]
    
    return(list(mgatk_se = mgatk_se, 
                species.seurat = species.seurat,
                downstream_dir = downstream_dir))
  }
  
  detect_mutations <- function(mgatk_se, downstream_dir) {
    message("Step 2: detect mutations...")
    
    source("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/variant_calling.R")
    
    mut_se <- call_mutations_mgatk(mgatk_se)
    
    saveRDS(mut_se, file.path(downstream_dir, "mut_se.RDS"))
    
    misc_df <- data.frame(rowData(mut_se))
    filter_df <- misc_df %>% 
      dplyr::filter(n_cells_conf_detected >= 1 & strand_correlation >= 0.65)
    
    return(list(mut_se = mut_se, filter_df = filter_df))
  }
  
  visualize_umap <- function(mut_se, filter_df, species.seurat, downstream_dir, species) {
    message("Step 3: plot UMAP...")
    
    mito_df <- data.frame(
      species.seurat@meta.data,
      species.seurat@reductions$umap@cell.embeddings,
      t(assays(mut_se)[["allele_frequency"]][as.character(filter_df$variant), ]),
      t(assays(mut_se)[["coverage"]][as.character(filter_df$variant), ])
    )
    
    plot_mito_umap <- function(mito_df, variant) {
      safe_colname <- make.names(variant)
      mito_df$variant_value <- mito_df[, safe_colname] * 100
      
      ggplot(mito_df %>% arrange(variant_value), 
             aes(x = umap_1, y = umap_2, color = variant_value)) +
        geom_point(size = 0.5) + 
        scale_color_gradientn(
          colors = c("lightgrey", "firebrick"),
          limits = c(0, 100)) +
        ggtitle(paste(variant)) +
        labs(color = "Frequency (%)") +
        theme_bw() +
        theme(
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5, color = "black")
        )
    }
    
    umap_dir <- file.path(downstream_dir, "umap_plots")
    if (!dir.exists(umap_dir)) {
      dir.create(umap_dir, recursive = TRUE)
    }
    message("UMAPPlots will be saved to: ", umap_dir)
    
    for (variant in filter_df$variant) {
      safe_variant <- gsub(">", "_", variant)
      safe_variant <- gsub(":", "-", safe_variant)
      
      pdf_file <- file.path(umap_dir, paste0(species, "_", safe_variant, ".umap.pdf"))
      
      pdf(pdf_file, width = 6, height = 5)
      print(plot_mito_umap(mito_df, variant))
      dev.off()
      
      message("Saved: ", basename(pdf_file))
    }
    
    return(umap_dir)
  }
  
  analyze_mutation_frequency <- function(species, downstream_dir) {
    message("Step 4: analyze mutation frequency...")
    
    base_path <- mgatk.out.path[species]
    hetero_path <- file.path(base_path, paste0(species, ".mgatk.cell_heteroplasmic_df.tsv.gz"))
    stats_path <- file.path(base_path, paste0(species, ".mgatk.variant_stats.tsv.gz"))
    
    hetero_df <- data.table::fread(hetero_path, sep = "\t", header = TRUE)
    stats_df <- data.table::fread(stats_path, sep = "\t", header = TRUE)
    
    data.table::setnames(hetero_df, names(hetero_df)[1], "barcode")
    data.table::setkey(hetero_df, "barcode")
    
    stats_df <- stats_df %>% 
      dplyr::filter(n_cells_conf_detected >= 1 & strand_correlation >= 0.65)
    
    bins <- seq(0, 1, by = 0.1)
    bin_labels <- sprintf("(%.1f-%.1f]", bins[-length(bins)], bins[-1])
    
    
    mut.freq.stat <- hetero_df %>%
      tidyr::pivot_longer(
        cols = -barcode,
        names_to = "mutation",
        values_to = "frequency"
      ) %>%
      dplyr::filter(mutation %in% stats_df$variant) %>%
      dplyr::mutate(frequency = ifelse(is.na(frequency), 0, frequency))
    
    mut.freq.stat$bin <- cut(
      mut.freq.stat$frequency,
      breaks = bins,
      include.lowest = FALSE,
      right = TRUE,
      labels = bin_labels
    )
    mut.freq.stat <- mut.freq.stat %>% 
      dplyr::mutate(bin = ifelse(is.na(bin), "0", as.character(bin)))
    
    hist_data <- mut.freq.stat %>%
      dplyr::group_by(mutation, bin) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop_last") %>%
      dplyr::mutate(percentage = (count / sum(count)) * 100) %>%
      dplyr::ungroup() %>%
      tidyr::complete(mutation, bin, fill = list(count = 0, percentage = 0))
    hist_data$bin <- factor(hist_data$bin, levels = c("0",bin_labels))
    saveRDS(hist_data, file.path(downstream_dir, "per.mutation.freq.distribution.RDS"))
    
    return(hist_data)
  }
  
  plot_mutation_histograms <- function(hist_data, downstream_dir) {
    message("Step 5: plot mutation-frequency histogram...")
    
    plot_single_histogram <- function(mutation_name, data) {
      plot_data <- data %>% 
        dplyr::filter(mutation == mutation_name)
      
      ggplot(plot_data, aes(x = bin, y = percentage)) +
        geom_col(fill = "steelblue", alpha = 0.8) +
        geom_text(
          aes(label = sprintf("%.2f%%", percentage)),
          vjust = -0.5, size = 3, color = "black"
        ) +
        labs(
          x = "Mutation Frequency in Each Cell",
          y = "Percentage of Cells (%)",
          title = paste("Mutation:", mutation_name)
        ) +
        theme_bw() +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.2)))
    }
    
    mutations <- unique(hist_data$mutation)
    
    pdf(file.path(downstream_dir, "all_mutations_histograms.pdf"), width = 8, height = 6)
    for (mut in mutations) {
      print(plot_single_histogram(mut, hist_data))
    }
    dev.off()
  }
  
  
  plot_low_vmr_scatter <- function(species, downstream_dir) {
    message("Plotting low-VMR variant scatter plot: ", species)
    
    base_path <- mgatk.out.path[species]
    stats_path <- file.path(base_path, paste0(species, ".mgatk.variant_stats.tsv.gz"))
    
    stats_df <- data.table::fread(stats_path, sep = "\t", header = TRUE)
    
    stats_df <- stats_df %>% 
      dplyr::filter(n_cells_conf_detected >= 1) %>%
      dplyr::mutate(
        log_vmr = log10(vmr),
        low_vmr = ifelse(
          log_vmr < -2 & strand_correlation > 0.65, 
          "Low VMR", 
          "Others"
        )
      )
    
    p <- ggplot(stats_df, aes(x = strand_correlation, y = log_vmr)) +
      geom_point(aes(color = low_vmr), alpha = 0.6, size = 2) +
      geom_label_repel(
        data = subset(stats_df, low_vmr == "Low VMR"),
        aes(label = variant),
        color = "red",
        box.padding = 0.5,
        segment.color = "red",
        size = 3,
        max.overlaps = 1000
      ) +
      geom_hline(yintercept = -1, linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0.65, linetype = "dashed", color = "black") +
      scale_color_manual(values = c("Low VMR" = "red", "Others" = "grey50")) +
      labs(
        x = "Strand Correlation",
        y = expression(log[10](VMR)),
        title = paste(species, "Possible Germline Variants"),
        color = "Variant Type"
      ) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.title = element_text(color = "black", face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
        legend.position = "top"
      )
    
    pdf_file <- file.path(downstream_dir, "Scatter.low.VMR.pdf")
    ggsave(pdf_file, plot = p, width = 8, height = 7)
    
    message("Scatter plot saved to: ", pdf_file)
  }
  
  message("\nStarting species analysis: ", species, "\n")
  
  data_prep <- prepare_data(species)
  message("Done Step 1: prepare data")
  
  mutation_results <- detect_mutations(
    data_prep$mgatk_se, 
    data_prep$downstream_dir
  )
  message("Done step2: variant calling")
  umap_dir <- visualize_umap(
    mutation_results$mut_se,
    mutation_results$filter_df,
    data_prep$species.seurat,
    data_prep$downstream_dir,
    species
  )

  message("All UMAP plots saved to: ", umap_dir)
  message("Done step3: UMAPvisualization")
  hist_data <- analyze_mutation_frequency(
    species, 
    data_prep$downstream_dir
  )
  message("Done step4: mutation-frequency analysis")
  plot_mutation_histograms(
    hist_data, 
    data_prep$downstream_dir
  )
  message("Done step5: plot histogram")
  plot_low_vmr_scatter(species, data_prep$downstream_dir)
  message("Done Step 6: low-VMR scatter plot")
  message("\n", species, " analysis complete! Results saved to: ", data_prep$downstream_dir, "\n")
}

species_to_analyze <- c(#"human", 
			"mouse", 
			"rat", 
			"chicken", 
			"cobra", 
			"sturgeon", 
			"stingray"
			)

for (sp in species_to_analyze) {
  run_mgatk_analysis(sp)
}
