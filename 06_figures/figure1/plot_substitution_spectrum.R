setwd("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/SNP_spectrum_newCleanBam/new_human_simple_spectrum/")
library(data.table)
library(dplyr)
library(ggplot2)
library(stringi)

# 1. SNP Spectrum ------------------------------------
plot_snp_spectrum_simple <- function(mgatk_path, 
                                     ref_allele_txt = "MT_refAllele.txt",
                                     sample_name, 
                                     called_variants) {
  ref_all <- fread(file.path(mgatk_path, ref_allele_txt))
  colnames(ref_all) <- c("pos", "ref")
  ref_all$ref <- toupper(ref_all$ref)
  
  ref_all_long <- rbind(ref_all, ref_all, ref_all, ref_all)
  ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = nrow(ref_all))
  ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt, ]
  
  # variant id
  ref_all_long$variant <- paste0(ref_all_long$pos, ref_all_long$ref, ">", ref_all_long$alt)
  ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
  
  reverse_complement <- function(s) {
    chartr("ATGC", "TACG", s)
  }
  
  ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C", "T"), "L", "H")
  
  ref_all_long$change_rc <- reverse_complement(ref_all_long$change)
  ref_all_long$group_change <- ifelse(ref_all_long$ref %in% c("C","T"),
                                      ref_all_long$change,
                                      ref_all_long$change_rc)
  
  ref_all_long$called <- ref_all_long$variant %in% called_variants
  
  total <- nrow(ref_all_long)
  total_called <- sum(ref_all_long$called)
  
  prop_df <- ref_all_long %>%
    group_by(group_change, strand) %>%
    dplyr::summarize(
      observed_prop_called = sum(called) / total_called,
      expected_prop = dplyr::n() / total,
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(fc_called = observed_prop_called / expected_prop)
  
  facet_order <- c("CA" = "C>A", "CG" = "C>G", "CT" = "C>T",
                   "TA" = "T>A", "TC" = "T>C", "TG" = "T>G")
  prop_df <- prop_df %>%
    dplyr::filter(group_change %in% names(facet_order))
  
  prop_df$group_change <- factor(prop_df$group_change,
                                 levels = names(facet_order))
  
  # Fold change plot
  p1 <- ggplot(prop_df, aes(x = strand, fill = strand, y = fc_called)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "black") +
    facet_wrap(~ group_change, nrow = 1,
               labeller = labeller(group_change = facet_order),
               strip.position = "bottom") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c("H" = "#A42227", "L" = "#003E7A")) +
    theme_classic(base_size = 12) +
    labs(y = "Substitution Rate \n (Observed / Expected)",
         x = NULL,
         title = paste(stri_trans_totitle(sample_name), "SNP Spectrum")) +
    geom_hline(yintercept = 1, linetype = 2, color = "black") +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(color = "black", face = "bold", size = 12),
      axis.title.x = element_blank(),  
      axis.text.x = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      strip.placement = "outside"
    )
  
  # Observed frequency plot
  p2 <- ggplot(prop_df, aes(x = strand, fill = strand, y = observed_prop_called)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8, color = "black") +
    facet_wrap(~ group_change, nrow = 1,
               labeller = labeller(group_change = facet_order),
               strip.position = "bottom") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = c("H" = "#A42227", "L" = "#003E7A")) +
    theme_classic(base_size = 12) +
    labs(y = "Observed Substitution Frequency",
         x = NULL,
         title = paste(stri_trans_totitle(sample_name), "SNP Spectrum")) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(color = "black", face = "bold", size = 12),
      axis.title.x = element_blank(), axis.text.x = element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
      strip.background = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      strip.placement = "outside"
    )
  
  return(list(fc_plot = p1, observed_plot = p2))
}

mgatk.out.path <- c(
  "human"="/PATH/TO/pbmc.mt/SNP/human/mgatk_newCleanBam/merge_pbmc_samples/final/",
  "mouse"="/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/",
  "rat"="/PATH/TO/pbmc.mt/SNP/rat/mgatk_newCleanBam/final/",
  "chicken"="/PATH/TO/pbmc.mt/SNP/chicken/mgatk_newCleanBam/final/",
  "cobra"="/PATH/TO/pbmc.mt/SNP/cobra/mgatk_newCleanBam/final/",
  "sturgeon"="/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk_newCleanBam/final/",
  "stingray"="/PATH/TO/pbmc.mt/SNP/stingray/mgatk_newCleanBam/final/"
)


variant.list.final <- readRDS("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Final_Variant_List_7species_260128.RDS")
#data <- read_tsv("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/human/variant_stats_final_filtered.tsv")
spectrum <- list()

i <- 1
spectrum[[i]] <- plot_snp_spectrum_simple(mgatk.out.path[i], 
                                          ref_allele_txt = "chrM_refAllele.txt",
                                          species[i], 
                                          variant.list.final[[i]]
                                          #unique(data$variant)
)



for (i in 2:length(species)){
  # data <- read_tsv(
  #   paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/",species[i],"/variant_stats_final_filtered.tsv")
  # )
  spectrum[[i]] <- plot_snp_spectrum_simple(mgatk.out.path[i], 
                                            ref_allele_txt = "MT_refAllele.txt",
                                            species[i], 
                                            variant.list.final[[i]]
                                            #unique(data$variant)
  )
}

ymax_fc <- max(sapply(spectrum, function(x) {
  max(ggplot_build(x[[1]])$data[[1]]$y, na.rm = TRUE)
}))
# > ymax_fc
# [1] 12.19683

ymax_obs <- max(sapply(spectrum, function(x) {
  max(ggplot_build(x[[2]])$data[[1]]$y, na.rm = TRUE)
}))
# > ymax_obs
# [1] 0.5336788

pdf("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/SNP_spectrum_newCleanBam/new_human_simple_spectrum/FCandRaw_simple_spectrum_SameYaxis_260128.pdf", width = 11, height = 3)
for (i in 1:length(species)){
  plot(
    (spectrum[[i]][[1]] + scale_y_continuous(limits = c(0, 12),
                                             expand = expansion(mult = c(0, 0.05)))) |
      (spectrum[[i]][[2]] + scale_y_continuous(limits = c(0, 0.55),
       expand = expansion(mult = c(0, 0.05))))
  )
}
dev.off()

pdf("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/SNP_spectrum_newCleanBam/new_human_simple_spectrum/FC_simple_spectrum_SameYaxis_260128.pdf", width = 6, height = 3)
for (i in 1:length(species)){
  plot(spectrum[[i]][[1]]+ scale_y_continuous(limits = c(0, 12),
                                              expand = expansion(mult = c(0, 0.05))))
}
dev.off()



# 2. Pie Chart ----------------------
plot_snp_pie_chart <- function(mgatk_path, 
                               ref_allele_txt = "MT_refAllele.txt",
                               sample_name, 
                               called_variants) {
  
  ref_all <- fread(file.path(mgatk_path, ref_allele_txt))
  colnames(ref_all) <- c("pos", "ref")
  ref_all$ref <- toupper(ref_all$ref)
  
  ref_all_long <- rbind(ref_all, ref_all, ref_all, ref_all)
  ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = nrow(ref_all))
  ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt, ]
  
  # variant id
  ref_all_long$variant <- paste0(ref_all_long$pos, ref_all_long$ref, ">", ref_all_long$alt)
  ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
  
  reverse_complement <- function(s) {
    chartr("ATGC", "TACG", s)
  }
  
  ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C", "T"), "L", "H")
  
  ref_all_long$change_rc <- reverse_complement(ref_all_long$change)
  ref_all_long$group_change <- ifelse(ref_all_long$ref %in% c("C","T"),
                                      ref_all_long$change,
                                      ref_all_long$change_rc)
  
  ref_all_long$called <- ref_all_long$variant %in% called_variants
  
  ref_all_long$complement_label <- sapply(ref_all_long$change, function(x) {
    ref <- substr(x, 1, 1)
    alt <- substr(x, 2, 2)
    comp_ref <- reverse_complement(ref)
    comp_alt <- reverse_complement(alt)
    paste0(ref, ">", alt, "/", comp_ref, ">", comp_alt)
  })
  
  ref_all_long$group <- "Others"
  ref_all_long$group[ref_all_long$complement_label %in% c("T>C/A>G", "A>G/T>C", "C>T/G>A", "G>A/C>T")] <- "Spontaneous deamination"
  ref_all_long$group[ref_all_long$complement_label %in% c("G>T/C>A", "C>A/G>T", "G>C/C>G", "C>G/G>C")] <- "ROS damage"
  
  total <- nrow(ref_all_long)
  total_called <- sum(ref_all_long$called)
  
  prop_df <- ref_all_long %>%
    group_by(group_change, strand) %>%
    dplyr::summarize(
      observed_prop_called = sum(called) / total_called,
      expected_prop = dplyr::n() / total,
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(fc_called = observed_prop_called / expected_prop)
  
  facet_order <- c("CA" = "C>A", "CG" = "C>G", "CT" = "C>T",
                   "TA" = "T>A", "TC" = "T>C", "TG" = "T>G")
  
  prop_df <- prop_df %>%
    dplyr::filter(group_change %in% names(facet_order))
  
  prop_df$group_change <- factor(prop_df$group_change,
                                 levels = names(facet_order))
  
  group_mapping <- data.frame(
    group_change = c("CA", "CG", "CT", "TA", "TC", "TG"),
    complement_label = c("C>A/G>T", "C>G/G>C", "C>T/G>A", 
                         "T>A/A>T", "T>C/A>G", "T>G/A>C"),
    group = c("ROS damage", "ROS damage", "Spontaneous deamination",
              "Others", "Spontaneous deamination", "Others")
  )
  
  prop_df <- prop_df %>%
    left_join(group_mapping, by = "group_change")
  
  group_fc <- prop_df %>%
    group_by(group) %>%
    dplyr::summarize(
      total_fc = sum(fc_called),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      corrected_prop = total_fc / sum(total_fc)
    )
  
  cat("=== Corrected pie-chart data:", sample_name, "===\n")
  print(group_fc)
  
  group_colors <- c("Spontaneous deamination" = "#AF95C4", 
                    "ROS damage" = "#5F3D69", 
                    "Others" = "#E3E3E3")
  group_fc$group <- factor(group_fc$group, levels = c("Spontaneous deamination", "ROS damage", "Others"))
  
  p <- ggplot(group_fc, aes(x = "", y = corrected_prop, fill = group)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = group_colors) +
    theme_void(base_size = 12) +
    labs(title = paste(stri_trans_totitle(sample_name)),
         fill = "Mutation Mechanisms") +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
  
  p <- p + geom_text(aes(label = paste0(round(corrected_prop * 100, 1), "%")),
                     position = position_stack(vjust = 0.5),
                     size = 4)
  
  return(p)
}
pie_charts <- list()
i <- 1
# data <- read_tsv("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/human/variant_stats_final_filtered.tsv")

pie_charts[[i]] <- plot_snp_pie_chart(mgatk.out.path[i],
                                      ref_allele_txt = "chrM_refAllele.txt",
                                      species[i],
                                      variant.list.final[[i]]
                                      #unique(data$variant)
                                      
                                      )

for (i in 2:length(species)){
  # data <- read_tsv(
  #   paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/",species[i],"/variant_stats_final_filtered.tsv")
  # )
  
  pie_charts[[i]] <- plot_snp_pie_chart(mgatk.out.path[i],
                                        ref_allele_txt = "MT_refAllele.txt",
                                        species[i],
                                        variant.list.final[[i]]
                                        #unique(data$variant)
                                        )
}

pdf("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/SNP_spectrum_newCleanBam/new_human_simple_spectrum/pie_charts_260128.pdf", width = 4, height = 3)
for (i in 1:length(species)){
  plot(pie_charts[[i]])
}
dev.off()
