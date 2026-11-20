setwd("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/SNP_spectrum_newCleanBam/new_human_simple_spectrum/")
library(data.table)
library(dplyr)
library(ggplot2)
library(stringi)
library(patchwork)

# 1. SNP Spectrum ------------------------------------
plot_snp_spectrum_simple <- function(mgatk_path, 
                                     ref_allele_txt = "MT_refAllele.txt",
                                     sample_name, 
                                     called_variants,
                                     low_coverage_positions = integer(0)) {
  ref_all <- fread(file.path(mgatk_path, ref_allele_txt))
  colnames(ref_all) <- c("pos", "ref")
  ref_all$ref <- toupper(ref_all$ref)
  ref_all$pos <- as.integer(ref_all$pos)
  # Expected and Observed must use exactly the same callable positions.
  low_coverage_positions <- unique(as.integer(low_coverage_positions))
  low_coverage_positions <- low_coverage_positions[!is.na(low_coverage_positions)]
  ref_all <- ref_all[!ref_all$pos %in% low_coverage_positions, ]
  if (nrow(ref_all) == 0L) stop("No callable reference positions for ", sample_name)
  called_variants <- unique(as.character(called_variants))
  called_positions <- as.integer(sub("^([0-9]+).*", "\\1", called_variants))
  called_variants <- called_variants[!called_positions %in% low_coverage_positions]
  
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
  
  return(list(fc_plot = p1, observed_plot = p2, prop_df = prop_df))
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

# Apply the very same per-species low-depth mask used for mutation filtering.
low.coverage.position <- readRDS(
  "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/coverage/results/new_human_data/low.coverage.position.7SpeciesTable.RDS"
)
species <- names(mgatk.out.path)
stopifnot(all(species %in% names(low.coverage.position)),
          all(species %in% names(variant.list.final)))
low.pos <- lapply(species, function(sp) {
  as.integer(unique(low.coverage.position[[sp]]$position))
})
names(low.pos) <- species

# Validate input, and record how much of each reference is callable.
audit <- lapply(species, function(sp) {
  ref_name <- if (sp == "human") "chrM_refAllele.txt" else "MT_refAllele.txt"
  r <- fread(file.path(mgatk.out.path[[sp]], ref_name), header = FALSE)
  names(r) <- c("pos", "ref")
  p <- as.integer(r$pos)
  bad <- intersect(p, low.pos[[sp]])
  vv <- variant.list.final[[sp]]
  vp <- as.integer(sub("^([0-9]+).*", "\\1", vv))
  overlap <- sum(vp %in% low.pos[[sp]])
  data.frame(species=sp,
             ref_positions=nrow(r),
             low_depth_rows=nrow(low.coverage.position[[sp]]),
             low_depth_unique_positions=length(bad),
             callable_positions=nrow(r)-length(bad),
             callable_fraction=(nrow(r)-length(bad))/nrow(r),
             retained_variants=length(vv),
             retained_variants_in_low_depth=overlap)
})
audit <- do.call(rbind, audit)
output_dir <- file.path(getwd(), "code", "callable_spectrum")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
fwrite(audit, file.path(output_dir, "callable_background_audit.tsv"), sep="\t")
print(audit)
if(any(audit$retained_variants_in_low_depth > 0)) {
  warning("Some retained variants are on low-coverage positions; they will be excluded from BOTH numerator and background. See audit TSV.")
}

#data <- read_tsv("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/human/variant_stats_final_filtered.tsv")
spectrum <- list()

i <- 1
spectrum[[i]] <- plot_snp_spectrum_simple(mgatk.out.path[i], 
                                          ref_allele_txt = "chrM_refAllele.txt",
                                          species[i], 
                                          variant.list.final[[species[i]]],
                                          low_coverage_positions = low.pos[[species[i]]]
                                          #unique(data$variant)
)



for (i in 2:length(species)){
  # data <- read_tsv(
  #   paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/",species[i],"/variant_stats_final_filtered.tsv")
  # )
  spectrum[[i]] <- plot_snp_spectrum_simple(mgatk.out.path[i], 
                                            ref_allele_txt = "MT_refAllele.txt",
                                            species[i], 
                                            variant.list.final[[species[i]]],
                                          low_coverage_positions = low.pos[[species[i]]]
                                            #unique(data$variant)
  )
}

# Audit the exact fold-enrichment changes in the 12 substitution/strand groups.
comparison_tables <- lapply(seq_along(species), function(i) {
  sp <- species[i]
  ref_name <- if (sp == "human") "chrM_refAllele.txt" else "MT_refAllele.txt"
  before <- plot_snp_spectrum_simple(mgatk.out.path[[sp]], ref_name, sp,
                                      variant.list.final[[sp]])$prop_df
  after <- spectrum[[i]]$prop_df
  comparison <- dplyr::full_join(
    dplyr::select(before, group_change, strand,
                  observed_before=observed_prop_called,
                  expected_before=expected_prop,
                  fold_before=fc_called),
    dplyr::select(after, group_change, strand,
                  observed_after=observed_prop_called,
                  expected_after=expected_prop,
                  fold_after=fc_called),
    by=c("group_change", "strand"))
  comparison$species <- sp
  comparison
})
comparison_tables <- dplyr::bind_rows(comparison_tables)
fwrite(comparison_tables,
       file.path(output_dir, "spectrum_before_vs_after.tsv"), sep="\t")


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

pdf(file.path(output_dir, "FCandRaw_simple_spectrum_SameYaxis_260128_callable.pdf"), width = 11, height = 3)
for (i in 1:length(species)){
  plot(
    (spectrum[[i]][[1]] + scale_y_continuous(limits = c(0, ceiling(ymax_fc * 1.05)),
                                             expand = expansion(mult = c(0, 0.05)))) |
      (spectrum[[i]][[2]] + scale_y_continuous(limits = c(0, ceiling(ymax_obs * 1.05 / 0.05) * 0.05),
       expand = expansion(mult = c(0, 0.05))))
  )
}
dev.off()

pdf(file.path(output_dir, "FC_simple_spectrum_SameYaxis_260128_callable.pdf"), width = 6, height = 3)
for (i in 1:length(species)){
  plot(spectrum[[i]][[1]]+ scale_y_continuous(limits = c(0, ceiling(ymax_fc * 1.05)),
                                              expand = expansion(mult = c(0, 0.05))))
}
dev.off()



# 2. Pie Chart ----------------------
plot_snp_pie_chart <- function(mgatk_path, 
                               ref_allele_txt = "MT_refAllele.txt",
                               sample_name, 
                               called_variants,
                               low_coverage_positions = integer(0)) {
  
  ref_all <- fread(file.path(mgatk_path, ref_allele_txt))
  colnames(ref_all) <- c("pos", "ref")
  ref_all$ref <- toupper(ref_all$ref)
  ref_all$pos <- as.integer(ref_all$pos)
  # Expected and Observed must use exactly the same callable positions.
  low_coverage_positions <- unique(as.integer(low_coverage_positions))
  low_coverage_positions <- low_coverage_positions[!is.na(low_coverage_positions)]
  ref_all <- ref_all[!ref_all$pos %in% low_coverage_positions, ]
  if (nrow(ref_all) == 0L) stop("No callable reference positions for ", sample_name)
  called_variants <- unique(as.character(called_variants))
  called_positions <- as.integer(sub("^([0-9]+).*", "\\1", called_variants))
  called_variants <- called_variants[!called_positions %in% low_coverage_positions]
  
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
                                      variant.list.final[[species[i]]],
                                      low_coverage_positions = low.pos[[species[i]]]
                                      #unique(data$variant)
                                      
                                      )

for (i in 2:length(species)){
  # data <- read_tsv(
  #   paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/",species[i],"/variant_stats_final_filtered.tsv")
  # )
  
  pie_charts[[i]] <- plot_snp_pie_chart(mgatk.out.path[i],
                                        ref_allele_txt = "MT_refAllele.txt",
                                        species[i],
                                        variant.list.final[[species[i]]],
                                      low_coverage_positions = low.pos[[species[i]]]
                                        #unique(data$variant)
                                        )
}

pdf(file.path(output_dir, "pie_charts_260128_callable.pdf"), width = 4, height = 3)
for (i in 1:length(species)){
  plot(pie_charts[[i]])
}
dev.off()


# Additional C>T/G>A comparison: raw 2-category counts and exposure ratios.
ct_audit <- do.call(rbind, lapply(species, function(sp) {
  rn <- if (sp == "human") "chrM_refAllele.txt" else "MT_refAllele.txt"
  rr <- fread(file.path(mgatk.out.path[[sp]], rn), header=FALSE)
  names(rr) <- c("pos","ref")
  rr$pos <- as.integer(rr$pos)
  rr$ref <- toupper(rr$ref)
  rrc <- rr[!rr$pos %in% low.pos[[sp]], ]
  vv <- variant.list.final[[sp]]
  vp <- as.integer(sub("^([0-9]+).*", "\\1", vv))
  vv <- vv[!vp %in% low.pos[[sp]]]
  nC <- sum(rrc$ref=="C")
  nG <- sum(rrc$ref=="G")
  nCT <- sum(grepl("^[0-9]+C>T$", vv))
  nGA <- sum(grepl("^[0-9]+G>A$", vv))
  data.frame(species=sp, nCT=nCT, nGA=nGA, callable_C=nC, callable_G=nG,
    L_H_rate_ratio=if(nGA>0 && nC>0) (nCT/nC)/(nGA/nG) else NA_real_)
}))
fwrite(ct_audit, file.path(output_dir, "CT_GA_callable_background.tsv"), sep="\t")
print(ct_audit)
