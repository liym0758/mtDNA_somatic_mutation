setwd("/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/strand.test/")
library(readr)
library(dplyr)
library(ggplot2)

mgatk_file <- "/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/mouse.mgatk.variant_stats.tsv.gz"

library(data.table)
library(ggplot2)
library(dplyr)
library(stringr)

# --- Step 1. Load mgatk sparse matrices for A, C, G, T and coverage ---
load_sparse_matrix <- function(file, base) {
  dt <- fread(cmd = paste("zcat", file), header = FALSE)
  colnames(dt) <- c("pos", "cell", "fwd", "rev")
  dt$allele <- base
  dt
}

# Define your file paths here
prefix <- "/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/mouse.mgatk"
files <- list(
  A = paste0(prefix, ".A.txt.gz"),
  C = paste0(prefix, ".C.txt.gz"),
  G = paste0(prefix, ".G.txt.gz"),
  T = paste0(prefix, ".T.txt.gz"),
  cov = paste0(prefix, ".coverage.txt.gz")
)

alleles <- lapply(names(files)[1:4], function(base) load_sparse_matrix(files[[base]], base))
allele_dt <- rbindlist(alleles)
allele_dt[, count := fwd + rev]
allele_dt$raw_strand_bias <- allele_dt$fwd/allele_dt$count


write.table(allele_dt, file="/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/cell_strand_position_count_data.tsv", sep = "\t")

strand_summary <- as.data.table(allele_dt)[, .(
  raw_fwd_total = sum(fwd,na.rm = TRUE),
  raw_rev_total = sum(rev,na.rm = TRUE),
  raw_count_total = sum((fwd+rev),na.rm = TRUE), 
  raw_mean_count_PerCell = mean((fwd+rev),na.rm = TRUE),
  raw_mean_bias = mean(raw_strand_bias)
), by = .(pos, allele)]


variant_stats <- fread("zcat /PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/mouse.mgatk.variant_stats.tsv.gz")
variant_stats <- variant_stats %>% dplyr::filter(n_cells_conf_detected >= 1)

variant_stats[, c("pos", "ref", "alt") := tstrsplit(gsub("([0-9]+)([ACGT])>([ACGT])", "\\1,\\2,\\3", variant), ",")]
variant_stats$pos <- as.integer(variant_stats$pos)

variant_stats <- merge(variant_stats, strand_summary, 
                       by.x = c("pos", "alt"), 
                       by.y = c("pos", "allele"),
                       all.x = F) %>% drop_na()


library(ggplot2)
library(data.table)

# python /PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/strand.test/binomial_CI.py

ci <- read_csv("/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/strand.test/binomial_CI_smooth.csv")

ci <- ci %>% dplyr::filter(depth <= (max(variant_stats$raw_mean_count_PerCell)+1))

fit_monotonic_spline <- function(x, y, increasing = TRUE) {
  xx <- seq(min(x), max(x), by = 1)
  fit <- smooth.spline(x, y, spar = 0.7)
  yy <- predict(fit, xx)$y
  
  if (increasing) {
    yy <- cummax(yy)
  } else {
    yy <- cummin(yy)
  }
  tibble(depth = xx, value = pmin(pmax(yy, 0), 1))
}

lower_95 <- fit_monotonic_spline(ci$depth, ci$lower_95, increasing = TRUE)
upper_95 <- fit_monotonic_spline(ci$depth, ci$upper_95, increasing = FALSE)
lower_99 <- fit_monotonic_spline(ci$depth, ci$lower_99, increasing = TRUE)
upper_99 <- fit_monotonic_spline(ci$depth, ci$upper_99, increasing = FALSE)

smooth_ci <- reduce(list(
  lower_95 %>% rename(lower_95 = value),
  upper_95 %>% rename(upper_95 = value),
  lower_99 %>% rename(lower_99 = value),
  upper_99 %>% rename(upper_99 = value)
), full_join, by = "depth")
smooth_ci[nrow(smooth_ci)+1,] <- list(1,0,1,0,1)
#variant_stats$strand_bias <- variant_stats$fwd_total/variant_stats$count_total

p2 <- ggplot(variant_stats, aes(x = raw_mean_count_PerCell, y = raw_mean_bias)) +
  geom_ribbon(data = smooth_ci, aes(x = depth, ymin = lower_99, ymax = upper_99),
              fill = "#6778AE", alpha = 0.3, inherit.aes = FALSE) +
  geom_ribbon(data = smooth_ci, aes(x = depth, ymin = lower_95, ymax = upper_95),
              fill = "#6778AE", alpha = 0.3, inherit.aes = FALSE) +
  
  geom_point(color = "black", size = 1, alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.2) +
  labs(title = "Raw Depth vs. Raw Bias (with Binomial CI)",
       x = "Raw Mean Count Per Variant Per Cell",
       y = "Raw Mean Bias",
       caption = "Blue bands: 95% and 99% confidence intervals based on binomial distribution") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        plot.caption = element_text(size = 8, color = "gray50"))


library(data.table)
library(dplyr)
library(ggplot2)
library(stringi)

plot_snp_spectrum <- function(mgatk_path, 
                              ref_allele_txt = "MT_refAllele.txt",
                              sample_name, 
                              called_variants) {
  ref_all <- fread(file = file.path(mgatk_path, ref_allele_txt))
  colnames(ref_all) <- c("pos", "ref")
  ref_all$ref <- toupper(ref_all$ref)
  l <- as.character(ref_all$ref)
  
  ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))
  ref_all <- ref_all[!grepl("N", ref_all$three),]
  
  ref_all_long <- rbind(ref_all, ref_all, ref_all, ref_all)
  ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = nrow(ref_all))
  ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt,]
  
  ref_all_long$variant <- paste0(ref_all_long$pos, ref_all_long$ref, ">", ref_all_long$alt)
  ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
  
  reverse_complement <- function(s) {
    chartr("ATGC", "TACG", s)
  }
  
  ref_all_long$change_rc <- reverse_complement(ref_all_long$change)
  ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C", "T"), "L", "H")
  ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
  ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
  ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)
  
  ref_all_long$called <- ref_all_long$variant %in% called_variants
  
  total <- nrow(ref_all_long)
  total_called <- sum(ref_all_long$called)
  prop_df <- ref_all_long %>%
    group_by(three_plot, group_change, strand) %>%
    summarize(
      observed_prop_called = sum(called) / total_called,
      expected_prop = n() / total,
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(fc_called = observed_prop_called / expected_prop)
  
  prop_df <- prop_df %>%
    mutate(
      change_plot = paste0(group_change, "_", three_plot),
      context_label = paste0(substr(three_plot, 1, 1), "_", substr(three_plot, 3, 3)),
      unique_id = paste0(group_change, "_", three_plot)
    )
  
  facet_order <- c("CA" = "C>A", "CG" = "C>G", "CT" = "C>T",
                   "TA" = "T>A", "TC" = "T>C", "TG" = "T>G")
  
  # Fold change plot
  p1 <- ggplot(prop_df, aes(x = unique_id, fill = strand, y = fc_called)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ group_change, nrow = 1, labeller = labeller(group_change = facet_order), scales = "free_x") +
    scale_x_discrete(
      name = "Trinucleotide Context",
      labels = function(x) {
        sapply(strsplit(x, "_"), function(s) paste0(substr(s[2], 1, 1), "_", substr(s[2], 3, 3)))
      }
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      limits = c(0, max(prop_df$fc_called) * 1.05)
    ) +
    scale_fill_manual(values = c("#E08D7A", "#859DCD")) +
    theme_bw(base_size = 12) +
    labs(y = "Substitution Rate (Observed / Expected)",
         title = paste(stri_trans_totitle(sample_name), "SNP Spectrum")) +
    geom_hline(yintercept = 1, linetype = 2, color = "black") +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(color = "black", face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 8),
      panel.spacing.x = unit(0.2, "lines"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  # Raw observed frequency plot
  p2 <- ggplot(prop_df, aes(x = unique_id, fill = strand, y = observed_prop_called)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ group_change, nrow = 1, labeller = labeller(group_change = facet_order), scales = "free_x") +
    scale_x_discrete(
      name = "Trinucleotide Context",
      labels = function(x) {
        sapply(strsplit(x, "_"), function(s) paste0(substr(s[2], 1, 1), "_", substr(s[2], 3, 3)))
      }
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      limits = c(0, max(prop_df$observed_prop_called) * 1.05)
    ) +
    scale_fill_manual(values = c("#E08D7A", "#859DCD")) +
    theme_bw(base_size = 12) +
    labs(y = "Observed Substitution Frequency",
         title = paste(stri_trans_totitle(sample_name), "SNP Spectrum")) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(color = "black", face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 10)),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 8),
      panel.spacing.x = unit(0.2, "lines"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  return(list(fc_plot = p1, observed_plot = p2))
}

# (2)Filter
max_bias <- 0.55
min_bias <- 0.45

variants_high_RawBias <- variant_stats %>% 
  dplyr::filter(raw_mean_bias >= max_bias | raw_mean_bias <= min_bias) %>% 
  dplyr::pull(variant)
length(variants_high_RawBias)

variants_low_RawBias <- variant_stats %>% 
  dplyr::filter(raw_mean_bias < max_bias & raw_mean_bias > min_bias) %>% 
  dplyr::pull(variant)
length(variants_low_RawBias)

spectrum_list1 <- plot_snp_spectrum(mgatk_path = "/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/",
                                    ref_allele_txt = "MT_refAllele.txt",
                                    sample_name = "mouse_High_RawBias", 
                                    called_variants = variants_high_RawBias
)
spectrum_list1[[1]] | spectrum_list1[[2]]

spectrum_list2 <- plot_snp_spectrum(mgatk_path = "/PATH/TO/pbmc.mt/SNP/mouse/mgatk_newCleanBam/final/",
                                    ref_allele_txt = "MT_refAllele.txt",
                                    sample_name = "mouse_Low_RawBias", 
                                    called_variants = variants_low_RawBias
)
spectrum_list2[[1]] | spectrum_list2[[2]]