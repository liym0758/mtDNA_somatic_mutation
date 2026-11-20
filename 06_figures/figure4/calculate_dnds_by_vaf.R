plot_dnds_by_vaf <- function(species, 
                             annovar_file, 
                             vaf_file, 
                             background_file, 
                             output_dir = NULL) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(tidyr)
    library(patchwork)
  })
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # --------------------------
  # 1) annotation: mutation_id -> class
  # --------------------------
  ann <- readr::read_tsv(annovar_file, show_col_types = FALSE) %>%
    dplyr::mutate(
      mutation_id = paste0(Start, Ref, ">", Alt),
      class = dplyr::case_when(
        ExonicFunc.ensGene == "nonsynonymous SNV" ~ "nonsynonymous",
        ExonicFunc.ensGene == "synonymous SNV" ~ "synonymous",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::distinct(mutation_id, class)
  
  # --------------------------
  # 2) vaf: mutation-level mean VAF per celltype_group -> bin
  # --------------------------
  vaf_raw <- readr::read_tsv(vaf_file, show_col_types = FALSE) %>%
    dplyr::mutate(
      mutation_id = mutation,
      celltype_group = dplyr::case_when(
        main_celltype2 %in% c("B cell", "T cell", "T/NK cell", "NK cell") ~ "T_B_NK",
        main_celltype2 == "Monocyte" ~ "Mono",
        TRUE ~ "Other"
      )
    ) %>%
    dplyr::filter(celltype_group %in% c("T_B_NK", "Mono")) %>%
    dplyr::filter(!is.na(frequency), is.finite(frequency))
  
  vaf_mut <- vaf_raw %>%
    dplyr::group_by(celltype_group, mutation_id) %>%
    dplyr::summarise(mean_vaf = mean(frequency, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(
      bin = dplyr::case_when(
        mean_vaf > 0   & mean_vaf <= 0.1 ~ "(0,0.1]",
        mean_vaf > 0.1 & mean_vaf <= 0.3 ~ "(0.1,0.3]",
        mean_vaf > 0.3 & mean_vaf <= 0.5 ~ "(0.3,0.5]",
        mean_vaf > 0.5 & mean_vaf <= 1   ~ "(0.5,1]",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(bin))
  
  combined_dnds_data <- vaf_mut %>%
    dplyr::left_join(ann, by = "mutation_id") %>%
    dplyr::filter(!is.na(class))
  
  readr::write_tsv(
    combined_dnds_data,
    file.path(output_dir, paste0(species, "_dnds_combined_data.tsv"))
  )
  
  # --------------------------
  # Background expectation (once)
  # --------------------------
  bg <- readr::read_csv(background_file, show_col_types = FALSE) %>%
    dplyr::filter(species == !!species)
  if (nrow(bg) != 1) stop("background_file in species=", species, " matched ", nrow(bg), " rows；expected 1 rows. ")
  
  NS_ratio_expected <- bg$non_synonymous_sites / bg$synonymous_sites
  p0 <- NS_ratio_expected / (1 + NS_ratio_expected)
  
  # helper: stars
  p_to_star <- function(p) {
    dplyr::case_when(
      is.na(p) ~ "",
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ ""
    )
  }
  
  # --------------------------
  # 3) per bin N/S + binom p-value (per point)
  # --------------------------
  process_dnds_by_celltype <- function(data, celltype) {
    bin_levels <- c("(0,0.1]", "(0.1,0.3]", "(0.3,0.5]", "(0.5,1]")
    
    df <- data %>%
      dplyr::filter(celltype_group == celltype) %>%
      dplyr::group_by(bin, class) %>%
      dplyr::summarise(count = dplyr::n_distinct(mutation_id), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = class, values_from = count, values_fill = 0) %>%
      tidyr::complete(bin = bin_levels, fill = list(nonsynonymous = 0, synonymous = 0))
    
    if (!"nonsynonymous" %in% names(df)) df$nonsynonymous <- 0
    if (!"synonymous" %in% names(df)) df$synonymous <- 0
    
    df %>%
      dplyr::mutate(
        total_mut = nonsynonymous + synonymous,
        NS_ratio_obs = dplyr::if_else(synonymous == 0, NA_real_, nonsynonymous / synonymous),
        NS_ratio_expected = NS_ratio_expected,
        dnds = NS_ratio_obs / NS_ratio_expected,
        species = species,
        celltype_group = celltype
      ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        p_binom = if (total_mut == 0) NA_real_ else stats::binom.test(nonsynonymous, total_mut, p = p0)$p.value,
        star = p_to_star(p_binom)
      ) %>%
      dplyr::ungroup()
  }
  
  df_T_B_NK <- process_dnds_by_celltype(combined_dnds_data, "T_B_NK")
  df_Mono   <- process_dnds_by_celltype(combined_dnds_data, "Mono")
  
  dnds_results <- dplyr::bind_rows(df_T_B_NK, df_Mono)
  readr::write_tsv(dnds_results, file.path(output_dir, paste0(species, "_dnds_results_withP.tsv")))
  
  # --------------------------
  # 4) Plot: same dot plot + star above points
  # --------------------------
  create_dnds_plot <- function(df, celltype) {
    df$bin <- factor(df$bin, levels = c("(0,0.1]", "(0.1,0.3]", "(0.3,0.5]", "(0.5,1]"))
    
    df_total <- df %>%
      dplyr::summarise(
        total_nonsyn = sum(nonsynonymous, na.rm = TRUE),
        total_syn = sum(synonymous, na.rm = TRUE),
        total_NS_ratio = ifelse(total_syn == 0, NA, total_nonsyn / total_syn),
        total_dnds = total_NS_ratio / unique(NS_ratio_expected)
      )
    
    df <- df %>%
      dplyr::mutate(
        total_mut_plot = pmin(total_mut, 100),
        legend_key = dplyr::case_when(
          celltype_group == "T_B_NK" ~ "T_B_NK",
          celltype_group == "Mono"   ~ "Mono"
        ),
        star_y = pmin(dnds + 0.18, 2.95)
      )
    
    ggplot(df, aes(x = bin, y = dnds, group = 1)) +
      geom_point(
        aes(size = total_mut_plot, fill = legend_key),
        shape = 21, color = "black", stroke = 1.1
      ) +
      geom_text(
        aes(y = star_y, label = star),
        size = 5, color = "black"
      ) +
      scale_fill_manual(
        values = c("T_B_NK" = "#8E90B9", "Mono" = "#7FB0B4"),
        breaks = c("T_B_NK", "Mono"),
        limits = c("T_B_NK", "Mono"),
        drop = FALSE
      ) +
      scale_size_continuous(
        limits = c(0, 100),
        range  = c(2.2, 7.0),
        breaks = c(20, 40, 60, 80, 100),
        labels = c("<=20", "40", "60", "80", ">=100")
      ) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
      geom_text(
        data = df_total,
        aes(x = 2, y = 3.5, label = sprintf("total dN/dS = %.2f", total_dnds)),
        inherit.aes = FALSE, size = 3.5
      ) +
      scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 0.5)) +
      theme_bw(base_size = 12) +
      labs(
        title = paste0(species, " - ", celltype),
        x = "Mean VAF bin (per mutation)",
        y = "dN/dS",
        fill = "Group",
        size = "# mutations in bin"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
  }
  
  p1 <- create_dnds_plot(df_T_B_NK, "T_B_NK") + theme(legend.position = "none")
  p2 <- create_dnds_plot(df_Mono, "Mono")
  
  combined_plot <- (p1 + p2) + patchwork::plot_layout(nrow = 1, guides = "collect") &
    theme(legend.position = "right")
  
  ggsave(
    filename = file.path(output_dir, paste0(species, "_dnds_by_vaf_SplitCellType_withStar.pdf")),
    plot = combined_plot,
    width = 8, height = 3
  )
  
  return(list(df_T_B_NK = df_T_B_NK, df_Mono = df_Mono, dnds_results = dnds_results))
}


# =========================
# Batch run: Dot plot ----------------------
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

species_list <- c('human','mouse','rat','chicken','cobra','sturgeon','stingray')
base_dir <- "/PATH/TO/pbmc.mt/SNP"

# annovar annotation file list
annovar_annotation_file_list <- c(
  "human" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/human/result/hg38_anno.hg38_chrM_multianno.txt",
  "mouse" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/mouse/result/mm39_anno.mm39_MT_multianno.txt",
  "rat" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/rat/result/rn7_MT_anno.rn7_MT_multianno.txt",
  "chicken" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/chicken/result/gg1_MT_anno.gg1_MT_multianno.txt",
  "cobra" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/cobra/result/na_MT_anno.na_MT_multianno.txt",
  "sturgeon" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/sturgeon/result/as_CNCB_MT_anno.as_CNCB_MT_multianno.txt",
  "stingray" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/stingray/result/pl_MT_anno.pl_MT_multianno.txt"
)

# background file (PAML sites summary)
background_file <- file.path(base_dir, "mgatk.downstream/Variant_Type_Anno/PAML/background_sites_summary.csv")

# VAF per-cell table folder (your fixedDepth_VAF output)
vaf_dir <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human"

# output dirs
out_dot <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/plots/dNdS_withP_20260220/DotPlot/"
out_bar <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/plots/dNdS_withP_20260220/BarPlot/"

dir.create(out_dot, recursive = TRUE, showWarnings = FALSE)
dir.create(out_bar, recursive = TRUE, showWarnings = FALSE)

# containers for optional combined tables
dot_res_all <- list()
bar_res_all <- list()

for (sp in species_list) {
  message("========== Processing species: ", sp, " ==========")
  
  annovar_file <- unname(annovar_annotation_file_list[sp])
  vaf_file <- file.path(vaf_dir, paste0(sp, "_detail_filtered_df_260128.tsv"))
  
  if (!file.exists(annovar_file)) {
    warning("Skip ", sp, ": annovar_file not found: ", annovar_file)
    next
  }
  if (!file.exists(vaf_file)) {
    warning("Skip ", sp, ": vaf_file not found: ", vaf_file)
    next
  }
  if (!file.exists(background_file)) {
    stop("background_file not found: ", background_file)
  }
  
  # ---- 1) Dot plot (binom stars per point)
  dot_list <- plot_dnds_by_vaf(
    species = sp,
    annovar_file = annovar_file,
    vaf_file = vaf_file,
    background_file = background_file,
    output_dir = out_dot
  )
  
  # collect dot results (TSV already written by function)
  if (!is.null(dot_list$dnds_results)) {
    dot_res_all[[sp]] <- dot_list$dnds_results %>%
      dplyr::mutate(plot_type = "dot")
  }
} 
 

# =========================
# Optional: write combined TSV across all species
# =========================
if (length(dot_res_all) > 0) {
  dot_all_df <- dplyr::bind_rows(dot_res_all)
  readr::write_tsv(dot_all_df, file.path(out_dot, "AllSpecies_dot_dnds_results_withP.tsv"))
}

if (length(bar_res_all) > 0) {
  bar_all_df <- dplyr::bind_rows(bar_res_all)
  readr::write_tsv(bar_all_df, file.path(out_bar, "AllSpecies_bar_dnds_results.tsv"))
}

# =========================
# Optional: merge PDFs (run in shell)
# =========================
# Dot plots:
# pdfunite /.../DotPlot_withStars/*_dnds_by_vaf_SplitCellType_withStar.pdf \
#   /.../DotPlot_withStars/AllSpecies_dnds_by_vaf_DotPlot_withStars.pdf

  
# =========================
# Site-level dN/dS (NO bins) + per-bar significance vs background ----------------
# Clustered bars: x=species, dodge by celltype_group
# ONE panel
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# -------------------------
# Config
# -------------------------
species_list <- c('human','mouse','rat','chicken','cobra','sturgeon','stingray')

annovar_annotation_file_list <- c(
  "human"    = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/human/result/hg38_anno.hg38_chrM_multianno.txt",
  "mouse"    = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/mouse/result/mm39_anno.mm39_MT_multianno.txt",
  "rat"      = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/rat/result/rn7_MT_anno.rn7_MT_multianno.txt",
  "chicken"  = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/chicken/result/gg1_MT_anno.gg1_MT_multianno.txt",
  "cobra"    = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/cobra/result/na_MT_anno.na_MT_multianno.txt",
  "sturgeon" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/sturgeon/result/as_CNCB_MT_anno.as_CNCB_MT_multianno.txt",
  "stingray" = "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/stingray/result/pl_MT_anno.pl_MT_multianno.txt"
)

vaf_dir <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/"
background_file <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/PAML/background_sites_summary.csv"

output_dir <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/plots/dNdS_withP_20260220/BarPlot/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------
# Helper: p -> stars
# -------------------------
p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  )
}

# -------------------------
# Helper: read ANNOVAR and get class per mutation_id
# -------------------------
read_annovar_class <- function(annovar_file) {
  readr::read_tsv(annovar_file, show_col_types = FALSE) %>%
    dplyr::mutate(
      mutation_id = paste0(Start, Ref, ">", Alt),
      class = dplyr::case_when(
        ExonicFunc.ensGene %in% c("nonsynonymous SNV") ~ "nonsynonymous",
        ExonicFunc.ensGene == "synonymous SNV" ~ "synonymous",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::distinct(mutation_id, class)
}

# -------------------------
# Helper: read per-cell mutation table and create site-level mean heteroplasmy
# NOTE: mean_heteroplasmy is mean over mutated cells (your current definition)
# -------------------------
read_site_meanHet <- function(vaf_file) {
  vaf_cell <- readr::read_tsv(vaf_file, show_col_types = FALSE) %>%
    dplyr::mutate(
      mutation_id = mutation,
      celltype_group = dplyr::case_when(
        main_celltype2 %in% c("B cell", "T cell", "T/NK cell", "NK cell") ~ "T_B_NK",
        main_celltype2 == "Monocyte" ~ "Mono",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(celltype_group)) %>%
    dplyr::filter(!is.na(frequency), is.finite(frequency))
  
  site_vaf <- vaf_cell %>%
    dplyr::group_by(celltype_group, mutation_id) %>%
    dplyr::summarise(
      n_cells_with_mut = dplyr::n(),
      mean_heteroplasmy = sum(frequency) / n_cells_with_mut,
      .groups = "drop"
    ) %>%
    dplyr::filter(mean_heteroplasmy > 0, mean_heteroplasmy <= 1)
  
  site_vaf
}

# -------------------------
# Compute per species dN/dS per celltype_group (NO bin)
# and per-bar p-value vs background using binom.test
# -------------------------
compute_species_dnds_noBin <- function(species, annovar_file, vaf_file, bg_df) {
  ann <- read_annovar_class(annovar_file)
  site_vaf <- read_site_meanHet(vaf_file)
  
  # Count unique sites per class (no bin)
  site_counts <- site_vaf %>%
    dplyr::left_join(ann, by = "mutation_id") %>%
    dplyr::filter(!is.na(class)) %>%
    dplyr::distinct(celltype_group, mutation_id, class) %>%
    dplyr::group_by(celltype_group, class) %>%
    dplyr::summarise(n_sites = dplyr::n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = class, values_from = n_sites, values_fill = 0)
  
  # Background expected ratio for this species
  bg_sub <- bg_df %>% dplyr::filter(species == !!species)
  if (nrow(bg_sub) != 1) {
    stop("Background file does not contain exactly one row for species: ", species,
         " (matched rows=", nrow(bg_sub), ")")
  }
  
  NS_ratio_expected <- bg_sub$non_synonymous_sites / bg_sub$synonymous_sites
  p0 <- NS_ratio_expected / (1 + NS_ratio_expected)
  
  dnds_df <- site_counts %>%
    dplyr::mutate(
      total_mut = nonsynonymous + synonymous,
      NS_ratio_obs = dplyr::if_else(synonymous == 0, NA_real_, nonsynonymous / synonymous),
      dnds = NS_ratio_obs / NS_ratio_expected,
      species = species
    )
  
  # Ensure both celltypes exist for each species (so bars align)
  dnds_df <- dnds_df %>%
    tidyr::complete(
      species = species,
      celltype_group = c("T_B_NK", "Mono"),
      fill = list(nonsynonymous = 0L, synonymous = 0L, total_mut = 0L)
    ) %>%
    dplyr::mutate(
      NS_ratio_obs = dplyr::if_else(synonymous == 0, NA_real_, nonsynonymous / synonymous),
      dnds = NS_ratio_obs / NS_ratio_expected
    )
  
  # Per-bar exact binomial test vs background p0
  dnds_df <- dnds_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      p_vs_bg = if (is.na(total_mut) || total_mut == 0) NA_real_
      else stats::binom.test(nonsynonymous, total_mut, p = p0)$p.value,
      star_vs_bg = p_to_star(p_vs_bg)
    ) %>%
    dplyr::ungroup()
  
  dnds_df
}

# -------------------------
# Read background once
# -------------------------
bg_df <- readr::read_csv(background_file, show_col_types = FALSE)

# -------------------------
# Loop all species and bind
# -------------------------
all_dnds <- list()

for (sp in species_list) {
  message("Processing species: ", sp)
  
  annovar_file <- unname(annovar_annotation_file_list[sp])
  vaf_file <- file.path(vaf_dir, paste0(sp, "_detail_filtered_df_260128.tsv"))
  
  if (!file.exists(annovar_file)) {
    warning("Missing annovar file: ", annovar_file, " -> skip ", sp)
    next
  }
  if (!file.exists(vaf_file)) {
    warning("Missing vaf file: ", vaf_file, " -> skip ", sp)
    next
  }
  
  one <- compute_species_dnds_noBin(
    species = sp,
    annovar_file = annovar_file,
    vaf_file = vaf_file,
    bg_df = bg_df
  )
  
  all_dnds[[sp]] <- one
}

plot_df <- dplyr::bind_rows(all_dnds) %>%
  dplyr::mutate(
    species = factor(str_to_title(species), levels = str_to_title(species_list)),
    celltype_group = factor(celltype_group, levels = c("T_B_NK", "Mono")),
    legend_key = dplyr::case_when(
      total_mut >= 20 & celltype_group == "T_B_NK" ~ "T_B_NK (>=20 muts)",
      total_mut >= 20 & celltype_group == "Mono"   ~ "Mono (>=20 muts)",
      TRUE                                         ~ "<20 muts"
    ),
    legend_key = factor(legend_key, levels = c("T_B_NK (>=20 muts)", "Mono (>=20 muts)", "<20 muts"))
  )

# Save underlying table (with p-values)
out_tsv <- file.path(output_dir, "all_species_sitelevel_dnds_noBin_withP.tsv")
readr::write_tsv(plot_df, out_tsv)

# -------------------------
# Plotting: single panel
# -------------------------
make_panel_noBin <- function(df, show_legend = TRUE) {
  dodge <- position_dodge(width = 0.75)
  
  ggplot(df, aes(x = species, y = dnds)) +
    geom_col(
      aes(group = celltype_group, fill = legend_key),
      position = dodge,
      width = 0.65,
      color = "black",
      show.legend = show_legend
    ) +
    # per-bar significance stars vs background
    geom_text(
      aes(
        group = celltype_group,
        label = star_vs_bg,
        y = dnds + 0.1
      ),
      position = dodge,
      size = 5,
      color = "black"
    ) +
    scale_fill_manual(
      values = c(
        "T_B_NK (>=20 muts)" = "#8E90B9",
        "Mono (>=20 muts)"   = "#7FB0B4",
        "<20 muts"           = "white"
      ),
      guide = if (show_legend) guide_legend(override.aes = list(color = "black")) else "none"
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    scale_y_continuous(
      limits = c(0, 1.8),
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    ) +
    labs(
      title = "Overall mtDNA dN/dS",
      x = "Species",
      y = "dN/dS",
      fill = NULL
    )
}

p <- make_panel_noBin(plot_df, show_legend = TRUE)

out_pdf <- file.path(output_dir, "AllSpecies_sitelevel_dNdS_noBin_cluster_barplot_withP.pdf")
ggsave(out_pdf, p, width = 8, height = 3)

message("Done.")
message("Output PDF:  ", out_pdf)
message("Output table:", out_tsv)