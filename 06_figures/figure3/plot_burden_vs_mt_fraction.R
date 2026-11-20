library(ggplot2)
library(dplyr)
library(cowplot)
mt.celltype.colors <- c('#B3659F','#8E90B9', '#CCE0F5','#3361A5','#7FB0B4')
names(mt.celltype.colors) <- c("B cell","T cell","NK cell","T/NK cell","Monocyte")

celltypes <- c("Monocyte","B cell","T cell","T/NK cell","NK cell")
species_list <- c("human","mouse","rat","chicken","cobra","sturgeon","stingray")

plot_data <- read_tsv("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/mutation_burden/scmtMPM_newHuman/PerCell_scmtMPM_260128.tsv")

plot_data2 <- plot_data %>%
  mutate(
    main_celltype2 = factor(main_celltype2, levels = celltypes),
    species = factor(species, levels = species_list)
  )



make_species_row <- function(sp) {
  df <- plot_data2 %>% filter(species == sp)
  
  ggplot(df, aes(x = mt.percent, y = n_mut)) +
    geom_point(
      aes(colour = main_celltype2),
      shape = 16,
      size = 2.5,
      alpha = 1
    ) +
    scale_colour_manual(values = mt.celltype.colors)+
    geom_smooth(
      method = "lm",
      se = FALSE,
      colour = "black",
      linewidth =1
    )+
    facet_wrap(~ main_celltype2, nrow = 1, drop = TRUE) +
    scale_fill_manual(values = mt.celltype.colors, drop = TRUE) +
    coord_cartesian(xlim = c(0, 80), ylim = c(0, 14)) +
    scale_x_continuous(limits = c(0, 80), expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(limits = c(0, 14),  expand = expansion(mult = c(0, 0.02))) +
    labs(
      title = str_to_title(sp),
      x = "MT%",
      y = "Number of mutation"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = NA, colour = "black"),
      strip.text = element_text(size=12, hjust = 0.5, face = "bold"),
      text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5, vjust = 0.5, size=12),
      axis.text.y = element_text(color = "black")
    )
  
}

library(patchwork)
library(dplyr)

max_panels <- 4  # 5

plots_with_n <- lapply(species_list, function(sp) {
  df_sp <- plot_data2 %>% filter(species == sp)
  n_panels <- df_sp %>%
    filter(!is.na(main_celltype2)) %>%
    distinct(main_celltype2) %>%
    nrow()
  
  p <- make_species_row(sp)
  
  list(p = p, n = n_panels)
})
names(plots_with_n) <- species_list

rows <- lapply(plots_with_n, function(x) {
  n <- x$n
  if (is.na(n) || n < 1) n <- 1
  
  (x$p + plot_spacer()) +
      plot_layout(widths = c(n, max_panels - n))
})

final_plot <- wrap_plots(rows, ncol = 1)

final_plot

pdf("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/mutation_burden/scmtMPM_newHuman/plots/MtPct_vs_nMut_dot.pdf", width = 8, height = 17)
print(final_plot)
dev.off()



library(ggplot2)
library(dplyr)
library(cowplot)
mt.celltype.colors <- c('#B3659F','#8E90B9', '#CCE0F5','#3361A5','#7FB0B4')
names(mt.celltype.colors) <- c("B cell","T cell","NK cell","T/NK cell","Monocyte")

celltypes <- c("Monocyte","B cell","T cell","T/NK cell","NK cell")
species_list <- c("human","mouse","rat","chicken","cobra","sturgeon","stingray")

plot_data <- read_tsv("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/mutation_burden/scmtMPM_newHuman/PerCell_scmtMPM_260128.tsv")

plot_data2 <- plot_data %>%
  mutate(
    main_celltype2 = factor(main_celltype2, levels = celltypes),
    species = factor(species, levels = species_list),
    log10_mtreads = log10(mt_reads)
  )



make_species_row <- function(sp) {
  df <- plot_data2 %>% dplyr::filter(species == sp)
  
  stat_df <- df %>%
    dplyr::filter(!is.na(mt.percent), !is.na(scmtMPM)) %>%
    dplyr::group_by(main_celltype2) %>%
    dplyr::summarise(
      n = dplyr::n(),
      r2 = {
        fit <- stats::lm(scmtMPM ~ mt.percent)
        summary(fit)$r.squared
      },
      p = {
        fit <- stats::lm(scmtMPM ~ mt.percent)
        stats::coef(summary(fit))[2, 4]
      },
      x = -Inf,
      y = Inf,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = paste0("R²=", sprintf("%.3f", r2),
                     "\nP=", format.pval(p, digits = 2, eps = 1e-3))
    )
  
  ggplot(df, aes(x = mt.percent, y = scmtMPM)) +
    geom_point(
      aes(colour = main_celltype2),
      shape = 16,
      size = 2.5,
      alpha = 1
    ) +
    scale_colour_manual(values = mt.celltype.colors) +
    geom_smooth(
      method = "lm",
      se = FALSE,
      colour = "black",
      linewidth = 1
    ) +
    geom_text(
      data = stat_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.10,
      size = 3.2,
      colour = "black"
    ) +
    facet_wrap(~ main_celltype2, nrow = 1, drop = TRUE) +
    labs(
      title = stringr::str_to_title(sp),
      x = "mt%",
      y = "scmtMPM"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = NA, colour = "black"),
      strip.text = element_text(size = 12, hjust = 0.5, face = "bold"),
      text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(color = "black"),
      axis.title.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black", angle = 0, hjust = 0.5, vjust = 0.5, size = 12),
      axis.text.y = element_text(color = "black")
    )
}

library(patchwork)
library(dplyr)

max_panels <- 4  # 5

plots_with_n <- lapply(species_list, function(sp) {
  df_sp <- plot_data2 %>% filter(species == sp)
  n_panels <- df_sp %>%
    filter(!is.na(main_celltype2)) %>%
    distinct(main_celltype2) %>%
    nrow()
  
  p <- make_species_row(sp)
  
  list(p = p, n = n_panels)
})
names(plots_with_n) <- species_list

rows <- lapply(plots_with_n, function(x) {
  n <- x$n
  if (is.na(n) || n < 1) n <- 1
  
  (x$p + plot_spacer()) +
      plot_layout(widths = c(n, max_panels - n))
})

final_plot <- wrap_plots(rows, ncol = 1)

#final_plot

pdf("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/mutation_burden/scmtMPM_newHuman/plots/mtpct_vs_scmtMPM_dot.pdf", width = 8, height = 17)
print(final_plot)
dev.off()