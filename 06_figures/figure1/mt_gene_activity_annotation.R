setwd("/PATH/TO/pbmc.mt/celltype.anno/")
library(ggrepel)
library(ggpubr)

my24colors <- c('#7EA475','#C4DAB0','#3B6E4F', '#7EA795', '#699B92', '#274B2F', 
                '#76B2BE','#C3D1DD','#6F8BB0', '#4F418E','#8291B0','#ACBDDE',
                '#DCE572','#CFC436','#E5DBB9', '#D1B56C', '#C6902E', '#AF651E', 
                '#A63376','#853280', '#B3659F','#D490B9', '#AF95C4', '#959492')



#########readRDS#############
#new_annotation
human <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/final_ATAC_Human.rds")
mouse <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/new_ATAC/final_ATAC_Mouse.rds")
rat <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/new_ATAC/final_ATAC_Rat.rds")
chicken <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/new_ATAC/final_ATAC_Chicken.rds")
cobra <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/new_ATAC/final_ATAC_Snake.rds")
sturgeon <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/new_ATAC/final_ATAC_XunFish.rds")
stingray <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/new_ATAC/final_ATAC_HongFish.rds")



objlist <- list(mouse=mouse,rat=rat,chicken=chicken,cobra=cobra,sturgeon=sturgeon,stingray=stingray,human=human)
species <- c('mouse','rat','chicken','cobra','sturgeon','stingray','human')

rm(mouse,rat,chicken,cobra,sturgeon,stingray,human)
gc()

objlist[[7]]@active.ident <- objlist[[7]]$Joint_CellType
all_celltypes <- c()
for (i in 1:length(objlist)){
  print(names(objlist)[i])
  #print(unique(objlist[[i]]@active.ident))
  all_celltypes <- c(all_celltypes, levels(objlist[[i]]@active.ident))
}
all_celltypes <- unique(all_celltypes)

# 1 coarse cell type ---------------
maintype_rename_mapping <- list(
  "B" = "B cell",
  "B cell" = "B cell",
  "B Memory" = "B cell",
  "B Naive" = "B cell",
  "Plsam" = "B cell",
  
  "CD4+ T" = "T cell",
  "CD4+ Naive T" = "T cell",
  "CD4+ Memory T" = "T cell",
  "CD4+ Memory" = "T cell",
  "CD4+ Naive" = "T cell",
  "CD8+ T" = "T cell",
  "CD8+ Memory" = "T cell",
  "CD8+ Naive" = "T cell",
  "delta+ T" = "T cell",
  "MAIT" = "T cell",
  "NKT" = "T cell",
  
  "NK" = "NK cell",
  "NK cell" = "NK cell",
  

  "T_NK" = "T/NK cell",

  
  "Monocyte" = "Myeloid cell",
  "CD14+ Monocyte" = "Myeloid cell",
  "CD16+ Monocyte" = "Myeloid cell",
  "CD14+ Mono" = "Myeloid cell",
  "CD16+ Mono" = "Myeloid cell",
  "Monocyte_1" = "Myeloid cell",
  "Monocyte_2" = "Myeloid cell",  

  
  "cDC" = "Myeloid cell",
  "pDC" = "Myeloid cell",
  "Neutrophil" = "Myeloid cell",
  
  
  "Platelet" = "Platelet",
  "Erythrocyte" = "Erythrocyte",
  
  "HSP" = "Others",
  "Others" = "Others"
)

main_type_order <- c("B cell","T cell","NK cell","T/NK cell","Myeloid cell","Platelet","Erythrocyte","Others")


# 3 sub cell type ------------------------------
subtype_rename_mapping <- list(
  "B" = "B cell",
  "B cell" = "B cell",
  "B Memory" = "Memory B",
  "B Naive" = "Naive B",
  "Plsam" = "Plasma B",
  
  "CD4+ T" = "CD4+ T",
  "CD4+ Naive T" = "CD4+ Naive T",
  "CD4+ Memory T" = "CD4+ Memory T",
  "CD4+ Memory" = "CD4+ Memory T",
  "CD4+ Naive" = "CD4+ Naive T",
  "CD8+ T" = "CD8+ T",
  "CD8+ Memory" = "CD8+ Memory T",
  "CD8+ Naive" = "CD8+ Naive T",
  "delta+ T" = "Delta T",
  "MAIT" = "MAIT",
  "NKT" = "NKT",
  
  "NK" = "NK cell",
  "NK cell" = "NK cell",
  

  "T_NK" = "T/NK cell",

  
  "Monocyte" = "Monocyte",
  "CD14+ Monocyte" = "CD14+ Monocyte",
  "CD16+ Monocyte" = "CD16+ Monocyte",
  "CD14+ Mono" = "CD14+ Monocyte",
  "CD16+ Mono" = "CD16+ Monocyte",
  "Monocyte_1" = "Monocyte_1",
  "Monocyte_2" = "Monocyte_2",  
  
  
  "cDC" = "cDC",
  "pDC" = "pDC",
  "Neutrophil" = "Neutrophil",
  
  
  "Platelet" = "Platelet",
  "Erythrocyte" = "Erythrocyte",
  
  "HSP" = "HSP",
  "Others" = "Others"
)

sub_type_order <- c(
    "HSP",
    
    "B cell",
    "Naive B",
    "Memory B",
    "Plasma B",
    
    "CD4+ T",
    "CD4+ Naive T",
    "CD4+ Memory T",
    "CD8+ T",
    "CD8+ Naive T",
    "CD8+ Memory T",
    "Delta T",
    "MAIT",
    "NKT",
    
    "T/NK cell",
    
    "NK cell",
    
    "Monocyte",
    "Monocyte_1",
    "Monocyte_2",
    "CD14+ Monocyte",
    "CD16+ Monocyte",

    
    "cDC",
    "pDC",
    
    "Neutrophil",
    
    "Platelet",
    "Erythrocyte",
    
    "Others"
  )

# Rename 
for (i in seq_along(objlist)) {
  # objlist[[i]]$barcode <- rownames(objlist[[i]]@meta.data)
  # objlist[[i]]$orig_celltype <- objlist[[i]]@active.ident
  
  current_ident <- objlist[[i]]$orig_celltype
  
  new_ident1 <- as.character(current_ident)
  mask1 <- new_ident1 %in% names(maintype_rename_mapping)
  new_ident1[mask1] <- maintype_rename_mapping[new_ident1[mask1]]
  
  new_ident2 <- as.character(current_ident)
  mask2 <- new_ident2 %in% names(subtype_rename_mapping)
  new_ident2[mask2] <- subtype_rename_mapping[new_ident2[mask2]]
  
  
  objlist[[i]]$main_celltype <- factor(
    new_ident1,
    levels = main_type_order
  )
  
  objlist[[i]]$sub_celltype <- factor(
    new_ident2,
    levels = sub_type_order
  )
}

for (i in 1:length(objlist)) {
  print(names(objlist)[i])
  print(unique(objlist[[i]]$main_celltype))
  print(unique(objlist[[i]]$sub_celltype))
}





meta.list_old <- readRDS("/PATH/TO/pbmc.mt/celltype.anno/seurat.meta/meta.list.mt.count.mtMask.RDS")

View(meta.list_old[[1]])

mt_cols <- list()
meta.list <- list()

for (i in 1:6){
  print(unique(meta.list_old[[i]]$species))
  mt_cols[[i]] <- meta.list_old[[i]][,c("barcode","total_reads","mt_reads","mt.percent","species")]
  
  meta.list[[i]] <- objlist[[i]]@meta.data
  meta.list[[i]] <- merge(meta.list[[i]],mt_cols[[i]], by = c("barcode"))
  names(meta.list)[i] <- unique(meta.list_old[[i]]$species)
  write.csv(meta.list[[i]], file = paste("./new_anno_ATAC/seurat_meta/",names(meta.list)[i],".meta.csv",sep = ''))
}

saveRDS(meta.list, "./new_anno_ATAC/seurat_meta/meta.list.RDS")


for (i in 1:7){print(nrow(objlist[[i]]@meta.data)-nrow(meta.list[[i]]))}
for (i in 1:7){
  print(names(objlist)[i])
  objlist[[i]]@meta.data <- meta.list[[i]]
  rownames(objlist[[i]]@meta.data) <- objlist[[i]]@meta.data$barcode
  
  saveRDS(objlist[[i]],file = paste0("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/RDS/",names(objlist)[i],".ATAC.RDS"))
}
saveRDS(objlist, file = "./new_anno_ATAC/RDS/objlist.RDS")


objlist <- readRDS("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/RDS/objlist.RDS")



maintype_rename_mapping2 <- list(
  "B cell" = "B cell",
  "Memory B" = "B cell",
  "Naive B" = "B cell",
  "Plasma B" = "B cell",
  
  "CD4+ T" = "T cell",
  "CD4+ Naive T" = "T cell",
  "CD4+ Memory T" = "T cell",
  "CD4+ Naive T" = "T cell",
  "CD8+ T" = "T cell",
  "CD8+ Memory T" = "T cell",
  "CD8+ Naive T" = "T cell",
  "Delta T" = "T cell",
  "MAIT" = "T cell",
  "NKT" = "T cell",
  
  "NK cell" = "NK cell",
  
  "T/NK cell" = "T/NK cell",
  
  
  "Monocyte" = "Monocyte",
  "CD14+ Monocyte" = "Monocyte",
  "CD16+ Monocyte" = "Monocyte",
  "Monocyte_1" = "Monocyte",
  "Monocyte_2" = "Monocyte",  
  
  
  "cDC" = "DC",
  "pDC" = "DC",
  "Neutrophil" = "Neutrophil",
  
  
  "Platelet" = "Platelet",
  "Erythrocyte" = "Erythrocyte",
  
  "HSP" = "Progenitor",
  "Others" = "Others"
)

main_type_order2 <- c("B cell","T cell","NK cell","T/NK cell","Monocyte","DC","Neutrophil","Platelet","Erythrocyte","Progenitor","Others")

#objlist <- list()
for (sp in species){
  #objlist[[sp]] <- readRDS(file = paste0("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/RDS/",sp,".ATAC.RDS"))

  
  current_ident <- objlist[[sp]]$sub_celltype
  
  new_ident1 <- as.character(current_ident)
  mask1 <- new_ident1 %in% names(maintype_rename_mapping2)
  new_ident1[mask1] <- maintype_rename_mapping2[new_ident1[mask1]]
  
  
  objlist[[sp]]$main_celltype2 <- factor(
    new_ident1,
    levels = main_type_order2
  )
  
  print(sp)
  print(unique(objlist[[sp]]$main_celltype2))

}

for (sp in species){
saveRDS(objlist[[sp]],file = paste0("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/RDS/",sp,".ATAC.RDS"))
}
saveRDS(objlist, file = "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/RDS/objlist.RDS")

meta.list <- list()
for (sp in species){
  print(sp)
  meta.list[[sp]] <- objlist[[sp]]@meta.data
  write.csv(meta.list[[sp]], file = paste("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/",sp,".meta.csv",sep = ''))
}

saveRDS(meta.list, "/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/meta.list.RDS")

sub_celltype_colors <- c(
  "HSP" = "#E63863",
  
  "B cell" = "#D490B9",
  "Naive B" = "#E5D2DD",
  "Memory B" = "#B3659F",
  "Plasma B" = "#853280",
  
  "CD4+ T" = "#ACBDDE",
  "CD4+ Naive T" = "#CCE0F5",
  "CD4+ Memory T" = "#8291B0",
  "CD8+ T" = "#AF95C4",
  "CD8+ Naive T" = "#CCC9E6",
  "CD8+ Memory T" = "#8C549C",
  
  "NK cell" = "#A63376",
  
  "T/NK cell" = "#8E90B9",
  "Delta T" = "#4F418E",
  "MAIT" = "#58A4C3",
  "NKT" = "#5F3D69",
  
  "Monocyte" = "#7EA795",
  "Monocyte_1" = "#7EA475",
  "Monocyte_2" = "#3B6E4F",
  "CD14+ Monocyte" = "#C4DAB0",
  "CD16+ Monocyte" = "#7EA475",
  
  "Neutrophil" = "#274B2F",
  
  "cDC" = "#DEA091",
  "pDC" = "#C6902E",
  
  "Platelet" = "#D1B56C",
  "Erythrocyte" = "#712820",
  "Others" = "#959492"
)

p <- list()
for (i in 1:7){
  umap_df <- Embeddings(objlist[[i]], reduction = "umap") %>% 
    as.data.frame() %>%
    setNames(c("UMAP_1", "UMAP_2")) %>%
    cbind(objlist[[i]]@meta.data[c("main_celltype", "sub_celltype")])
  
  label_coords <- umap_df %>%
    group_by(main_celltype) %>%
    summarise(
      UMAP_1 = median(UMAP_1),
      UMAP_2 = median(UMAP_2)
    )
  
  p[[i]] <- DimPlot(objlist[[i]], 
                     label = FALSE,
                     repel = TRUE, 
                     reduction = "umap",
                     group.by = "sub_celltype",
                     cols = sub_celltype_colors, 
                     pt.size = 0.3)+
    theme(text = element_text(size = 20), strip.text = element_text (size = 20))+
    ggtitle(paste("ATAC UMAP of",str_to_title(species[i]),"PBMC"))+ 
    geom_text_repel(
      data = label_coords,
      aes(x = UMAP_1, y = UMAP_2, label = main_celltype),
      size = 5,
      color = "black",
    )

}

pdf("./new_anno_ATAC/umaps/celltype.umap.pdf",width=8, height=6)
for (i in 1:7){plot(p[[i]])}
dev.off()

# mt.celltype.colors <- c('#CCC9E6', '#CCE0F5','#9CBAD8','#7FB0B4', '#C4DAB0')
# names(mt.celltype.colors) <- c("B cell","T cell","NK cell","T/NK cell","Myeloid cell")
mt.celltype.colors <- c('#B3659F','#8E90B9', '#CCE0F5','#3361A5','#7FB0B4')
names(mt.celltype.colors) <- c("B cell","T cell","NK cell","T/NK cell","Monocyte")

for (i in 1:7){
  p[[i]] <- VlnPlot(
    object = objlist[[i]],
    features = c("mt.percent"),
    group.by = "main_celltype2",
    pt.size = 0)+ggtitle(paste(str_to_title(paste("Mitochondiral% of cell types in",species[i]),"PBMC")))
}


combined_data <- data.frame()
for (i in 1:7) {
  species_data <- p[[i]][[1]][["data"]]
  species_data$species <- rep(species[i],nrow(species_data))
  combined_data <- rbind(combined_data, species_data)

}
colnames(combined_data) <- c("mt.percent","Cell_Type","Species")
combined_data$Species <- factor(combined_data$Species, levels = c(species[7],species[1:6]))
combined_data <- combined_data %>% dplyr::filter(Cell_Type %in% c("B cell","T cell","NK cell","T/NK cell","Monocyte"))
combined_data$Cell_Type <- factor(combined_data$Cell_Type,levels = c("B cell","T cell","NK cell","T/NK cell","Monocyte"))

all_species <- c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray")
all_cell_types <- c("B cell","T cell","NK cell","T/NK cell","Monocyte")
combined_data$Species <- str_to_title(combined_data$Species)


df <- combined_data %>%
  mutate(
    Species = factor(Species, levels = all_species),
    Cell_Type = factor(Cell_Type, levels = all_cell_types)
  ) %>%
  complete(Species, Cell_Type)


#df <- read.csv(file = '/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/mt.percent/mt.pct.vln.MainCelltype2.data.csv', row.names = 1)

df$Species <- factor(df$Species, levels = c('Human','Mouse','Rat','Chicken','Cobra','Sturgeon','Stingray'))

df_clean <- df %>% dplyr::filter(!is.na(mt.percent))

celltype_medians <- df_clean %>%
  group_by(Cell_Type) %>%
  dplyr::summarise(median_mt = mean(mt.percent, na.rm = TRUE)) %>%
  arrange(median_mt)

df_clean$Cell_Type <- factor(df_clean$Cell_Type, levels = celltype_medians$Cell_Type)

create_signif_layer <- function(species_name, data) {
  species_data <- data %>% filter(Species == species_name)
  
  available_celltypes <- as.character(unique(species_data$Cell_Type))
  
  if (species_name %in% c("Human", "Mouse", "Rat", "Chicken")) {
    valid_comparisons <- list()
    if (length(available_celltypes) >= 2) {
      all_pairs <- combn(available_celltypes, 2, simplify = FALSE)
      
      for (pair in all_pairs) {
        cell1 <- pair[1]
        cell2 <- pair[2]
        
        cell1_data <- species_data %>% filter(as.character(Cell_Type) == cell1)
        cell2_data <- species_data %>% filter(as.character(Cell_Type) == cell2)
        
        if (nrow(cell1_data) >= 2 && nrow(cell2_data) >= 2) {
          valid_comparisons <- c(valid_comparisons, list(c(cell1, cell2)))
        }
      }
    }
  } else {
    valid_comparisons <- list()
    
    if ("Monocyte" %in% available_celltypes && "B cell" %in% available_celltypes) {
      monocyte_data <- species_data %>% filter(as.character(Cell_Type) == "Monocyte")
      bcell_data <- species_data %>% filter(as.character(Cell_Type) == "B cell")
      
      if (nrow(monocyte_data) >= 2 && nrow(bcell_data) >= 2) {
        valid_comparisons <- c(valid_comparisons, list(c("Monocyte", "B cell")))
      }
    }
    
    if ("T/NK cell" %in% available_celltypes && "B cell" %in% available_celltypes) {
      tnk_data <- species_data %>% filter(as.character(Cell_Type) == "T/NK cell")
      bcell_data <- species_data %>% filter(as.character(Cell_Type) == "B cell")
      
      if (nrow(tnk_data) >= 2 && nrow(bcell_data) >= 2) {
        valid_comparisons <- c(valid_comparisons, list(c("T/NK cell", "B cell")))
      }
    }
    
    if ("T/NK cell" %in% available_celltypes && "Monocyte" %in% available_celltypes) {
      tnk_data <- species_data %>% filter(as.character(Cell_Type) == "T/NK cell")
      monocyte_data <- species_data %>% filter(as.character(Cell_Type) == "Monocyte")
      
      if (nrow(tnk_data) >= 2 && nrow(monocyte_data) >= 2) {
        valid_comparisons <- c(valid_comparisons, list(c("T/NK cell", "Monocyte")))
      }
    }
  }
  
  if (length(valid_comparisons) > 0) {
    return(
      geom_signif(
        data = species_data,
        comparisons = valid_comparisons,
        test = "wilcox.test",
        test.args = list(alternative = "two.sided"),
        map_signif_level = function(p) {
          ifelse(p < 0.001, "***", 
                 ifelse(p < 0.01, "**", 
                        ifelse(p < 0.05, "*", "ns")))
        },
        textsize = 3,
        vjust = 0.2,
        step_increase = 0.1,
        na.rm = TRUE
      )
    )
  } else {
    return(NULL)
  }
}

species_list <- unique(df_clean$Species)

signif_layers <- lapply(species_list, function(sp) {
  layer <- create_signif_layer(sp, df_clean)
  if (!is.null(layer)) {
    cat("Species:", sp, "- Comparisons:", 
        paste(sapply(layer$comparisons, function(x) paste(x, collapse=" vs ")), 
              collapse=", "), "\n")
  } else {
    cat("Species:", sp, "- No valid comparisons\n")
  }
  layer
})

signif_layers <- signif_layers[!sapply(signif_layers, is.null)]

p <- ggplot(df_clean, aes(x = Cell_Type, y = mt.percent, fill = Cell_Type)) +
  geom_violin(
    trim = TRUE,
    scale = "width",
    color = "black",
    linewidth = 0.3,
    alpha = 0.8,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~ Species,
    ncol = 1,
    strip.position = "top",
    scales = "free_x",
    drop = FALSE
  ) +
  scale_fill_manual(
    name = "Cell Type",
    values = mt.celltype.colors,
    na.translate = FALSE,
    guide = guide_legend(nrow = 1)
  ) +
  labs(
    x = "Cell Types",
    y = "Mitochondrial Percent (%)",
    title = "MT% Distribution"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    strip.background = element_rect(fill = NA, color = NA),
    strip.text = element_text(
      face = "bold",
      size = 12,
      color = "black",
      margin = margin(1, 0, 1, 0, "mm")
    ),
    axis.text.x = element_text(
      color = "black",
      size = 12
    ),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5,
      margin = margin(b = 10)
    ),
    plot.margin = margin(t = 70, r = 5, b = 5, l = 5, unit = "pt"),
    panel.spacing.y = unit(1.2, "lines")
  ) +
  scale_x_discrete(
    drop = FALSE,
    labels = function(x) str_wrap(x, width = 15)
  )

if (length(signif_layers) > 0) {
  p <- p + signif_layers
}

#print(p)

pdf("./new_anno_ATAC/mt.percent/mt.pct.vln.SigTest.MainCelltype2.7panel.pdf",width = 6, height = 15)
plot(p)
dev.off()

write.csv(df_clean, file = '/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/mt.percent/mt.pct.vln.MainCelltype2.data.csv', quote = F)


# 
# p <- ggplot(df, aes(x = Cell_Type, y = mt.percent, fill = Cell_Type)) +
#   geom_violin(
#     trim = TRUE,
#     scale = "width",
#     color = "black",
#     linewidth = 0.3,
#     alpha = 0.8,
#   ) +
#   facet_wrap(
#     ~ Species,
#     ncol = 1,
#   ) +
#   scale_fill_manual(
#     name = "Cell Type",
#     values = mt.celltype.colors,
#     guide = guide_legend(nrow = 1)
#   ) +
#   labs(
#     x = "Cell Types",
#     y = "Mitochondrial Percent (%)",
#     title = "MT% Distribution"
#   ) +
#   theme_bw(base_size = 12) +
#   NoLegend()+
#   theme(
#     panel.grid = element_blank(),
#     panel.border = element_rect(color = "black", linewidth = 0.5),
#     strip.background = element_rect(fill = NA, color = NA),
#     strip.text = element_text(
#       face = "bold",
#       size = 12,
#       color = "black",
#     ),
#     axis.text.x = element_text(
#       color = "black",
#       size = 12#,
#       #angle = 45,
#       #hjust = 1,
#       #vjust = 1
#     ),
#     axis.text.y = element_text(color = "black", size = 12),
#     axis.title = element_text(face = "bold", size = 14),
#     plot.title = element_text(
#       face = "bold",
#       size = 12,
#       hjust = 0.5,
#       margin = margin(b = 10)
#     )
#   )+
# 
# pdf("./new_anno_ATAC/mt.percent/mt.pct.vln.MainCelltype2.7panel.pdf",width = 6, height = 13)
# plot(p)
# dev.off()
# 
# write.csv(df, file = './new_anno_ATAC/mt.percent/mt.pct.vln.MainCelltype2.data.csv')


species.colors <- c("#625D9E", "#CCC9E6","#DCC1DD","#F2E2B2","#C5DEBA", "#AAD4F8", "#5184B2")
names(species.colors) <- c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray")

df_celltype <- df %>%
  group_by(Cell_Type) %>% 
  summarise(
    mt_mean = mean(mt.percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mt_mean)

df_celltype$Cell_Type <- factor(df_celltype$Cell_Type, levels = all_cell_types)

df_subgroup_species <- df %>%
  group_by(Cell_Type, Species) %>%
  summarise(
    subgroup_mean = mean(mt.percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Cell_Type = factor(Cell_Type, levels = levels(df_celltype$Cell_Type)))

df_celltype$Cell_Type <- factor(df_celltype$Cell_Type, levels = celltype_medians$Cell_Type)

p1 <- ggplot() +
  geom_col(
    data = df_celltype,
    aes(x = Cell_Type, y = mt_mean),
    width = 0.6,
    fill = mt.celltype.colors[match(celltype_medians$Cell_Type,names(mt.celltype.colors))],
    color = "black",
    alpha = 0.8
  ) +
  
  geom_point(
    data = df_subgroup_species,
    aes(
      x = Cell_Type,
      y = subgroup_mean,
      fill = Species
    ),
    shape = 21,
    color = "black",
    size = 2,
    position = position_jitter(width = 0.2),
    stroke = 0.5
  ) +
  scale_fill_manual(
    values = species.colors,
    name = "Species"
  ) +
  geom_text(
    data = df_celltype,
    aes(
      x = Cell_Type,
      y = mt_mean,
      label = sprintf("%.1f%%", mt_mean)
    ),
    vjust = 0,
    size = 3,
    color = "black"
  ) +
  labs(
    x = "",
    y = "Mt% Mean"
  ) +
  theme_prism() +
  theme(panel.grid = element_blank())+
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )+
  theme(axis.text.x = element_text(color = "black",size = 10),
        axis.text.y = element_text(color = "black",size = 10))+
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))

pdf("./new_anno_ATAC/mt.percent/mt.pct.top.bar.SigTest.MainCelltype2.pdf",width = 7, height = 2)
plot(p1)
dev.off()



# species.colors <- c("#625D9E", "#CCC9E6","#DCC1DD","#F2E2B2","#C5DEBA", "#AAD4F8", "#5184B2")
# names(species.colors) <- c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray")
# 
# df_celltype <- df %>%
#   group_by(Cell_Type) %>% 
#   summarise(
#     mt_mean = mean(mt.percent, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
# df_celltype$Cell_Type <- factor(df_celltype$Cell_Type, levels = c("B cell","T cell","NK cell","T/NK cell","Monocyte","Other myeloid cell"))
# 
# df_subgroup_species <- df %>%
#   group_by(Cell_Type, Species) %>%
#   summarise(
#     subgroup_mean = mean(mt.percent, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   mutate(Cell_Type = factor(Cell_Type, levels = levels(df_celltype$Cell_Type)))
# 
# p1 <- ggplot() +
#   geom_col(
#     data = df_celltype,
#     aes(x = Cell_Type, y = mt_mean),
#     width = 0.6,
#     fill = "#BAD2DA",
#     alpha = 0.8
#   ) +
#   geom_point(
#     data = df_subgroup_species,
#     aes(
#       x = Cell_Type,
#       y = subgroup_mean,
#     ),
#     size = 2,
#   ) +
#   scale_fill_manual(
#     values = species.colors,
#     name = "Species"
#   ) +
#   geom_text(
#     data = df_celltype,
#     aes(
#       x = Cell_Type,
#       y = mt_mean,
#       label = sprintf("%.1f%%", mt_mean)
#     ),
#     vjust = 0.1,
#     size = 3,
#     color = "black"
#   ) +
#   labs(
#     x = "",
#     y = "Mt% Mean"
#   ) +
#   theme_bw() +
#   theme(panel.grid = element_blank())+
#   theme(
#   )+
#   theme(axis.text.x = element_text(color = "black",size = 10),
#         axis.text.y = element_text(color = "black",size = 10))+
# 
# pdf("./new_anno_ATAC/mt.percent/mt.pct.top.bar.MainCelltype2.pdf",width = 7, height = 2)
# plot(p1)
# dev.off()

df_species <- df %>%
  group_by(Species) %>%
  summarise(
    mt_mean = mean(mt.percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mt_mean)
df_species$Species <- factor(df_species$Species,levels = rev(c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray")))

df_subgroup <- df %>%
  group_by(Species, Cell_Type) %>%
  summarise(
    subgroup_mean = mean(mt.percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Species = factor(Species, levels = rev(c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray"))))


p2 <- ggplot() +
  geom_col(
    data = df_species,
    aes(x = mt_mean, y = Species),
    width = 0.2,
    fill = "#8291B0",#"#9CBAD8",
    alpha = 1
  ) +
  geom_point(
    data = df_subgroup,
    aes(
      x = subgroup_mean,
      y = Species,
      fill = Cell_Type,
    ),
    position = position_jitterdodge(
      jitter.width = 0.3,
      dodge.width = 0.4
    ),
    size = 2,shape = 21,
    color = "grey40"
  ) +
  scale_fill_manual(
    values = mt.celltype.colors,
    name = "Cell Type"
  ) +
  geom_text(
    data = df_species,
    aes(
      x = mt_mean,
      y = Species,
      label = sprintf("%.1f%%", mt_mean)
    ),
    hjust = -0.1,
    size = 3,
    color = "black"
  ) +
  labs(
    x = "Mt% Mean",
    y = "",
    color = "Cell Type"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.2)))+
  theme_bw() +
  theme(panel.grid = element_blank())+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )+
  theme(axis.text.x = element_text(color = "black",size = 10),
        axis.text.y = element_text(color = "black",size = 10))

pdf("./new_anno_ATAC/mt.percent/mt.pct.right.bar.MainCelltype2.pdf",width = 4, height = 13)
plot(p2)
dev.off()

write.csv(df,file="./new_anno_ATAC/mt.percent/mt.pct.vln.MainCelltype2.data.csv")

barcode.list <- list()
for (i in 1:7){
  print(names(meta.list)[i])
  barcode.list[[i]] <- meta.list[[i]]$barcode
  write.table(barcode.list[[i]], file = paste0("/PATH/TO/pbmc.mt/celltype.anno/new_anno_ATAC/seurat_meta/barcodes/",names(meta.list)[i],".barcodes.txt"), row.names = F, quote = F)
}

