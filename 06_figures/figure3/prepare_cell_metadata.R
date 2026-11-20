setwd("/PATH/TO/pbmc.mt/celltype.anno/")
my24colors <- c('#7EA475','#C4DAB0','#3B6E4F', '#7EA795', '#699B92', '#274B2F', 
                '#76B2BE','#C3D1DD','#6F8BB0', '#4F418E','#8291B0','#ACBDDE',
                '#DCE572','#CFC436','#E5DBB9', '#D1B56C', '#C6902E', '#AF651E', 
                '#A63376','#853280', '#B3659F','#D490B9', '#AF95C4', '#959492')



#########readRDS#############
human <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/final_ATAC_Human.rds")
mouse <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/final_ATAC_Mouse.rds")
rat <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/Final_ATAC_seurat_Rat.rds")
chicken <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/final_ATAC_Chicken.rds")
cobra <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/final_ATAC_Snake.rds")
sturgeon <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/final_ATAC_XunFish.rds")
stingray <- readRDS("/PATH/TO/external_pbmc/PBMC_IMMU/results/final_ATAC_HongFish.rds")

objlist <- list(mouse=mouse,rat=rat,chicken=chicken,cobra=cobra,sturgeon=sturgeon,stingray=stingray)
species <- c('C57','rat','chicken','cobra','sturgeon','stingray','human')
species2 <- c('mouse','rat','chicken','cobra','sturgeon','stingray','human')
objlist[[7]] <- human
names(objlist)[7] <- "human"

meta.list <- list()
for (i in 1:length(objlist)){
  #objlist[[i]]$barcode <- rownames(objlist[[i]]@meta.data)
  meta.list[[i]] <- objlist[[i]]@meta.data
  write.csv(meta.list[[i]], file = paste("./seurat.meta/",names(objlist)[i],".meta.csv",sep = ''))
}
saveRDS(meta.list, "./seurat.meta/meta.list.RDS")

library(Seurat)
library(tidyverse)

meta.list <- readRDS("./seurat.meta/meta.list.RDS")


mt_counts <- list()
for (i in 1:length(objlist)){
  mt_counts[[i]] <- read.table(file = paste("/PATH/TO/pbmc.mt/mt.bam/bam.QC/orig.bams/",species[i],"/mt_counts.txt",sep = ''), col.names = c("barcode","mt.count")) 
  names(mt_counts)[i] <- species[i]
}


library(tidyverse)
cellranger_metrics <- list()
# cellranger_metrics[[1]] <- as.data.frame(read_csv("/PATH/TO/pbmc.mt/SNP/mouse/cellranger/C57_mouse/outs/singlecell.csv")) %>%
#   dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)

# cellranger_metrics[[2]] <- read_csv("/PATH/TO/pbmc.mt/SNP/rat/cellranger/Rat/outs/singlecell.csv") %>%
#   dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
# cellranger_metrics[[3]] <- read_csv("/PATH/TO/pbmc.mt/SNP/chicken/cellranger/Chicken/outs/singlecell.csv") %>%
#   dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
# cellranger_metrics[[4]] <- read_csv("/PATH/TO/pbmc.mt/SNP/cobra/cellranger/Cobra/outs/singlecell.csv") %>%
#   dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
# cellranger_metrics[[5]] <- read_csv("/PATH/TO/pbmc.mt/SNP/sturgeon/cellranger/Sturgeon_NCBI_mt/outs/singlecell.csv") %>%
#   dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
# cellranger_metrics[[6]] <- read_csv("/PATH/TO/pbmc.mt/SNP/stingray/cellranger/mt.mask/Stingray/outs/singlecell.csv") %>%
#   dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
# cellranger_metrics[[7]] <- read_csv("/PATH/TO/external_data/project/human_PBMC/lib5_NhPBMC/NhPBMC_joint/outs/per_barcode_metrics.csv") %>%
#   dplyr::select(barcode = atac_barcode, total_reads = atac_raw_reads, mt_reads = atac_mitochondrial_reads)

monkey.singlecell <- read.csv("/PATH/TO/external_data/1_Evolution/Data/Monkey/real_monkey-PBMC-C1/outs/singlecell.csv")
monkey.singlecell <- monkey.singlecell[-1,]
sum(monkey.singlecell$mitochondrial)/sum(monkey.singlecell$total)*100
#0.001813262
rm(monkey.singlecell)

cellranger_metrics[[1]] <- as.data.frame(read_csv("/PATH/TO/pbmc.mt/SNP/mouse/cellranger/mt.mask/C57_mouse/outs/singlecell.csv")) %>%
  dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
cellranger_metrics[[2]] <- read_csv("/PATH/TO/pbmc.mt/SNP/rat/cellranger/Rat_mt_mask/outs/singlecell.csv") %>%
  dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
cellranger_metrics[[3]] <- read_csv("/PATH/TO/pbmc.mt/SNP/chicken/cellranger/Chicken_mt_mask/outs/singlecell.csv") %>%
  dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
cellranger_metrics[[4]] <- read_csv("/PATH/TO/pbmc.mt/SNP/cobra/cellranger/Cobra_mt_mask/outs/singlecell.csv") %>%
  dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
cellranger_metrics[[5]] <- read_csv("/PATH/TO/pbmc.mt/SNP/sturgeon/cellranger/Sturgeon_mt_mask/outs/singlecell.csv") %>%
  dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
cellranger_metrics[[6]] <- read_csv("/PATH/TO/pbmc.mt/SNP/stingray/cellranger/mt.mask/Stingray/outs/singlecell.csv") %>%
  dplyr::select(barcode = barcode, total_reads = total, mt_reads = mitochondrial)
cellranger_metrics[[7]] <- read_csv("/PATH/TO/pbmc.mt/SNP/human/cellranger/Human_mt_mask/outs/per_barcode_metrics.csv") %>%
  dplyr::select(barcode = atac_barcode, total_reads = atac_raw_reads, mt_reads = atac_mitochondrial_reads)
for (i in 1:6){
  meta.list[[i]] <- meta.list[[i]] %>%
    left_join(cellranger_metrics[[i]], by = "barcode")
  # meta.list[[i]] <- meta.list[[i]] %>%
  #   left_join(mt_counts[[i]], by = "barcode") 
  meta.list[[i]]$mt.percent <- meta.list[[i]]$mt_reads/meta.list[[i]]$total_reads*100
  meta.list[[i]]$species <- rep(species2[i],nrow(meta.list[[i]]))
}
meta.list[[7]]$total_reads <- meta.list[[7]]$atac_raw_reads
meta.list[[7]]$mt_reads <- meta.list[[7]]$atac_mitochondrial_reads
meta.list[[7]]$mt.percent <- meta.list[[7]]$mt_reads/meta.list[[7]]$total_reads*100
meta.list[[7]]$species <- rep(species2[7],nrow(meta.list[[7]]))


for (i in 1:7){
  rownames(meta.list[[i]]) <- meta.list[[i]]$barcode
}
saveRDS(meta.list, "/PATH/TO/pbmc.mt/celltype.anno/seurat.meta/meta.list.mt.count.mtMask.RDS")

objlist[[7]]$total_reads <- objlist[[7]]$atac_raw_reads
objlist[[7]]$mt_reads <- objlist[[7]]$atac_mitochondrial_reads
objlist[[7]]$mt.percent <- objlist[[7]]$mt_reads/objlist[[7]]$total_reads*100

for (i in 1:7){
  objlist[[i]]$nuclear_reads <- objlist[[i]]$total_reads-objlist[[i]]$mt_reads
  objlist[[i]]$mt.percent2 <- objlist[[i]]$mt_reads/objlist[[i]]$nuclear_reads*100
}


library_summary.list <- list()
for (i in 1:7){
  library_summary.list[[i]] <- data.frame(species = meta.list[[i]]$species, total.reads = meta.list[[i]]$total_reads, mt.reads = meta.list[[i]]$mt_reads, mt.percent = meta.list[[i]]$mt.percent)
}
library_summary <- rbind(library_summary.list[[1]],library_summary.list[[2]],library_summary.list[[3]],library_summary.list[[4]],library_summary.list[[5]],library_summary.list[[6]],library_summary.list[[7]])

summary_df <- library_summary %>%
  drop_na() %>% 
  group_by(species) %>%
  summarise(
    mt.reads = sum(mt.reads),
    total.reads = sum(total.reads)
  ) %>%
  mutate(
    mt.percent = mt.reads / total.reads * 100
  )

summary_df$species <- factor(summary_df$species, levels = c('human','mouse','rat','chicken','cobra','sturgeon','stingray'))
p <- ggplot(summary_df, 
       aes(x = species,  
           y = mt.percent)) +
  geom_col(fill = "#69b3a2", color = "black",size = 0.2) +
  geom_text(aes(label = paste(round(mt.percent,digits = 2),"%")), 
            vjust = -0.5, size = 3.5) +
  labs(
    x = "Species",
    y = "Mitochondrial Percentage (%)",
    title = "Mitochondrial Reads Percentage") +
  theme_classic2() +
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"))

pdf('./mt.percent/preQC.total.mt.pct.BarPlt.pdf',width = 5, height = 5)
plot(p)
dev.off()

p <- ggplot(library_summary, aes(x = species, y = mt.percent, fill = species)) +
  geom_violin(trim = FALSE, alpha = 0.7, adjust = 1.5,color = 'grey20',size = 0.2) + 
  stat_summary(
    fun = mean, 
    geom = "crossbar",
    width = 0.6,
    size = 0.1,
    color = "black",
    show.legend = FALSE
  ) +
  geom_text(
    data = summary_df,
    aes(x = species, y = mt.percent, label = paste(round(mt.percent, 1),'%',sep = '')),
    vjust = -1,
    color = "black",
    size = 3
  ) +
  scale_fill_manual(values = my24colors[2:24]) +
  labs(
    title = "Mitochondrial% of 7 Species",
    x = "Species",
    y = "Mitochondrial%"
  ) +
  theme_classic2() +
  theme(panel.grid = element_blank())+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )+
  theme(axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf('./mt.percent/preQC.mt.percent.7spe.vlnplt.pdf',width = 5, height = 5)
plot(p)
dev.off()

objlist <- readRDS("/PATH/TO/pbmc.mt/celltype.anno/seurat.meta/objlist.orig.RDS")
for (i in 1:7){print(nrow(objlist[[i]]@meta.data)-nrow(meta.list[[i]]))}
for (i in 1:7){
  objlist[[i]]@meta.data <- meta.list[[i]]
}

for (i in 1:7){objlist[[i]]$celltype1 <- objlist[[i]]@active.ident}

rename_mapping <- list(
  "B" = "B cell",
  "B cell" = "B cell",
  "NK" = "NK cell",
  "NK cell" = "NK cell",
  "T" = "T cell",
  "T cell" = "T cell",
  "T/NK" = "T/NK cell",
  "Orthers" = "Others",
  "Otrhers" = "Others",
  "Platelets" = "Platelet",
  "Platelet" = "Platelet",
  "Neutrophil" = "Neutrophil", 
  "Myeloid cell" = "Myeloid cell", 
  "Monocyte" = "Monocyte",
  "Erythrocyte" = "Erythrocyte"
)
for (i in seq_along(objlist)) {
  current_ident <- objlist[[i]]$celltype1
  
  new_ident <- as.character(current_ident)
  mask <- new_ident %in% names(rename_mapping)
  new_ident[mask] <- rename_mapping[new_ident[mask]]
  
  original_levels <- levels(current_ident)
  new_levels <- ifelse(
    original_levels %in% names(rename_mapping),
    rename_mapping[original_levels],
    original_levels
  )
  unique_levels <- unique(new_levels)
  levels_order <- na.omit(match(c("B cell","T cell","T/NK cell","NK cell","Monocyte","Myeloid cell","Neutrophil","Platelet","Erythrocyte","Others"),unique_levels))
  objlist[[i]]$celltype1 <- factor(
    new_ident,
    levels = unique_levels[levels_order]
  )
}

for (i in 1:length(objlist)) {
  print(unique(objlist[[i]]$celltype1))
}

for (i in 1:7){objlist[[i]]$celltype <- objlist[[i]]@active.ident}

rename_mapping <- list(
  "B" = "B cell",
  "B cell" = "B cell",
  "NK" = "NK cell",
  "NK cell" = "NK cell",
  "T" = "T cell",
  "T cell" = "T cell",
  "T/NK" = "T/NK cell",
  "Orthers" = "Others",
  "Otrhers" = "Others",
  "Platelets" = "Platelet",
  "Platelet" = "Platelet",
  "Neutrophil" = "Others", 
  "Myeloid cell" = "Myeloid cell", 
  "Monocyte" = "Myeloid cell",
  "Erythrocyte" = "Erythrocyte"
)
for (i in seq_along(objlist)) {
  current_ident <- objlist[[i]]$celltype
  
  new_ident <- as.character(current_ident)
  mask <- new_ident %in% names(rename_mapping)
  new_ident[mask] <- rename_mapping[new_ident[mask]]
  
  original_levels <- levels(current_ident)
  new_levels <- ifelse(
    original_levels %in% names(rename_mapping),
    rename_mapping[original_levels],
    original_levels
  )
  
  objlist[[i]]$celltype <- factor(
    new_ident,
    levels = c("B cell","T cell","T/NK cell","NK cell","Myeloid cell","Platelet","Erythrocyte","Others")
  )
}

for (i in 1:length(objlist)) {
  print(unique(objlist[[i]]$celltype))
}


unique(c(levels(objlist[[1]]$celltype),levels(objlist[[2]]$celltype),levels(objlist[[3]]$celltype),levels(objlist[[4]]$celltype),levels(objlist[[5]]$celltype),levels(objlist[[6]]$celltype),levels(objlist[[7]]$celltype)))

#celltype.colors <- c('#CCC9E6','#ACBDDE','#7EA475','#C4DAB0','#AF95C4','#8291B0', '#7FB0B4', '#CCE0F5','#DCC1DD', '#959492')
#names(celltype.colors) <- c("B cell","T cell","T/NK cell","NK cell","Monocyte","Myeloid cell","Neutrophil","Platelet","Erythrocyte","Others")

# celltype.colors <- c('#CCE0F5','#C4DAB0','#7EA475','#7FB0B4', '#CCC9E6', '#DCC1DD')
# names(celltype.colors) <- c("B cell","T cell","NK cell","T/NK cell","Myeloid cell","Others")
mt.celltype.colors <- c('#CCC9E6', '#CCE0F5','#9CBAD8','#7FB0B4', '#C4DAB0', '#7EA475','#DCC1DD','#959492')
names(mt.celltype.colors) <- c("B cell","T cell","NK cell","T/NK cell","Myeloid cell","Platelet","Erythrocyte","Others")

p1 <- list()
for (i in 1:7){
  p1[[i]] <- DimPlot(objlist[[i]], label = TRUE, repel = TRUE, reduction = "umap",group.by = "celltype",pt.size = 0.3, cols = celltype.colors)+
    theme(text = element_text(size = 20), strip.text = element_text (size = 20))+ggtitle(paste("UMAP of",str_to_title(species2[i]),"PBMC"))
}
pdf("./celltype/new2.0.celltype.umap.pdf",width=8, height=6)
for (i in 1:7){plot(p1[[i]])}
dev.off()

lymphoid <-  c("B cell" ,"B cells",  "S0X5+ B cells", "B Naive", "B Memory", 
               "CD8+ T", "CD8+ T cells","CD4+ T", "T cell",  "CD4+ T cells", "TCR gamma/delta+ T-cells",  "CD4+ Naive", "CD8+ Memory", "CD4+ Memory", 
               "CD8+ Naive", "T-NK", "NKL","NK","MAIT")
myeloid <-  c("Monocyte", "Neotrophil", "Neutrophil", "Monocytes", 
              "CD16+ Mono", "CD14+ Mono", "cDC", "pDC")
others <-  c("Platelets/Erythrocytes", "Thrombocyte", "Erythrocyte", "Otrhers", 
             "Plsam", "HSP")
predicted.id.colors <- c(my24colors,rev(my36colors)[c(2:5,7:9,6,19,23,22,24:36)])
names(predicted.id.colors) <- c(lymphoid,myeloid,others)
p2 <- list()
for (i in 1:7){
  idx <- na.omit(match(c(lymphoid,myeloid,others),unique(objlist[[i]]$predicted.id)))
  objlist[[i]]$predicted.id <- factor(objlist[[i]]$predicted.id, levels = unique(objlist[[i]]$predicted.id)[idx])
  p2[[i]] <- DimPlot(objlist[[i]], label = TRUE, repel = TRUE, reduction = "umap",group.by = "predicted.id",pt.size = 0.3, cols = predicted.id.colors)+
    theme(text = element_text(size = 20), strip.text = element_text (size = 20))+ggtitle(paste("UMAP of",str_to_title(species2[i]),"PBMC"))
}
pdf("./celltype/sub.celltype.umap.pdf",width=8, height=6)
for (i in 1:7){plot(p2[[i]])}
dev.off()

for (i in 1:7){
  p[[i]] <- VlnPlot(
    object = objlist[[i]],
    features = c("mt.percent"),
    group.by = "celltype",
    pt.size = 0,
    cols=celltype.colors)+ggtitle(paste(str_to_title(paste("Mitochondiral% of cell types in",species2[i]),"PBMC")))
}
pdf("./mt.percent/new.preQC.corse.celltype.mt.mask.pdf",width = 8, height = 5)
for (i in 1:7){plot(p[[i]])}
dev.off()

## 
combined_data <- data.frame()
for (i in 1:7) {
  species_data <- p[[i]][[1]][["data"]]
  species_data$species <- rep(species2[i],nrow(species_data))
  combined_data <- rbind(combined_data, species_data)
}
colnames(combined_data) <- c("mt.percent","Cell_Type","Species")
combined_data$Species <- factor(combined_data$Species, levels = c(species2[7],species2[1:6]))
#combined_data$Cell_Type <- factor(combined_data$Cell_Type,levels = c("B cell","T cell","T/NK cell","NK cell","Monocyte","Myeloid cell","Neutrophil","Platelet","Erythrocyte","Others"))
#combined_data$Cell_Type <- factor(combined_data$Cell_Type,levels = c("B cell","T cell","NK cell","T/NK cell","Myeloid cell","Others"))
combined_data$Cell_Type <- factor(combined_data$Cell_Type,levels = c("B cell","T cell","NK cell","T/NK cell","Myeloid cell","Platelet","Erythrocyte","Others"))

vln1 <- ggplot(combined_data, aes(x = Species, y = mt.percent, fill = Cell_Type)) +
  geom_violin(trim = FALSE, scale = "width") +
  labs(title = "MT% Distribution by Species",
       x = "Species",
       y = "MT%") +
  scale_fill_manual(values = celltype.colors) +
  theme_classic2() +
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )+
  theme(axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


species.colors <- c('#CCC9E6','#ACBDDE','#7EA475','#C4DAB0','#AF95C4','#8291B0', '#7FB0B4', '#CCE0F5','#DCC1DD', '#959492')
names(species.colors) <- c('human','mouse','rat','chicken','cobra','sturgeon','stingray')
vln2 <- ggplot(combined_data, aes(x = Cell_Type, y = mt.percent, fill = Species)) +
  geom_violin(trim = FALSE, scale = "width") +
  labs(title = "MT% Distribution by Cell Types",
       x = "Cell Types",
       y = "MT%") +
  scale_fill_manual(values = species.colors) +
  theme_classic2() +
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )+
  theme(axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf("./mt.percent/MT.pct.merged.celltype.vlnplt.pdf", width = 18, height = 5)
plot(vln1)
plot(vln2)
dev.off()

library(ggplot2)
library(tidyr)

combined_data <- data.frame()
for (i in 1:7) {
  species_data <- p[[i]][[1]][["data"]]
  species_data$species <- rep(species2[i],nrow(species_data))
  combined_data <- rbind(combined_data, species_data)
}
colnames(combined_data) <- c("mt.percent","Cell_Type","Species")
combined_data$Species <- factor(combined_data$Species, levels = c(species2[7],species2[1:6]))

combined_data$Cell_Type <- factor(combined_data$Cell_Type,levels = c("B cell","T cell","NK cell","T/NK cell","Myeloid cell","Platelet","Erythrocyte","Others"))

all_species <- c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray")
all_cell_types <- c("B cell","T cell","NK cell","T/NK cell","Myeloid cell","Platelet","Erythrocyte","Others")
combined_data$Species <- str_to_title(combined_data$Species)


df <- combined_data %>%
  mutate(
    Species = factor(Species, levels = all_species),
    Cell_Type = factor(Cell_Type, levels = all_cell_types)
  ) %>%
  complete(Species, Cell_Type)



p <- ggplot(df, aes(x = Cell_Type, y = mt.percent, fill = Cell_Type)) +
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
  NoLegend()+
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(
      face = "bold",
      size = 10,
      color = "black",
      margin = margin(1, 0, 1, 0, "mm")
    ),
    axis.text.x = element_text(
      color = "black",
      size = 10#,
      #angle = 45,
      #hjust = 1,
      #vjust = 1
    ),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title = element_text(face = "bold", size = 12),
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5,
      margin = margin(b = 10)
    )
  )
pdf("./mt.percent/mt.mask/mt.pct.vln.7panel.mtMask.pdf",width = 6, height = 12)
plot(p)
dev.off()

write.csv(df, file = './mt.percent/mt.mask/mt.pct.vln.mtMask.data.csv')

##################################################
species.colors <- c("#625D9E", "#CCC9E6","#DCC1DD","#F2E2B2","#C5DEBA", "#AAD4F8", "#5184B2")
names(species.colors) <- c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray")

df_celltype <- df %>%
  group_by(Cell_Type) %>% 
  summarise(
    mt_mean = mean(mt.percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mt_mean)
df_celltype$Cell_Type <- factor(df_celltype$Cell_Type, levels = c("B cell","T cell","NK cell","T/NK cell","Monocyte","Myeloid cell","Neutrophil","Platelet","Erythrocyte","Others"))

df_subgroup_species <- df %>%
  group_by(Cell_Type, Species) %>%
  summarise(
    subgroup_mean = mean(mt.percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Cell_Type = factor(Cell_Type, levels = levels(df_celltype$Cell_Type)))

p1 <- ggplot() +
  geom_col(
    data = df_celltype,
    aes(x = Cell_Type, y = mt_mean),
    width = 0.6,
    fill = "#BAD2DA",
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
    color = "grey40",
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
    vjust = 0.1,
    size = 3,
    color = "black"
  ) +
  # stat_pvalue_manual(
  #   data = stat_data1,
  #   label = "p.adj.signif",
  #   tip.length = 0.01,
  #   bracket.size = 0.2,
  #   label.size = 2,
  # ) +
  labs(
    x = "",
    y = "Mt% Mean"
  ) +
  shared_theme

pdf("./mt.percent/mt.mask/mt.pct.top.bar2.mtMask.pdf",width = 6, height = 2)
plot(p1)
dev.off()

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
    fill = "#9CBAD8",
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
  shared_theme

pdf("./mt.percent/mt.mask/mt.pct.right.bar.mtMask.pdf",width = 3.5, height = 12)
plot(p2)
dev.off()

write.csv(df,file="./mt.percent/mt.mask/mt.pct.vln.data.csv")


###################################
for (i in 1:7){
  p[[i]] <- VlnPlot(
    object = objlist[[i]],
    features = c("mt.percent"),
    group.by = "predicted.id",
    pt.size = 0,
    cols=predicted.id.colors)+ggtitle(paste(species2[i],"- mt% of each cell type"))
}
pdf("./mt.percent/preQC.celltype.mt.percent.pdf",width = 8, height = 5)
for (i in 1:7){plot(p[[i]])}
dev.off()

combined_data <- data.frame()
for (i in 1:7) {
  species_data <- p[[i]][[1]][["data"]]
  species_data$species <- rep(species2[i],nrow(species_data))
  combined_data <- rbind(combined_data, species_data)
}
colnames(combined_data) <- c("mt.percent","Cell_Type","Species")
combined_data$Species <- factor(combined_data$Species, levels = c(species2[7],species2[1:6]))
combined_data$Cell_Type <- factor(combined_data$Cell_Type,levels = c(lymphoid,myeloid,others))

vln1 <- ggplot(combined_data, aes(x = Species, y = mt.percent, fill = Cell_Type)) +
  geom_violin(trim = FALSE, scale = "width") +
  labs(title = "MT% Distribution by Species",
       x = "Species",
       y = "MT%") +
  scale_fill_manual(values = predicted.id.colors) +
  theme_classic2() +
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )+
  theme(axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))


vln2 <- ggplot(combined_data, aes(x = Cell_Type, y = mt.percent, fill = Species)) +
  geom_violin(trim = FALSE, scale = "width") +
  labs(title = "MT% Distribution by Cell Types",
       x = "Cell Types",
       y = "MT%") +
  scale_fill_manual(values = species.colors) +
  theme_classic2() +
  theme(panel.grid = element_blank())+
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )+
  theme(axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf("./mt.percent/MT.pct.merged.sub.celltype.vlnplt.pdf", width = 20, height = 5)
plot(vln1)
plot(vln2)
dev.off()
saveRDS(objlist, "/PATH/TO/pbmc.mt/celltype.anno/seurat.meta/objlist.RDS")




##########Seurat QC#############
objlist.filter <- list()
for (i in 1:7){
  objlist.filter[[i]] <- subset(
    x = objlist[[i]],
    subset = mt.percent > 5
  )
}
saveRDS(objlist.filter, "/PATH/TO/pbmc.mt/celltype.anno/seurat.meta/objlist.filter.RDS")
cell.num <- data.frame(
  preQC = c(nrow(objlist[[1]]@meta.data),nrow(objlist[[2]]@meta.data),nrow(objlist[[3]]@meta.data),nrow(objlist[[4]]@meta.data),nrow(objlist[[5]]@meta.data),nrow(objlist[[6]]@meta.data),nrow(objlist[[7]]@meta.data)),
  afterQC = c(nrow(objlist.filter[[1]]@meta.data),nrow(objlist.filter[[2]]@meta.data),nrow(objlist.filter[[3]]@meta.data),nrow(objlist.filter[[4]]@meta.data),nrow(objlist.filter[[5]]@meta.data),nrow(objlist.filter[[6]]@meta.data),nrow(objlist.filter[[7]]@meta.data)),
  species = species2
)
df_long <- pivot_longer(cell.num, cols = c(preQC, afterQC), names_to = "QC_type", values_to = "value")
df_long$QC_type <- factor(df_long$QC_type,levels <- c('preQC','afterQC'))
df_long$species <- factor(df_long$species, levels = species2[c(7,1:6)])
p <- ggplot(df_long, aes(x = species, y = value, fill = QC_type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black",size = 0.2) +
  scale_fill_manual(values = my24colors[3:24]) +
  labs(title = "Cell number before and after QC (mt% > 5)",
       x = "Species",
       y = "# of cells",
       fill = "QC Type") +
  geom_text(aes(label = value), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3, color = "black") +
  theme_classic2() +
  theme(panel.grid = element_blank())+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )+
  theme(axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12))
  
pdf('./mt.percent/afterQC.cell.num.pdf', width = 8, height = 5)
plot(p)
dev.off()

library_summary.list <- list()
for (i in 1:7){
  library_summary.list[[i]] <- data.frame(
    species = objlist.filter[[i]]$species, 
    total.reads = objlist.filter[[i]]$total_reads, 
    mt.reads = objlist.filter[[i]]$mt_reads, 
    mt.percent = objlist.filter[[i]]$mt.percent)
}
library_summary <- rbind(library_summary.list[[1]],library_summary.list[[2]],library_summary.list[[3]],library_summary.list[[4]],library_summary.list[[5]],library_summary.list[[6]],library_summary.list[[7]])

p <- ggplot(library_summary, aes(x = species, y = mt.percent, fill = species)) +
  geom_violin(trim = FALSE, alpha = 0.7, adjust = 1.5,color = 'grey20',size = 0.2) + 
  stat_summary(
    fun = mean, 
    geom = "crossbar",
    width = 0.6,
    size = 0.1,
    color = "black",
    show.legend = FALSE
  ) +
  geom_text(
    data = summary_df,
    aes(x = species, y = mt.percent, label = paste(round(mt.percent, 1),'%',sep = '')),
    vjust = -1.3,
    color = "black",
    size = 3
  ) +
  scale_fill_manual(values = my24colors[2:24]) +
  labs(
    title = "Mitochondrial% of 7 Species (afterQC)",
    x = "Species",
    y = "Mitochondrial%"
  ) +
  theme_classic2() +
  theme(panel.grid = element_blank())+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )+
  theme(axis.text.x = element_text(color = "black",size = 12),
        axis.text.y = element_text(color = "black",size = 12))+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf('./mt.percent/afterQC.mt.percent.7spe.vlnplt.pdf',width = 5, height = 4)
plot(p)
dev.off()

for (i in 1:7){
  p[[i]] <- VlnPlot(
    object = objlist.filter[[i]],
    features = c("mt.percent"),
    group.by = "predicted.id",
    pt.size = 0,
    cols=my24colors[2:24])+ggtitle(paste(species2[i],"- mt% of each cell type"))
}
pdf("./mt.percent/afterQC.celltype.mt.percent.pdf",width = 8, height = 5)
for (i in 1:7){plot(p[[i]])}
dev.off()
saveRDS(objlist.filter, "/PATH/TO/pbmc.mt/celltype.anno/seurat.meta/objlist.filter.RDS")


###########save##############
saveRDS(objlist, "/PATH/TO/pbmc.mt/celltype.anno/seurat.meta/objlist.RDS")

