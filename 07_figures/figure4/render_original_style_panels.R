suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(ggplot2);library(stringr)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
task <- normalizePath(args[1],winslash='/',mustWork=TRUE)
prior <- normalizePath(args[2],winslash='/',mustWork=TRUE)
out <- file.path(task,'Supplementary_candidates','figure4_original_style_species')
dir.create(out,showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(task,'Tables'),showWarnings=FALSE,recursive=TRUE)

# The 56 rows and the original P values are read unchanged from the legacy output.
# Only display text changes: existing stars remain, P>=0.05 becomes ns, empty bins stay blank.
dot <- readr::read_tsv(file.path(prior,'remote_outputs','tables','old_AllSpecies_dot_dnds_results_withP.tsv'),show_col_types=FALSE)
stopifnot(nrow(dot)==56L,all(dot$total_mut==dot$nonsynonymous+dot$synonymous))
dot <- dot %>% mutate(star_display=if_else(is.na(p_binom),'',if_else(p_binom>=0.05,'ns',star)))
readr::write_tsv(dot,file.path(task,'Tables','figure4_56_original_results_with_display_labels.tsv'))
bar <- readr::read_tsv(file.path(prior,'remote_outputs','tables','old_all_species_sitelevel_dnds_noBin_withP.tsv'),show_col_types=FALSE)
stopifnot(nrow(bar)==14L)
readr::write_tsv(bar,file.path(task,'Tables','figure4_14_original_overall_results.tsv'))

# Adapted directly from create_dnds_plot() in the original dNdS_withP.R.
create_dnds_plot <- function(df,celltype,species) {
  df$bin <- factor(df$bin,levels=c('(0,0.1]','(0.1,0.3]','(0.3,0.5]','(0.5,1]'))
  df_total <- df %>% summarise(total_nonsyn=sum(nonsynonymous,na.rm=TRUE),total_syn=sum(synonymous,na.rm=TRUE),total_NS_ratio=ifelse(total_syn==0,NA,total_nonsyn/total_syn),total_dnds=total_NS_ratio/unique(NS_ratio_expected))
  df <- df %>% mutate(total_mut_plot=pmin(total_mut,100),legend_key=case_when(celltype_group=='T_B_NK'~'T_B_NK',celltype_group=='Mono'~'Mono'),star_y=pmin(dnds+0.18,2.95))
  ggplot(df,aes(x=bin,y=dnds,group=1))+
    geom_point(aes(size=total_mut_plot,fill=legend_key),shape=21,color='black',stroke=1.1)+
    geom_text(aes(y=star_y,label=star_display),size=5,color='black',na.rm=TRUE)+
    scale_fill_manual(values=c('T_B_NK'='#8E90B9','Mono'='#7FB0B4'),breaks=c('T_B_NK','Mono'),limits=c('T_B_NK','Mono'),drop=FALSE)+
    scale_size_continuous(limits=c(0,100),range=c(2.2,7.0),breaks=c(20,40,60,80,100),labels=c('<=20','40','60','80','>=100'))+
    geom_hline(yintercept=1,linetype='dashed',color='grey40')+
    geom_text(data=df_total,aes(x=2,y=3.5,label=sprintf('total dN/dS = %.2f',total_dnds)),inherit.aes=FALSE,size=3.5)+
    scale_y_continuous(limits=c(0,3),breaks=seq(0,3,0.5))+
    theme_bw(base_size=12)+
    labs(title=paste0(species,' - ',celltype),x='Mean VAF bin (per mutation)',y='dN/dS',fill='Group',size='# mutations in bin')+
    theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
}
species <- c('human','mouse','rat','chicken','cobra','sturgeon','stingray')
for(sp in species){
  d <- dot %>% filter(species==sp)
  p1 <- create_dnds_plot(d %>% filter(celltype_group=='T_B_NK'),'T_B_NK',sp)+theme(legend.position='none')
  p2 <- create_dnds_plot(d %>% filter(celltype_group=='Mono'),'Mono',sp)
  ggsave(file.path(out,paste0(sp,'_T_B_NK.pdf')),p1,width=4,height=3)
  ggsave(file.path(out,paste0(sp,'_Mono.pdf')),p2,width=4,height=3)
}

# Original VAF bin counts and proportions; only the y-unit and label change.
vaf <- readr::read_tsv(file.path(prior,'remote_outputs','tables','fig3_vaf_bins.tsv'),show_col_types=FALSE) %>% filter(filter_mode=='legacy_all_positive')
stopifnot(all(abs(vaf$percent-100*vaf$numerator/vaf$denominator)<1e-9))
readr::write_tsv(vaf,file.path(task,'Tables','figure3_vaf_all_positive_bins.tsv'))
vout <- file.path(task,'Supplementary_candidates','figure3_vaf_original_style')
dir.create(vout,showWarnings=FALSE,recursive=TRUE)
bin_levels <- levels(cut(0.05,breaks=seq(0,1,.1),include.lowest=TRUE))
type_levels <- c('synonymous','nonsynonymous','rRNA','tRNA','intergenic')
for(sp in species){
  x <- vaf %>% filter(species==sp,mutation_type %in% type_levels) %>% mutate(vaf_bin=factor(vaf_bin,levels=bin_levels),mutation_type=factor(mutation_type,levels=type_levels)) %>% tidyr::complete(mutation_type,vaf_bin,fill=list(numerator=0L,percent=0))
  stopifnot(all(x$percent>=0),all(abs(x %>% group_by(mutation_type) %>% summarise(z=sum(percent)) %>% pull(z)-100)<1e-7))
  p <- ggplot(x,aes(x=vaf_bin,y=percent,color=mutation_type,group=mutation_type))+
    geom_line(linewidth=1)+geom_point(size=1.5)+theme_bw(base_size=12)+
    theme(panel.grid=element_blank(),text=element_text(color='black'),axis.text.x=element_text(color='black',angle=45,hjust=1),axis.text.y=element_text(color='black'),axis.title.x=element_text(color='black'),axis.title.y=element_text(color='black'),plot.title=element_text(color='black',hjust=.5),legend.title=element_text(color='black'),legend.text=element_text(color='black'))+
    labs(title=str_to_title(sp),x='Variant Allele Frequency (VAF)',y='Proportion of carrier cell–variant observations (%)',color='Type')+
    scale_color_manual(values=c('synonymous'='#E78AC3','nonsynonymous'='#FC8D62','rRNA'='#8DA0CB','tRNA'='#A6D854','intergenic'='#66C2A5'))
  ggsave(file.path(vout,paste0(sp,'_VAF_all_positive_axis_corrected.pdf')),p,width=6,height=3.5)
}

# Assemble the equal-width Figure 4 from the same validated tables and panel function.
suppressPackageStartupMessages({library(cowplot); library(png)})
out_dir <- Sys.getenv('MITOEVO_FIG4_OUTPUT_ROOT', unset=file.path(task, 'Figure4_equal_width'))
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
colors <- c(T_B_NK='#8E90B9', Mono='#7FB0B4')
# Repeated sidecar legend for each species pair. This keeps both lineage plotting
# areas equal while retaining a local key beside every pair. Text is >=12 pt.
legend_panel <- ggplot() +
  annotate("text", x=.02, y=.98, label="# mutations", hjust=0, size=4.3, color="black") +
  annotate("text", x=.02, y=.91, label="in bin", hjust=0, size=4.3, color="black") +
  annotate("point", x=.15, y=c(.81,.72,.63,.54,.45), size=c(2.2,3.4,4.6,5.8,7.0),
           shape=21, fill="white", color="black", stroke=.8) +
  annotate("text", x=.34, y=c(.81,.72,.63,.54,.45),
           label=c("<=20","40","60","80",">=100"), hjust=0, size=4.3, color="black") +
  annotate("text", x=.02, y=.34, label="Group", hjust=0, size=4.3, color="black") +
  annotate("point", x=.15, y=c(.22,.11), size=4.3, shape=21,
           fill=c(colors["T_B_NK"],colors["Mono"]), color="black", stroke=.8) +
  annotate("text", x=.34, y=c(.22,.11), label=c("T_B_NK","Mono"), hjust=0, size=4.3, color="black") +
  coord_cartesian(xlim=c(0,1),ylim=c(0,1),clip="off") +
  theme_void() + theme(text=element_text(color="black"))

panels <- setNames(vector("list", length(species)*2L), as.vector(rbind(paste0(species,"_T_B_NK"),paste0(species,"_Mono"))))
for (sp in species) {
  d <- filter(dot, species==sp)
  panels[[paste0(sp,"_T_B_NK")]] <- create_dnds_plot(filter(d,celltype_group=="T_B_NK"),"T_B_NK",sp) + theme(legend.position="none", axis.text=element_text(size=12,color="black"), axis.title=element_text(size=13,color="black"), plot.title=element_text(size=14,color="black",hjust=.5))
  panels[[paste0(sp,"_Mono")]] <- create_dnds_plot(filter(d,celltype_group=="Mono"),"Mono",sp) + theme(legend.position="none", axis.text=element_text(size=12,color="black"), axis.title=element_text(size=13,color="black"), plot.title=element_text(size=14,color="black",hjust=.5))
}

# Replot overall panel h from the already validated 14-row summary table. All
# values, bars, colors, legend categories, and original significance labels are retained.
overall_df <- bar %>%
  mutate(
    species=factor(species,levels=c("Human","Mouse","Rat","Chicken","Cobra","Sturgeon","Stingray")),
    celltype_group=factor(celltype_group,levels=c("T_B_NK","Mono")),
    legend_key=factor(legend_key,levels=c("T_B_NK (>=20 muts)","Mono (>=20 muts)","<20 muts"))
  )
stopifnot(nrow(overall_df)==14L)
dodge <- position_dodge(width=.75)
overall_panel <- ggplot(overall_df,aes(x=species,y=dnds)) +
  geom_col(aes(group=celltype_group,fill=legend_key),position=dodge,width=.65,color="black") +
  geom_text(aes(group=celltype_group,label=star_vs_bg,y=dnds+.1),position=dodge,size=5,color="black") +
  scale_fill_manual(values=c("T_B_NK (>=20 muts)"="#8E90B9","Mono (>=20 muts)"="#7FB0B4","<20 muts"="white"),
                    guide=guide_legend(override.aes=list(color="black"))) +
  geom_hline(yintercept=1,linetype="dashed",color="grey40") +
  scale_y_continuous(limits=c(0,1.8)) +
  labs(title="Overall mtDNA dN/dS",x="Species",y="dN/dS",fill=NULL) +
  theme_bw(base_size=14) +
  theme(
    text=element_text(color="black"),panel.grid=element_blank(),
    axis.text=element_text(size=12,color="black"),
    axis.title=element_text(size=14,color="black"),
    plot.title=element_text(size=16,face="bold",hjust=.5,color="black"),
    legend.text=element_text(size=12,color="black"),
    legend.title=element_text(size=12,color="black")
  )

# Canvas follows the original 1300 x 1240 page and panel ordering. Both lineage
# plots occupy identical 270 x 202.5 pt boxes; tick labels are 12 pt and all text is black.
page_w <- 1300; page_h <- 1240
canvas <- ggdraw() + theme(plot.background=element_rect(fill="white",color=NA))
positions <- list(c(0,3),c(0,2),c(0,1),c(0,0),c(1,3),c(1,2),c(1,1))
for (i in seq_along(species)) {
  col <- positions[[i]][1]; row <- positions[[i]][2]
  x0 <- if (col==0) 25 else 665
  y0 <- 55 + row*290
  plot_y <- y0 + 28.75
  plot_h <- 202.5
  canvas <- canvas +
    draw_plot(panels[[paste0(species[i],"_T_B_NK")]], x=x0/page_w, y=plot_y/page_h, width=270/page_w, height=plot_h/page_h) +
    draw_plot(panels[[paste0(species[i],"_Mono")]], x=(x0+278)/page_w, y=plot_y/page_h, width=270/page_w, height=plot_h/page_h) +
    draw_plot(legend_panel, x=(x0+548)/page_w, y=(y0+20)/page_h, width=82/page_w, height=230/page_h) +
    draw_label(paste0(letters[i],"  ",tools::toTitleCase(species[i])), x=(x0+3)/page_w, y=(y0+261)/page_h,
               hjust=0, vjust=0, size=12, color="black", fontface="bold")
}
# Overall panel h keeps the original grouped-bar structure and p-value labels, with legible black text.
canvas <- canvas + draw_plot(overall_panel,x=665/page_w,y=71.5/page_h,width=605/page_w,height=226.875/page_h)
canvas <- canvas + draw_label("h  Overall dN/dS", x=668/page_w, y=316/page_h,
                               hjust=0, vjust=0, size=12, color="black", fontface="bold")

pdf_out <- file.path(out_dir,"Figure4_equal_width_candidate.pdf")
png_out <- file.path(out_dir,"Figure4_equal_width_candidate.png")
ggsave(pdf_out, canvas, width=page_w/72, height=page_h/72, units="in", device=cairo_pdf, bg="white")
ggsave(png_out, canvas, width=page_w/72, height=page_h/72, units="in", dpi=160, bg="white")
cat("PDF:",pdf_out,"\nPNG:",png_out,"\nRows:",nrow(dot),"\n")

