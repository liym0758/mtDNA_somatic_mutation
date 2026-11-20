#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
script_arg <- grep('^--file=',commandArgs(),value=TRUE)
if(length(script_arg)!=1L)stop('Run with Rscript')
code <- dirname(normalizePath(sub('^--file=','',script_arg),winslash='/',mustWork=TRUE))
out <- Sys.getenv('MITOEVO_FIG1_OUTPUT_ROOT')
if(!nzchar(out))stop('Set MITOEVO_FIG1_OUTPUT_ROOT')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
key <- read.delim(file.path(code,'global_umap_number_key.tsv'),check.names=FALSE)
stopifnot(identical(key$number,sort(unique(key$number))))
key$y <- rev(seq_len(nrow(key)))
groups <- do.call(rbind,lapply(split(key,key$group),function(d)
 data.frame(group=d$group[1],ymin=min(d$y),ymax=max(d$y),y=mean(d$y))))
cols <- setNames(key$color,key$sub_celltype)
p <- ggplot(key,aes(y=y)) +
 geom_point(aes(x=.10,color=sub_celltype),size=4.8) +
 geom_text(aes(x=.18,label=paste(number,sub_celltype)),hjust=0,size=3.65,color='black') +
 geom_segment(data=groups,aes(x=.95,xend=.95,y=ymin-.3,yend=ymax+.3),
  inherit.aes=FALSE,color='grey60',linewidth=.35) +
 geom_text(data=groups,aes(x=1.0,y=y,label=group),inherit.aes=FALSE,
  hjust=0,size=3.5,fontface='bold',color='black') +
 scale_color_manual(values=cols,guide='none') +
 coord_cartesian(xlim=c(.02,1.5),ylim=c(.2,nrow(key)+.8),clip='off') +
 labs(title='ATAC UMAP cell type legend') +
 theme_void(base_size=11) +
 theme(plot.title=element_text(hjust=.5,face='bold',size=14),
  plot.margin=margin(12,10,10,10,'pt'))
ggsave(file.path(out,'global_umap_number_legend.pdf'),p,width=6,height=8,units='in',device='pdf',useDingbats=FALSE)
write.table(key[,c('number','sub_celltype','color','group')],
 file.path(out,'global_umap_number_key.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
