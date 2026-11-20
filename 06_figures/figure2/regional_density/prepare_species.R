#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('usage: prepare_species.R species job_dir')
sp<-args[1];job<-args[2];root<-'/PATH/TO/pbmc.mt'
types<-c('T cell','B cell','NK cell','T/NK cell','Monocyte')
meta_path<-if(sp=='human')file.path(root,'SNP/human/clean_bam/merge_bam/human.PBMConly.sample.barcode.celltype.mtPercent.withArchR.csv') else file.path(root,'celltype.anno/new_anno_ATAC/seurat_meta',paste0(sp,'.meta.csv'))
meta<-fread(meta_path)
stopifnot(all(c('barcode','main_celltype2')%in%names(meta)))
m<-meta[main_celltype2%chin%types,.(barcode,main_celltype2)]
if(anyDuplicated(m$barcode))stop('duplicate metadata barcode: ',sp)
idx<-setNames(0:4,types);m[,celltype_index:=unname(idx[main_celltype2])]
# The Human VAF Histogram metadata already carries the PBMC sample prefix,
# so its barcode matches merged.coverage.txt.gz exactly.
m[,coverage_barcode:=barcode]
fwrite(m[,.(barcode,coverage_barcode,main_celltype2,celltype_index)],
       file.path(job,paste0(sp,'.metadata_mapping.tsv')),sep='\t')
fwrite(m[,.(coverage_barcode,celltype_index)],
       file.path(job,paste0(sp,'.barcode_map.tsv')),sep='\t',col.names=FALSE)
counts<-m[,.(metadata_cells=.N),by=.(celltype=main_celltype2)]
counts<-merge(data.table(celltype=types),counts,by='celltype',all.x=TRUE)
counts[is.na(metadata_cells),metadata_cells:=0L]
fwrite(counts,file.path(job,paste0(sp,'.metadata_counts.tsv')),sep='\t')
cat(sp,' metadata_cells=',nrow(m),'\n')
