suppressPackageStartupMessages({library(data.table);library(dplyr)})
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==2)
task<-normalizePath(args[1],winslash='/',mustWork=TRUE)
prior<-normalizePath(args[2],winslash='/',mustWork=TRUE)
vafdir<-'/PATH/TO/local_simulation/server_file/VAF_histogram/fixedDepth_VAF-new_human-Reorder'
species<-c('human','mouse','rat','chicken','cobra','sturgeon','stingray')
vrows<-list();hrows<-list();qrows<-list()
for(sp in species){
 x<-as.data.table(readRDS(file.path(vafdir,paste0(sp,'_Reorder_celltype_VAF_data_260128.RDS')))[[1]])
 wanted<-fread(file.path(prior,'remote_outputs','tables',paste0(sp,'_final_variants.tsv')))$variant
 x<-x[mutation %in% wanted & is.finite(frequency) & frequency>0]
 stopifnot(!anyDuplicated(x[,.(mutation,barcode)]),all(x$frequency<=1),uniqueN(x$mutation)==length(wanted))
 v<-x[,.(mean_vaf=mean(frequency),n_carrier_observations=.N,n_distinct_cells=uniqueN(barcode),pct_conf_detected_cell=unique(pct_conf_detected_cell)[1],position=as.integer(pos[1])),by=mutation]
 v[,species:=sp];vrows[[sp]]<-v
 x[,bin:=cut(frequency,breaks=seq(0,1,.1),include.lowest=FALSE,right=TRUE,labels=sprintf('(%.1f-%.1f]',seq(0,.9,.1),seq(.1,1,.1)))]
 h<-x[,.(observations=.N,unique_variants=uniqueN(mutation),unique_cells=uniqueN(barcode)),by=.(celltype=celltype_resolution,bin)]
 h[,`:=`(species=sp,celltype_total=sum(observations)),by=celltype]
  ctvar<-x[,.(celltype_variant_total=uniqueN(mutation)),by=.(celltype=celltype_resolution)]
  h<-merge(h,ctvar,by='celltype',all.x=TRUE,sort=FALSE)
 h[,within_celltype_percent:=100*observations/celltype_total]
 hrows[[sp]]<-h
 qrows[[sp]]<-data.table(species=sp,n_observations=nrow(x),n_variants=nrow(v),n_zero=sum(x$frequency==0),n_unique_cell_mutation=uniqueN(x,by=c('mutation','barcode')))
}
fwrite(rbindlist(vrows),file.path(task,'Tables','figure3_lollipop_original_variant_summary.tsv'),sep='\t')
fwrite(rbindlist(hrows),file.path(task,'Tables','figure2e_original_histogram_counts.tsv'),sep='\t')
fwrite(rbindlist(qrows),file.path(task,'Tables','original_vaf_input_validation.tsv'),sep='\t')
print(rbindlist(qrows))
