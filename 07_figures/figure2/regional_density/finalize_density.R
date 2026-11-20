#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
root<-'/PATH/TO/pbmc.mt'
job<-file.path(root,'code/regional_density_final_five_celltypes_20261107_v1')
species<-c('human','mouse','rat','chicken','cobra','sturgeon','stingray')
regions<-c('CDS','rRNA','tRNA','Noncoding')
types<-c('T cell','B cell','NK cell','T/NK cell','Monocyte')
lookup<-fread(file.path(job,'Inputs/position_region_original_mask.tsv'))
lookup[,`:=`(position=as.integer(position),legacy_low_coverage_mask=as.logical(legacy_low_coverage_mask))]
stopifnot(setequal(unique(lookup$species),species),!anyDuplicated(lookup[,.(species,position)]),
          all(lookup$region%chin%regions))
lollipop<-fread(file.path(job,'Inputs/original_lollipop_events.tsv'))
lollipop[,position:=as.integer(position)]
all_sites<-list();all_events<-list();coverage_counts<-list();mask_removed<-list();recon<-list();event_recon<-list()
for(sp in species){
  message('Finalize ',sp)
  map<-fread(file.path(job,paste0(sp,'.metadata_mapping.tsv')))
  matched<-fread(file.path(job,paste0(sp,'.matched_barcodes.tsv')))
  matched[,celltype:=types[celltype_index+1L]]
  if(sp=='human'&&nrow(matched)==0L)stop('Human coverage barcode intersection is empty; verify prefixed metadata')
  matched[,metadata_barcode:=coverage_barcode]
  if(anyDuplicated(matched$metadata_barcode))stop('multiple coverage barcodes map to one metadata cell: ',sp)
  matched_map<-merge(matched,map[,.(barcode,main_celltype2)],by.x=c('metadata_barcode','celltype'),
                     by.y=c('barcode','main_celltype2'),all=FALSE)
  if(nrow(matched_map)!=nrow(matched))stop('coverage barcode failed to map back to metadata: ',sp)
  cc<-fread(file.path(job,paste0(sp,'.metadata_counts.tsv')))
  mc<-matched_map[,.(coverage_matched_cells=.N),by=.(celltype)]
  cc<-merge(cc,mc,by='celltype',all.x=TRUE);cc[is.na(coverage_matched_cells),coverage_matched_cells:=0L]
  cc[,`:=`(species=sp,unmatched_metadata_cells=metadata_cells-coverage_matched_cells)]
  coverage_counts[[sp]]<-cc

  covpos<-fread(file.path(job,paste0(sp,'.positive_positions.tsv')))
  covpos[,position:=as.integer(position)]
  lk<-lookup[species==sp]
  unmasked<-lk[legacy_low_coverage_mask==FALSE]
  positive_positions<-covpos[n_positive_coverage_records>0,unique(position)]
  Ltab<-unmasked[position%chin%positive_positions,.N,by=region]
  setnames(Ltab,'N','L')
  Ltab<-merge(data.table(region=regions),Ltab,by='region',all.x=TRUE)
  Ltab[is.na(L),L:=0L]

  vf<-file.path(root,'SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human',paste0(sp,'_detail_filtered_df_260128.tsv'))
  vaf<-fread(vf,select=c('mutation','pos','barcode','frequency','main_celltype2'))
  vaf<-vaf[main_celltype2%chin%types & !is.na(frequency) & frequency>0]
  vaf[,metadata_barcode:=barcode]
  matched_meta<-matched_map[,.(metadata_barcode,main_celltype2=celltype)]
  vaf<-merge(vaf,matched_meta,by=c('metadata_barcode','main_celltype2'),all=FALSE)
  vaf[, c("position","species") := .(as.integer(pos),sp)]
  event<-vaf[,.(n_positive_observations=.N,n_positive_cells=uniqueN(barcode),
                carrying_cell_types=paste(sort(unique(main_celltype2)),collapse=';')),
             by=.(species,position,mutation)]
  event<-merge(event,lk[,.(position,region,legacy_low_coverage_mask)],by='position',all.x=TRUE)
  if(anyNA(event$region))stop('VAF mutation has no original region annotation: ',sp)
  masked<-event[legacy_low_coverage_mask==TRUE]
  if(nrow(masked)){mask_removed[[sp]]<-masked}
  event<-event[legacy_low_coverage_mask==FALSE]
  all_events[[sp]]<-event
  site<-event[,.(N=uniqueN(position)),by=.(species,region)]
  count<-merge(data.table(species=sp,region=regions),site,by=c('species','region'),all.x=TRUE)
  count[is.na(N),N:=0L]
  count<-merge(count,Ltab,by='region',all.x=TRUE)
  count[,density_sites_per_kb:=ifelse(L>0,N/L*1000,NA_real_)]
  all_sites[[sp]]<-count

  lp<-unique(lollipop[species==sp,.(species,position)])
  ns<-unique(event[,.(species,position)])
  both<-merge(lp,ns,by=c('species','position'))
  lonly<-lp[!ns,on=.(species,position)];lonly[,status:='lollipop_only']
  vonly<-ns[!lp,on=.(species,position)];vonly[,status:='vaf_five_celltype_only']
  diff<-rbindlist(list(lonly,vonly),fill=TRUE)
  if(nrow(diff))recon[[sp]]<-diff
  lr<-data.table(species=sp,lollipop_events=uniqueN(lollipop[species==sp,mutation]),
                 lollipop_positions=nrow(lp),vaf_five_type_events=uniqueN(event$mutation),
                 vaf_five_type_positions=uniqueN(ns$position),shared_positions=nrow(both),
                 lollipop_only_positions=nrow(lonly),vaf_five_type_only_positions=nrow(vonly))
  recon_summary<-lr
  recon_summary[,`:=`(N_CDS=count[region=='CDS',N],N_rRNA=count[region=='rRNA',N],
    N_tRNA=count[region=='tRNA',N],N_Noncoding=count[region=='Noncoding',N])]
  recon[[paste0(sp,'_summary')]]<-recon_summary
  evl<-unique(lollipop[species==sp,.(species,mutation,position)])
  evn<-unique(event[,.(species,mutation,position)])
  evboth<-merge(evl,evn,by=c('species','mutation','position'))
  evlo<-evl[!evn,on=.(species,mutation,position)];evlo[,status:='lollipop_event_only']
  evno<-evn[!evl,on=.(species,mutation,position)];evno[,status:='vaf_five_celltype_event_only']
  if(nrow(evlo)||nrow(evno))event_recon[[sp]]<-rbindlist(list(evlo,evno),fill=TRUE)
}
summary<-rbindlist(all_sites);summary<-summary[order(match(species,c('human','mouse','rat','chicken','cobra','sturgeon','stingray')),match(region,regions))]
all_event_dt<-rbindlist(all_events,fill=TRUE)
all_cov<-rbindlist(coverage_counts,fill=TRUE)
recon_list<-recon[!grepl('_summary$',names(recon))]
recon_summary<-rbindlist(recon[grepl('_summary$',names(recon))],fill=TRUE)
recon_detail<-if(length(recon_list))rbindlist(recon_list,fill=TRUE)else data.table(species=character(),position=integer(),status=character())
event_recon_detail<-if(length(event_recon))rbindlist(event_recon,fill=TRUE)else data.table(species=character(),mutation=character(),position=integer(),status=character())
mask_dt<-if(length(mask_removed))rbindlist(mask_removed,fill=TRUE)else data.table(species=character(),position=integer(),mutation=character(),region=character(),legacy_low_coverage_mask=logical())
fwrite(summary,file.path(job,'Tables/regional_density_N_L.tsv'),sep='\t')
fwrite(all_event_dt,file.path(job,'Tables/vaf_five_celltype_mutation_sites.tsv'),sep='\t')
fwrite(all_cov,file.path(job,'Tables/five_celltype_coverage_barcode_counts.tsv'),sep='\t')
fwrite(recon_summary,file.path(job,'Tables/lollipop_position_reconciliation_summary.tsv'),sep='\t')
fwrite(recon_detail,file.path(job,'Tables/lollipop_position_reconciliation_differences.tsv'),sep='\t')
fwrite(event_recon_detail,file.path(job,'Tables/lollipop_variant_reconciliation_differences.tsv'),sep='\t')
fwrite(mask_dt,file.path(job,'Tables/vaf_records_excluded_by_original_mask.tsv'),sep='\t')

summary[,species:=as.character(species)]
summary[,region:=factor(region,levels=regions)]
summary[,panel_species:=factor(species,levels=c('human','mouse','rat','chicken','cobra','sturgeon','stingray'),labels=c('Human','Mouse','Rat','Chicken','Cobra','Sturgeon','Stingray'))]
summary[,coverage_flag:=as.character(panel_species)=='Sturgeon' & as.character(region)=='Noncoding' & L<100]
summary[,bar_label:=ifelse(coverage_flag,paste0('Coverage limited\nN=',N,'; L=',L),paste0('N=',N,'\nL=',L))]
cols<-c(CDS='#625D9E',rRNA='#B3659F',tRNA='#C7B1D4',Noncoding='#D3D3D3')
ymax<-max(summary$density_sites_per_kb,na.rm=TRUE)
p<-ggplot(summary,aes(x=region,y=density_sites_per_kb,fill=region))+
  geom_col(color='black',linewidth=.35,width=.76,na.rm=TRUE)+
  geom_text(aes(label=bar_label),vjust=-.2,size=2.15,lineheight=.9,na.rm=TRUE)+
  facet_wrap(~panel_species,ncol=2)+scale_fill_manual(values=cols,drop=FALSE)+
  scale_x_discrete(drop=FALSE)+scale_y_continuous(expand=expansion(mult=c(0,.22)))+
  coord_cartesian(ylim=c(0,ymax*1.22),clip='off')+
  labs(x=NULL,y='Regional mutation density (mutated positions/kb)',
       title='Regional mutation density across seven vertebrates',
       subtitle='N = unique positive-VAF mutation positions in five immune cell types; L = unmasked positions with at least one positive-depth coverage record',
       caption='Sturgeon Noncoding is flagged when L < 100 bp; its observed N and density are retained.')+
  theme_bw(base_size=9)+theme(panel.background=element_rect(fill='white',colour=NA),
    plot.background=element_rect(fill='white',colour=NA),panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),strip.background=element_rect(fill='white',colour='black'),
    strip.text=element_text(face='bold'),axis.text.x=element_text(angle=35,hjust=1),
    legend.position='none',plot.title=element_text(face='bold',hjust=.5),
    plot.subtitle=element_text(size=8),plot.caption=element_text(size=7,hjust=0))
ggsave(file.path(job,'Figures/regional_mutation_density_five_celltypes.pdf'),p,width=9.2,height=8.0,device=cairo_pdf)
ggsave(file.path(job,'Figures/regional_mutation_density_five_celltypes.png'),p,width=9.2,height=8.0,dpi=320,bg='white')
cat('regions=',nrow(summary),'events=',nrow(all_event_dt),'position differences=',nrow(recon_detail),
    'event differences=',nrow(event_recon_detail),'old-mask positive records=',nrow(mask_dt),'\n')
print(summary)
