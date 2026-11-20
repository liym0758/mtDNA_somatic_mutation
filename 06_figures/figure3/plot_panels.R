suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(grid)})
root<-Sys.getenv('MITOEVO_FIG2_RUN_ROOT')
stopifnot(nzchar(root), dir.exists(file.path(root,'source')))
if(dir.exists(file.path(root,'plots_spacing_revision')))stop('Refusing existing panel output')
dir.create(file.path(root,'logs'),showWarnings=FALSE,recursive=TRUE)
src<-file.path(root,'source'); out<-file.path(root,'plots_spacing_revision'); tab<-file.path(root,'tables')
dir.create(out,showWarnings=FALSE,recursive=TRUE);dir.create(tab,showWarnings=FALSE,recursive=TRUE)
logfile<-file.path(root,'logs','plot_panels.log')
sink(logfile,split=TRUE);on.exit(sink(),add=TRUE)

sps<-c('human','mouse','rat','chicken','cobra','sturgeon','stingray')
groups<-c('Monocyte','B cell','T cell','T/NK cell','NK cell')
main_groups<-c('Monocyte','Lymphoid')
metrics<-c('mt_percent','scmtMPM')
colors<-c('Monocyte'='#7FB0B4','B cell'='#B3659F','T cell'='#8E90B9',
          'NK cell'='#CCE0F5','T/NK cell'='#3361A5','Lymphoid'='#8E90B9')
spcolors<-setNames(c('#625D9E','#CCC9E6','#DCC1DD','#F2E2B2','#C5DEBA','#AAD4F8','#5184B2'),
                   tools::toTitleCase(sps))
theme_fig<-function()theme_bw(base_size=11)+theme(panel.grid=element_blank(),
 panel.border=element_rect(colour='black',linewidth=.45),strip.background=element_blank(),
 strip.text=element_text(face='bold',colour='black',size=10.5),axis.text=element_text(colour='black',size=9.5),
 axis.title=element_text(face='bold',size=10),plot.title=element_text(face='bold',size=12,hjust=.5),
 legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=8),
 panel.spacing.y=unit(.65,'lines'),plot.margin=margin(5,8,5,5))
star_label<-function(p)ifelse(p<.001,'***',ifelse(p<.01,'**',ifelse(p<.05,'*','')))
save_panel<-function(p,name,width,height){
 ggsave(file.path(out,paste0(name,'.pdf')),p,width=width,height=height,device=grDevices::cairo_pdf,bg='white')
 ggsave(file.path(out,paste0(name,'.png')),p,width=width,height=height,dpi=300,bg='white')
}
read_tsv<-function(name)read.delim(file.path(src,name),check.names=FALSE,stringsAsFactors=FALSE)

cells<-read_tsv('all_eligible_cells.tsv')
curves<-read_tsv('weighted_violin_density_curves.tsv')
wsummary<-read_tsv('weighted_group_summaries.tsv')
rank<-read_tsv('full_equal_weighted_rank.tsv')
species_means<-read_tsv('species_means_and_L_over_M_ratios.tsv')
stopifnot(nrow(cells)==52128L,all(c('mt_percent_3000reads','scmtMPM_read20','main_celltype2')%in%names(cells)),
          setequal(unique(cells$species),sps),nrow(rank)==14L)

# a/c: retain the original per-species normalization; median only the seven species.
bars<-species_means
stopifnot(nrow(bars)==14L,all(bars$species%in%sps),all(bars$metric%in%metrics))
median_summary<-bind_rows(lapply(metrics,function(metric){
 z<-bars[bars$metric==metric,]
 data.frame(metric=metric,n_species=length(unique(z$species)),
  median_relative_monocyte=median(z$relative_mean_monocyte),
  median_relative_lymphoid=median(z$relative_mean_lymphoid),
  median_difference_L_minus_M=median(z$relative_mean_lymphoid)-median(z$relative_mean_monocyte),
  direction=if(median(z$relative_mean_lymphoid)>median(z$relative_mean_monocyte))'Lymphoid higher'else'Lymphoid lower')
}))
stopifnot(nrow(median_summary)==2L,
 median_summary$median_relative_monocyte[median_summary$metric=='mt_percent']>0.82,
 median_summary$median_relative_monocyte[median_summary$metric=='mt_percent']<0.84,
 median_summary$median_relative_lymphoid[median_summary$metric=='mt_percent']>1.16,
 median_summary$median_relative_lymphoid[median_summary$metric=='mt_percent']<1.18,
 median_summary$median_relative_monocyte[median_summary$metric=='scmtMPM']>1.01,
 median_summary$median_relative_monocyte[median_summary$metric=='scmtMPM']<1.03,
 median_summary$median_relative_lymphoid[median_summary$metric=='scmtMPM']>0.97,
 median_summary$median_relative_lymphoid[median_summary$metric=='scmtMPM']<0.99)
write.table(median_summary,file.path(tab,'cross_species_seven_species_median_bars.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(bars[,c('species','metric','mean_monocyte','mean_lymphoid','relative_mean_monocyte','relative_mean_lymphoid','lymphoid_over_monocyte')],
 file.path(tab,'cross_species_species_dots.tsv'),sep='\t',row.names=FALSE,quote=FALSE)

summary_plot<-function(metric){
 z<-bars[bars$metric==metric,];z$species_title<-factor(tools::toTitleCase(z$species),levels=names(spcolors))
 dotdata<-bind_rows(data.frame(species_title=z$species_title,group='Monocyte',value=z$relative_mean_monocyte),
                    data.frame(species_title=z$species_title,group='B/T/NK',value=z$relative_mean_lymphoid))
 dotdata$group<-factor(dotdata$group,levels=c('Monocyte','B/T/NK'))
 b<-data.frame(group=factor(c('Monocyte','B/T/NK'),levels=c('Monocyte','B/T/NK')),
                value=c(median(z$relative_mean_monocyte),median(z$relative_mean_lymphoid)))
 b$label_y<-vapply(as.character(b$group),function(g)max(dotdata$value[dotdata$group==g]),numeric(1))+.07*max(dotdata$value)
 set.seed(20261030)
 ggplot()+geom_col(data=b,aes(group,value),width=.58,fill=c(colors['Monocyte'],colors['Lymphoid']),
    colour='black',linewidth=.35)+
  geom_point(data=dotdata,aes(group,value,fill=species_title),shape=21,colour='black',size=2.7,stroke=.4,
    position=position_jitter(width=.11,height=0))+
  geom_text(data=b,aes(group,label_y,label=sprintf('%.2f',value)),vjust=0,size=3.7)+
  scale_fill_manual(values=spcolors,drop=FALSE)+scale_x_discrete(limits=c('Monocyte','B/T/NK'))+
  scale_y_continuous(limits=c(0,NA),breaks=scales::breaks_extended(n=4),expand=expansion(mult=c(0,.19)))+
  guides(fill=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=2.6)))+theme_fig()+
  theme(axis.title.x=element_blank(),legend.position='bottom',legend.box='horizontal')+
  labs(title=if(metric=='mt_percent')'Across-species mt%'else'Across-species scmtMPM',
       y='Normalized mean')
}
save_panel(summary_plot('mt_percent'),'Fig2a_mt_summary',4.1,3.7)
save_panel(summary_plot('scmtMPM'),'Fig2c_scmtMPM_summary',4.1,3.7)

# b/d: frozen all-cell subtype-weighted density, weighted quartiles, and existing rank tests.
curves<-curves[curves$version=='weighted'&curves$display,]
wsummary<-wsummary[wsummary$scheme=='equal_subtype',]
rank$species<-factor(rank$species,levels=sps)
write.table(rank[,c('species','metric','theta','standard_error','p_two_sided','p_holm','repeats','tie_probability')],
 file.path(tab,'weighted_rank_tests_used.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(wsummary[,c('species','metric','plot_group','n','mean','q05','q25','median','q75','q95','weighted_zero_fraction')],
 file.path(tab,'weighted_violin_group_summaries.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
zero<-wsummary[wsummary$species%in%c('rat','sturgeon')&wsummary$metric=='scmtMPM',
 c('species','plot_group','n','weighted_zero_fraction','raw_zero_cells')]
write.table(zero,file.path(tab,'Rat_Sturgeon_scmtMPM_zero_fractions.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
geo<-list()
for(sp in sps)for(metric in metrics){
 z<-curves[curves$species==sp&curves$metric==metric,]
 stopifnot(all(c('Monocyte','Lymphoid')%in%z$plot_group))
 mx<-max(z$density);if(mx<=0)mx<-1
 for(g in main_groups){
  d<-z[z$plot_group==g,];x<-match(g,main_groups);half<-.29*d$density/mx
  geo[[length(geo)+1L]]<-data.frame(species=sp,metric=metric,plot_group=g,
    px=c(x-half,rev(x+half)),py=c(d$y,rev(d$y)))
 }
}
geometry<-bind_rows(geo);geometry$species<-factor(geometry$species,levels=sps)
geometry$plot_group<-factor(geometry$plot_group,levels=main_groups)
wsummary$species<-factor(wsummary$species,levels=sps)
wsummary$plot_group<-factor(wsummary$plot_group,levels=main_groups)
wsummary$x<-match(as.character(wsummary$plot_group),main_groups)
main_brackets<-function(metric,ss){
 a<-rank[rank$metric==metric,];a<-as.data.frame(a);a$species<-factor(a$species,levels=sps)
 a<-a[a$p_holm<.05,]
 if(!nrow(a))return(a)
 a$ymax<-vapply(as.character(a$species),function(sp)max(ss$q95[ss$species==sp]),numeric(1))
 a$span<-pmax(a$ymax,.5);a$y<-a$ymax+.14*a$span;a$bottom<-a$y-.025*a$span;a$label<-star_label(a$p_holm)
 a$text_y<-a$y+.045*a$span;a$x2<-2
 a
}
main_violin<-function(metric){
 z<-geometry[geometry$metric==metric,];ss<-wsummary[wsummary$metric==metric,]
 lab<-setNames(tools::toTitleCase(sps),sps);a<-main_brackets(metric,ss)
 p<-ggplot()+geom_polygon(data=z,aes(px,py,fill=plot_group,group=interaction(species,plot_group)),
    colour='grey40',alpha=.9,linewidth=.4)+
  geom_rect(data=ss,aes(xmin=x-.045,xmax=x+.045,ymin=q25,ymax=q75),fill='white',colour='grey25',linewidth=.35)+
  geom_segment(data=ss,aes(x=x-.045,xend=x+.045,y=median,yend=median),colour='grey25',linewidth=.4)+
  facet_wrap(~species,ncol=1,scales='free_y',strip.position='top',labeller=as_labeller(lab))+
  scale_fill_manual(values=colors)+scale_x_continuous(breaks=1:2,labels=c('Monocyte','B/T/NK'),limits=c(.4,2.6))+
  scale_y_continuous(limits=c(0,NA),breaks=scales::breaks_extended(n=3),expand=expansion(mult=c(0,.13)))+
  theme_fig()+theme(legend.position='none',axis.title.x=element_blank(),axis.title.y=element_text(size=10),
    axis.text.x=element_text(size=9),strip.placement='outside',strip.text.x=element_text(size=9.5,face='bold'),
    panel.spacing.y=unit(.48,'lines'))+
  labs(title=NULL,y=if(metric=='mt_percent')'mt%'else'scmtMPM',x=NULL)
 if(nrow(a))p<-p+geom_segment(data=a,aes(x=1,xend=2,y=y,yend=y),inherit.aes=FALSE,linewidth=.35)+
   geom_segment(data=a,aes(x=1,xend=1,y=bottom,yend=y),inherit.aes=FALSE,linewidth=.35)+
   geom_segment(data=a,aes(x=2,xend=2,y=bottom,yend=y),inherit.aes=FALSE,linewidth=.35)+
   geom_text(data=a,aes(x=1.5,y=text_y,label=label),inherit.aes=FALSE,size=4,fontface='bold')
 p
}
save_panel(main_violin('mt_percent'),'Fig2b_mt_violin',2.7,10.5)
save_panel(main_violin('scmtMPM'),'Fig2d_scmtMPM_violin',2.7,10.5)

# Supplementary panels: all eligible cells, old class groups, ordinary Wilcoxon and Holm.
class_groups<-function(sp){if(sp%in%c('human','mouse','rat','chicken'))c('B cell','T cell','NK cell')else c('B cell','T/NK cell')}

# Supplementary cell-type summaries: normalize each species' subtype means to
# the mean of its available cell-type means, then summarize species equally.
celltype_levels<-c('Monocyte','B cell','T cell','NK cell','T/NK cell')
supp_summary_rows<-list();ss_i<-0L
for(metric in metrics)for(sp in sps){
 mcol<-if(metric=='mt_percent')'mt_percent_3000reads'else'scmtMPM_read20'
 available<-c('Monocyte',class_groups(sp))
 species_means<-vapply(available,function(g){
   ii<-which(cells$species==sp&cells$main_celltype2==g)
   if(!length(ii))stop(paste('Missing summary group',sp,g,metric))
   mean(cells[[mcol]][ii])
 },numeric(1))
 baseline<-mean(species_means)
 ss_i<-ss_i+1L
 supp_summary_rows[[ss_i]]<-data.frame(species=sp,metric=metric,plot_group=available,
    n_cells=vapply(available,function(g)sum(cells$species==sp&cells$main_celltype2==g),integer(1)),
    mean_value=unname(species_means),species_baseline=baseline,
    normalized_mean=unname(species_means)/baseline)
}
supp_summary_species<-bind_rows(supp_summary_rows)
supp_summary_species$plot_group<-factor(supp_summary_species$plot_group,levels=celltype_levels)
supp_summary_species$species_title<-factor(tools::toTitleCase(supp_summary_species$species),levels=names(spcolors))
supp_summary_medians<-supp_summary_species%>%group_by(metric,plot_group)%>%
 summarize(n_species=n_distinct(species),median_normalized_mean=median(normalized_mean),.groups='drop')
species_normalization_check<-supp_summary_species%>%group_by(species,metric)%>%
 summarize(mean_normalized_mean=mean(normalized_mean),.groups='drop')
stopifnot(all(abs(species_normalization_check$mean_normalized_mean-1)<1e-12))
expected_species_n<-c('Monocyte'=7L,'B cell'=7L,'T cell'=4L,'NK cell'=4L,'T/NK cell'=3L)
for(metric in metrics){
 z<-supp_summary_medians[supp_summary_medians$metric==metric,]
 stopifnot(identical(as.integer(z$n_species),
                     as.integer(expected_species_n[as.character(z$plot_group)])))
}
write.table(supp_summary_species,file.path(tab,'supplementary_celltype_summary_species_means.tsv'),
 sep='\t',row.names=FALSE,quote=FALSE)
write.table(supp_summary_medians,file.path(tab,'supplementary_celltype_summary_medians.tsv'),
 sep='\t',row.names=FALSE,quote=FALSE)
supp_summary_plot<-function(metric){
 d<-supp_summary_species[supp_summary_species$metric==metric,]
 b<-supp_summary_medians[supp_summary_medians$metric==metric,]
 b$label_y<-vapply(as.character(b$plot_group),function(g)
   max(d$normalized_mean[d$plot_group==g]),numeric(1))+.055*max(d$normalized_mean)
 b$plot_group<-factor(b$plot_group,levels=celltype_levels)
 ggplot()+geom_col(data=b,aes(plot_group,median_normalized_mean,fill=plot_group),
    width=.62,colour='black',linewidth=.35)+
  geom_point(data=d,aes(plot_group,normalized_mean,colour=species_title),
    size=2.5,position=position_jitter(width=.09,height=0))+
  geom_text(data=b,aes(plot_group,label_y,label=sprintf('%.2f',median_normalized_mean)),
    vjust=0,size=3.4)+
  scale_fill_manual(values=colors,drop=TRUE)+scale_colour_manual(values=spcolors,drop=FALSE)+
  scale_x_discrete(limits=celltype_levels,drop=TRUE)+
  scale_y_continuous(limits=c(0,NA),breaks=scales::breaks_extended(n=4),
    expand=expansion(mult=c(0,.18)))+theme_fig()+
  guides(fill='none',colour=guide_legend(nrow=2,byrow=TRUE,override.aes=list(size=2.5)))+
  theme(axis.title.x=element_blank(),legend.position='bottom',legend.box='horizontal')+
  labs(title=if(metric=='mt_percent')'Across-species mt% by cell type'else'Across-species scmtMPM by cell type',
       y='Normalized mean')
}
save_panel(supp_summary_plot('mt_percent'),'FigS_mt_celltype_summary',6.2,4.1)
save_panel(supp_summary_plot('scmtMPM'),'FigS_scmtMPM_celltype_summary',6.2,4.1)

testrows<-list();displayrows<-list();auditrows<-list();k<-kd<-ka<-0L
for(metric in metrics)for(sp in sps){
 mcol<-if(metric=='mt_percent')'mt_percent_3000reads'else'scmtMPM_read20'
 imono<-which(cells$species==sp&cells$main_celltype2=='Monocyte')
 mono<-cells[[mcol]][imono]
 for(g in class_groups(sp)){
  ig<-which(cells$species==sp&cells$main_celltype2==g)
  val<-cells[[mcol]][ig]
  if(!length(mono)||!length(val))stop(paste('Missing group',sp,g,metric))
  p<-if(length(unique(c(mono,val)))==1)1 else suppressWarnings(wilcox.test(val,mono,alternative='two.sided',exact=FALSE,correct=TRUE)$p.value)
  k<-k+1L;testrows[[k]]<-data.frame(species=sp,metric=metric,other_group=g,n_monocyte=length(mono),n_other=length(val),
    mean_monocyte=mean(mono),mean_other=mean(val),median_monocyte=median(mono),median_other=median(val),
    wilcoxon_p_two_sided=p,p_holm=NA_real_,holm_family='all 18 species-by-subtype comparisons within metric',
    method='two-sided ordinary Wilcoxon rank-sum; tie correction; continuity correction; normal approximation')
  q<-unname(quantile(val,c(.05,.95),type=7));keep<-val>=q[1]&val<=q[2]
  a<-cells[ig[keep],c('species','barcode','main_celltype2')]
  a$plot_group<-a$main_celltype2;a$metric<-metric;a$value<-val[keep]
  displayrows[[length(displayrows)+1L]]<-a
  auditrows[[length(auditrows)+1L]]<-data.frame(species=sp,metric=metric,plot_group=g,n_test=length(val),n_display=sum(keep),
    p05=q[1],p95=q[2],statistics_extremes_removed=0)
  if(g==class_groups(sp)[1]){
    vm<-cells[[mcol]][imono];qm<-quantile(vm,c(.05,.95),type=7);keepm<-vm>=qm[1]&vm<=qm[2]
    am<-cells[imono[keepm],c('species','barcode','main_celltype2')];am$plot_group<-'Monocyte';am$metric<-metric;am$value<-vm[keepm]
    displayrows[[length(displayrows)+1L]]<-am
    auditrows[[length(auditrows)+1L]]<-data.frame(species=sp,metric=metric,plot_group='Monocyte',n_test=length(vm),n_display=sum(vm>=qm[1]&vm<=qm[2]),
      p05=qm[1],p95=qm[2],statistics_extremes_removed=0)
  }
 }
}
class_tests<-bind_rows(testrows)
class_tests$p_holm<-ave(class_tests$wilcoxon_p_two_sided,class_tests$metric,FUN=function(p)p.adjust(p,'holm'))
stopifnot(all(table(class_tests$metric)==18L),all(!is.na(class_tests$p_holm)))
class_display<-bind_rows(displayrows)
class_audit<-bind_rows(auditrows)
write.table(class_tests,file.path(tab,'split_celltype_wilcoxon_all_eligible.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(class_audit,file.path(tab,'split_celltype_P5_P95_display_audit.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(class_display,file.path(tab,'split_celltype_P5_P95_display_data.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
class_tests$species<-factor(class_tests$species,levels=sps)
class_display$species<-factor(class_display$species,levels=sps)
class_display$plot_group<-factor(class_display$plot_group,levels=groups)

class_brackets<-function(metric,display){
  a<-class_tests[class_tests$metric==metric,]
 if(!nrow(a))return(a)
 a$species<-factor(a$species,levels=sps)
 a$x2<-vapply(seq_len(nrow(a)),function(i){
   shown<-groups[groups%in%as.character(display$plot_group[display$species==a$species[i]])]
   match(a$other_group[i],shown)
 },numeric(1))
 a$ymax<-vapply(as.character(a$species),function(sp)max(display$value[display$species==sp]),numeric(1))
 a$span<-pmax(a$ymax,.5)
 a$stack<-a$x2-2
 a$y<-a$ymax+a$span*(.12+.11*a$stack);a$bottom<-a$y-.02*a$span
 a$text_y<-a$y+.035*a$span
 a$label<-ifelse(a$p_holm<.05,star_label(a$p_holm),'ns')
 a
}
class_violin<-function(metric){
 z<-class_display[class_display$metric==metric,];a<-class_brackets(metric,z)
 p<-ggplot(z,aes(plot_group,value,fill=plot_group))+geom_violin(width=.68,trim=TRUE,colour='grey40',alpha=.9,linewidth=.35)+
   geom_boxplot(width=.085,fill='white',colour='grey25',outlier.shape=NA,linewidth=.3)+
   facet_wrap(~species,ncol=1,scales='free',strip.position='top',labeller=as_labeller(setNames(tools::toTitleCase(sps),sps)))+
   scale_fill_manual(values=colors,drop=TRUE)+scale_x_discrete(drop=TRUE)+
   scale_y_continuous(limits=c(0,NA),breaks=scales::breaks_extended(n=3),expand=expansion(mult=c(0,.14)))+
   theme_fig()+theme(legend.position='none',axis.title.x=element_blank(),axis.title.y=element_text(size=10),
     axis.text.x=element_text(angle=30,hjust=1,size=8.5),strip.placement='outside',
     strip.text.x=element_text(size=9.5,face='bold'),panel.spacing.y=unit(.35,'lines'))+
   labs(title=NULL,y=if(metric=='mt_percent')'mt%'else'scmtMPM',x=NULL)
 if(nrow(a))p<-p+geom_segment(data=a,aes(x=1,xend=x2,y=y,yend=y),inherit.aes=FALSE,linewidth=.3)+
   geom_segment(data=a,aes(x=1,xend=1,y=bottom,yend=y),inherit.aes=FALSE,linewidth=.3)+
   geom_segment(data=a,aes(x=x2,xend=x2,y=bottom,yend=y),inherit.aes=FALSE,linewidth=.3)+
   geom_text(data=a,aes(x=(1+x2)/2,y=text_y,label=label),inherit.aes=FALSE,size=3.6,fontface='bold')
 p
}
save_panel(class_violin('mt_percent'),'FigS_mt_split_violin',3.4,11.0)
save_panel(class_violin('scmtMPM'),'FigS_scmtMPM_split_violin',3.4,13.5)

cat('PANELS_COMPLETE\n')
cat('eligible_cells=',nrow(cells),'\n',sep='')
cat('supplementary_wilcoxon_tests_per_metric=',nrow(class_tests)/2,'; Holm within metric\n',sep='')
cat('median_summary:\n');print(median_summary)
cat('main_zero_fractions_rat_sturgeon:\n');print(zero)
