#!/usr/bin/env Rscript
# Figure 4 CI and binned VAF-trend extension. Reads validated all_positive results;
# never modifies the source filtering, site counts, mean VAF, or PAML background.
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})

arg <- grep('^--file=', commandArgs(), value=TRUE)
if(length(arg)!=1L) stop('Run with Rscript plot_fig4_ci_trend.R')
out <- dirname(normalizePath(sub('^--file=','',arg), winslash='/', mustWork=TRUE))
code <- Sys.getenv('MITOEVO_SOURCE_CODE_ROOT')
if(!nzchar(code) || !dir.exists(code)) stop('Set MITOEVO_SOURCE_CODE_ROOT to the local source code directory')
src <- function(...) file.path(code,...)
path <- list(
  pooled=src('mitoevo_targeted_analysis','Figure4_pooled_bin_comparisons_20261105_v1','pooled_lineages_four_vaf_bin_counts.tsv'),
  sites=src('mitoevo_targeted_analysis','Tables','coding_site_all_unique.tsv'),
  grouped=src('manuscript_story_task1','Tables','figure4_56_original_results_with_display_labels.tsv'),
  grouped_overall=src('manuscript_story_task1','Tables','figure4_14_original_overall_results.tsv'),
  prior=src('manuscript_revision_20261104_v1','tables','fig4_dnds_CI.tsv'),
  background=src('manuscript_revision_20261104_v1','remote_outputs','tables','background_sites_summary.csv'),
  lineage_sites=src('manuscript_revision_20261104_v1','remote_outputs','tables','fig4_site_vaf_both_filters.tsv'))
for(p in path) if(!file.exists(p)) stop('Missing input: ',p)
sp_order <- c('human','mouse','rat','chicken','cobra','sturgeon','stingray')
sp_title <- c('Human','Mouse','Rat','Chicken','Cobra','Sturgeon','Stingray')
bin_order <- c('(0,0.1]','(0.1,0.3]','(0.3,0.5]','(0.5,1]')
midpoints <- c(0.05,0.20,0.40,0.75)

pooled_raw <- read.delim(path$pooled,check.names=FALSE)
site <- read.delim(path$sites,check.names=FALSE)
grouped_raw <- read.delim(path$grouped,check.names=FALSE)
grouped_overall_raw <- read.delim(path$grouped_overall,check.names=FALSE)
prior <- read.delim(path$prior,check.names=FALSE)
bg <- read.csv(path$background,check.names=FALSE)
lineage_sites <- read.delim(path$lineage_sites,check.names=FALSE)
bg$expected_NS <- bg$non_synonymous_sites/bg$synonymous_sites

# Provenance checks on all_positive carrier observations and mutation identity.
stopifnot(nrow(site)==977L, all(site$lineage=='All_unique'),
          all(site$mean_vaf>0 & site$mean_vaf<=1),
          all(site$n_observations>0),
          all(site$mean_vaf==site$mean_vaf | is.na(site$mean_vaf)==FALSE),
          !anyDuplicated(site[c('species','mutation')]),
          all(site$bin==bin_order[findInterval(site$mean_vaf,c(0,.1,.3,.5,1),rightmost.closed=TRUE,left.open=TRUE)]))
ls <- subset(lineage_sites,filter_mode=='legacy_all_positive' &
               mutation_type %in% c('nonsynonymous','synonymous'))
stopifnot(!anyDuplicated(ls[c('species','mutation','lineage')]),
          all(ls$mean_vaf>0 & ls$mean_vaf<=1),
          all(ls$n_observations>0))
key <- paste(site$species,site$mutation,sep='|')
ls <- ls[paste(ls$species,ls$mutation,sep='|') %in% key,]
agg_vaf <- aggregate(cbind(vaf_sum=ls$mean_vaf*ls$n_observations,
                           n_observations=ls$n_observations),
                     by=list(species=ls$species,mutation=ls$mutation),FUN=sum)
joined <- merge(site,agg_vaf,by=c('species','mutation'),suffixes=c('_saved','_from_lineages'))
stopifnot(nrow(joined)==nrow(site),
          all(joined$n_observations_saved==joined$n_observations_from_lineages),
          all(abs(joined$mean_vaf-joined$vaf_sum_from_lineages/joined$n_observations_from_lineages)<1e-10))

# Clopper-Pearson interval for p=N/(N+S), then odds(p)/PAML N:S opportunity ratio.
cp_dnds <- function(N,S,bg_ratio) {
  n <- N+S
  if(n==0) return(c(ci_low=NA_real_,ci_high=NA_real_))
  ci <- binom.test(N,n,conf.level=.95)$conf.int
  odds <- function(p) ifelse(p>=1,Inf,p/(1-p))
  c(ci_low=odds(ci[1])/bg_ratio,ci_high=odds(ci[2])/bg_ratio)
}
add_ci <- function(d) {
  ci <- t(vapply(seq_len(nrow(d)),function(i) cp_dnds(d$N[i],d$S[i],d$expected_NS[i]),
                 numeric(2)))
  d$ci_low <- ci[,1]; d$ci_high <- ci[,2]
  d$status <- ifelse(d$n==0,'empty',ifelse(d$S==0,'undefined_S_zero','defined'))
  d$midpoint <- ifelse(d$bin=='All',NA_real_,midpoints[match(d$bin,bin_order)])
  d$recomputed_dnds <- ifelse(d$S>0,d$N/d$S/d$expected_NS,NA_real_)
  ok <- d$status=='defined'
  stopifnot(all(abs(d$recomputed_dnds[ok]-d$dnds[ok])<1e-9),
            all(is.na(d$dnds[!ok]) | is.infinite(d$dnds[!ok])))
  d$dnds[d$status!='defined'] <- NA_real_
  d
}

# A: one unique position-REF-ALT site per species, across target immune cells.
stopifnot(nrow(pooled_raw)==28L,all(pooled_raw$bin==bin_order[pooled_raw$bin_number]))
pooled <- data.frame(species=pooled_raw$species,lineage='All_unique',bin=pooled_raw$bin,
                     N=pooled_raw$N_obs,S=pooled_raw$S_obs,
                     n=pooled_raw$n_unique_coding_sites,
                     expected_NS=pooled_raw$PAML_N_over_S_background,
                     dnds=pooled_raw$dnds_original_background,
                     n_positive_observations=pooled_raw$n_positive_cell_variant_observations)
pooled_prior <- subset(prior,filter_mode=='legacy_all_positive' & lineage=='All_unique' & bin!='All')
stopifnot(nrow(pooled_prior)==28L)
pooled <- merge(pooled,pooled_prior[c('species','bin','N','S','p_neutral')],
                by=c('species','bin'),suffixes=c('','_prior'),sort=FALSE)
stopifnot(all(pooled$N==pooled$N_prior),all(pooled$S==pooled$S_prior))
pooled$N_prior <- NULL; pooled$S_prior <- NULL
pooled <- add_ci(pooled)

pooled_overall_prior <- subset(prior,filter_mode=='legacy_all_positive' & lineage=='All_unique' & bin=='All')
stopifnot(nrow(pooled_overall_prior)==7L)
pooled_overall <- data.frame(species=pooled_overall_prior$species,lineage='All_unique',bin='All',
                             N=pooled_overall_prior$N,S=pooled_overall_prior$S,
                             n=pooled_overall_prior$n,expected_NS=pooled_overall_prior$expected_NS_ratio,
                             dnds=pooled_overall_prior$dnds,p_neutral=pooled_overall_prior$p_neutral)
pooled_overall <- add_ci(pooled_overall)

# B: original 56 lineage-bin rows and 14 original lineage-overall rows.
stopifnot(nrow(grouped_raw)==56L,!anyDuplicated(grouped_raw[c('species','celltype_group','bin')]))
grouped <- data.frame(species=grouped_raw$species,lineage=grouped_raw$celltype_group,
                      bin=grouped_raw$bin,N=grouped_raw$nonsynonymous,
                      S=grouped_raw$synonymous,n=grouped_raw$total_mut,
                      expected_NS=grouped_raw$NS_ratio_expected,dnds=grouped_raw$dnds,
                      p_neutral=grouped_raw$p_binom)
grouped <- add_ci(grouped)
ls$bin_check <- as.character(cut(ls$mean_vaf,breaks=c(0,.1,.3,.5,1),labels=bin_order))
for(i in seq_len(nrow(grouped))) {
  z <- ls[ls$species==grouped$species[i] & ls$lineage==grouped$lineage[i] &
            ls$bin_check==grouped$bin[i],]
  stopifnot(sum(z$mutation_type=='nonsynonymous')==grouped$N[i],
            sum(z$mutation_type=='synonymous')==grouped$S[i])
}
grouped_overall <- data.frame(species=tolower(grouped_overall_raw$species),
                              lineage=grouped_overall_raw$celltype_group,bin='All',
                              N=grouped_overall_raw$nonsynonymous,
                              S=grouped_overall_raw$synonymous,
                              n=grouped_overall_raw$total_mut,
                              dnds=grouped_overall_raw$dnds,
                              p_neutral=grouped_overall_raw$p_vs_bg)
grouped_overall$expected_NS <- bg$expected_NS[match(grouped_overall$species,bg$species)]
stopifnot(nrow(grouped_overall)==14L,all(!is.na(grouped_overall$expected_NS)))
grouped_overall <- add_ci(grouped_overall)

# Cross-check all tabulated counts and the unchanged PAML opportunities.
for(sp in sp_order) {
  b <- subset(pooled,species==sp); z <- subset(site,species==sp)
  o <- subset(pooled_overall,species==sp)
  stopifnot(sum(b$N)==sum(z$is_nonsyn),sum(b$N)==o$N,
            sum(b$S)==sum(1-z$is_nonsyn),sum(b$S)==o$S,
            all(abs(b$expected_NS-bg$expected_NS[match(sp,bg$species)])<1e-10))
  for(lin in c('T_B_NK','Mono')) {
    g <- subset(grouped,species==sp & lineage==lin)
    go <- subset(grouped_overall,species==sp & lineage==lin)
    stopifnot(sum(g$N)==go$N,sum(g$S)==go$S,
              all(abs(g$expected_NS-bg$expected_NS[match(sp,bg$species)])<1e-10))
  }
}

# One binned logistic trend per species and analysis group. Empty bins are excluded.
fit_trend <- function(d) {
  d <- d[d$n>0,]
  ans <- data.frame(slope=NA_real_,se=NA_real_,odds_ratio_per_0.1_VAF=NA_real_,
                    p_trend=NA_real_,n_effective_bins=nrow(d),n_sites=sum(d$n),
                    N_total=sum(d$N),S_total=sum(d$S),status='')
  if(nrow(d)<2) {ans$status <- 'fewer_than_2_nonempty_bins'; return(ans)}
  if(sum(d$n)<10) {ans$status <- 'fewer_than_10_sites_unreliable'; return(ans)}
  if(sum(d$N)==0 || sum(d$S)==0) {ans$status <- 'single_mutation_class'; return(ans)}
  warnings <- character()
  fit <- tryCatch(withCallingHandlers(
    glm(cbind(N,S)~midpoint,family=binomial(),data=d),
    warning=function(w) {warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')}),
    error=function(e) NULL)
  if(is.null(fit)) {ans$status <- 'glm_error'; return(ans)}
  cf <- summary(fit)$coefficients
  if(nrow(cf)<2 || !isTRUE(fit$converged) ||
     any(!is.finite(cf[2,c('Estimate','Std. Error','Pr(>|z|)')])) ||
     any(grepl('separation|fitted probabilities|did not converge',warnings,ignore.case=TRUE)) ||
     abs(cf[2,'Estimate'])>50 || cf[2,'Std. Error']>20 ||
     any(fitted(fit)<1e-7 | fitted(fit)>1-1e-7)) {
    ans$status <- paste0('unstable_or_separated',if(length(warnings)) paste0(': ',paste(unique(warnings),collapse=' | ')) else '')
    return(ans)
  }
  ans$slope <- cf[2,'Estimate']; ans$se <- cf[2,'Std. Error']
  ans$odds_ratio_per_0.1_VAF <- exp(0.1*ans$slope)
  ans$p_trend <- cf[2,'Pr(>|z|)']; ans$status <- 'estimated'
  ans
}
trend_rows <- list(); k <- 0L
for(analysis in c('All_unique','T_B_NK','Mono')) for(sp in sp_order) {
  d <- if(analysis=='All_unique') subset(pooled,species==sp) else subset(grouped,species==sp & lineage==analysis)
  k <- k+1L
  trend_rows[[k]] <- cbind(data.frame(analysis=analysis,species=sp),fit_trend(d))
}
trends <- do.call(rbind,trend_rows)
trends$q_BH_across_7_species <- NA_real_
for(analysis in unique(trends$analysis)) {
  ix <- which(trends$analysis==analysis & !is.na(trends$p_trend))
  if(length(ix)) trends$q_BH_across_7_species[ix] <- p.adjust(trends$p_trend[ix],method='BH',n=7)
}
trends$BH_planned_family_size <- 7L

write.csv(pooled,file.path(out,'Figure4A_pooled_bin_CI.csv'),row.names=FALSE,na='')
write.csv(pooled_overall,file.path(out,'Figure4A_pooled_overall_CI.csv'),row.names=FALSE,na='')
write.csv(grouped,file.path(out,'Figure4B_lineage_bin_CI.csv'),row.names=FALSE,na='')
write.csv(grouped_overall,file.path(out,'Figure4B_lineage_overall_CI.csv'),row.names=FALSE,na='')
write.csv(trends,file.path(out,'Figure4_VAF_bin_trends.csv'),row.names=FALSE,na='')

purple <- '#8E90B9'; teal <- '#7FB0B4'
theme_fig <- theme_bw(base_size=9.5)+theme(panel.grid=element_blank(),
  plot.title=element_text(face='bold',size=10.5,hjust=0),
  plot.subtitle=element_text(size=8.5),axis.text=element_text(color='black'),
  plot.margin=margin(6,8,7,6))
fmt_p <- function(p) if(is.na(p)) 'NA' else if(p<.001) '<0.001' else sprintf('%.3f',p)
get_trend <- function(sp,lin) trends[trends$species==sp & trends$analysis==lin,][1,]

plot_bin <- function(d,sp,lin,title,fill) {
  d <- d[match(bin_order,d$bin),]
  d$x <- seq_len(4); d$count_label <- sprintf('N=%d, S=%d',d$N,d$S)
  d$status_label <- ifelse(d$status=='empty','No coding\nvariants',
                           ifelse(d$status=='undefined_S_zero','NA (S=0)',''))
  defined <- subset(d,status=='defined')
  defined$neutral_label <- ifelse(defined$p_neutral<.001,'***',
                                  ifelse(defined$p_neutral<.01,'**',
                                         ifelse(defined$p_neutral<.05,'*','ns')))
  empty <- subset(d,status=='empty')
  undefined <- subset(d,status=='undefined_S_zero')
  ci <- subset(d,n>0)
  cap <- 2.55
  ci$plot_low <- pmin(ci$ci_low,cap-.14)
  ci$plot_high <- pmin(ci$ci_high,cap)
  normal <- ci[is.finite(ci$ci_high) & ci$ci_high<=cap & ci$ci_low<=cap,]
  clipped <- ci[!is.finite(ci$ci_high) | ci$ci_high>cap | ci$ci_low>cap,]
  clipped$x_ci <- ifelse(clipped$status=='undefined_S_zero',clipped$x-.12,clipped$x)
  tr <- get_trend(sp,lin)
  subtitle <- sprintf('P-trend = %s (%d bins)',fmt_p(tr$p_trend),tr$n_effective_bins)
  ggplot(d,aes(x=x))+
    geom_hline(yintercept=1,linetype='dashed',color='grey40',linewidth=.4)+
    geom_errorbar(data=normal,aes(ymin=ci_low,ymax=ci_high),width=.13,linewidth=.48,color='grey25')+
    geom_segment(data=clipped,aes(x=x_ci,xend=x_ci,y=plot_low,yend=plot_high),
                 arrow=grid::arrow(length=grid::unit(.075,'in'),type='closed'),
                 linewidth=.48,color='grey25')+
    geom_point(data=defined,aes(y=dnds,size=pmin(n,100)),shape=21,fill=fill,color='black',stroke=.7)+
    geom_text(data=defined,aes(x=x+.15,y=pmin(dnds+.17,2.62),label=neutral_label),
              size=2.8,color='black')+
    scale_size_continuous(range=c(2.2,5),limits=c(0,100),guide='none')+
    geom_text(aes(y=2.81,label=count_label),size=2.5)+
    geom_text(data=empty,aes(y=.32,label=status_label),size=2.25,lineheight=.9,color='grey25')+
    geom_label(data=undefined,aes(y=.32,label=status_label),size=2.25,
               label.size=0,fill='white',color='grey25')+
    scale_x_continuous(breaks=1:4,labels=bin_order,limits=c(.63,4.37))+
    scale_y_continuous(limits=c(0,3),breaks=c(0,1,2,3))+
    labs(title=title,subtitle=subtitle,x='Mean VAF bin (per mutation)',y='dN/dS')+
    theme_fig+theme(axis.text.x=element_text(angle=42,hjust=1,size=7.2))
}

star <- function(p) ifelse(is.na(p),'',ifelse(p<.001,'***',ifelse(p<.01,'**',ifelse(p<.05,'*','ns'))))
plot_overall <- function(d,grouped=FALSE) {
  d$spx <- match(d$species,sp_order)
  d$x <- d$spx + if(grouped) ifelse(d$lineage=='T_B_NK',-.18,.18) else 0
  d$fill <- if(grouped) ifelse(d$lineage=='T_B_NK',purple,teal) else purple
  d$count_label <- sprintf('N=%d\nS=%d',d$N,d$S)
  d$star <- star(d$p_neutral)
  cap <- 1.8
  normal <- d[is.finite(d$ci_high) & d$ci_high<=cap,]
  clipped <- d[!is.finite(d$ci_high) | d$ci_high>cap,]
  clipped$plot_low <- pmin(clipped$ci_low,cap-.12)
  p <- ggplot(d,aes(x=x,y=dnds))+
    geom_hline(yintercept=1,linetype='dashed',color='grey40',linewidth=.4)+
    geom_col(aes(fill=lineage),width=if(grouped) .32 else .65,color='black',linewidth=.35)+
    geom_errorbar(data=normal,aes(ymin=ci_low,ymax=ci_high),width=if(grouped) .08 else .13,
                  linewidth=.42,color='grey20')+
    geom_segment(data=clipped,aes(x=x,xend=x,y=plot_low,yend=cap),
                 arrow=grid::arrow(length=grid::unit(.065,'in'),type='closed'),
                 linewidth=.42,color='grey20')+
    geom_text(aes(y=1.97,label=star),size=if(grouped) 2.55 else 3.0)+
    geom_text(aes(y=2.19,label=count_label),size=if(grouped) 2.02 else 2.15,lineheight=.87)+
    scale_x_continuous(breaks=1:7,labels=sp_title,limits=c(.58,7.42))+
    scale_y_continuous(limits=c(0,2.33),breaks=c(0,.5,1,1.5,2))+
    scale_fill_manual(values=c('All_unique'=purple,'T_B_NK'=purple,'Mono'=teal),
                      breaks=if(grouped) c('T_B_NK','Mono') else 'All_unique',
                      labels=if(grouped) c('T_B_NK','Mono') else 'Pooled')+
    labs(title='h  Overall mtDNA dN/dS',x='Species',y='dN/dS',fill=NULL)+theme_fig+
    theme(axis.text.x=element_text(angle=45,hjust=1,size=7.1),
          legend.position=if(grouped) 'bottom' else 'none',
          legend.text=element_text(size=8),legend.key.size=grid::unit(.27,'in'))
  p
}

foot <- ggdraw()+
  draw_label('95% CI: Clopper-Pearson N/(N+S), transformed using the original PAML N:S opportunity ratio.',
             x=0,y=.70,hjust=0,size=7.6)+
  draw_label('Upward arrows indicate a CI above the plotted range. Stars/ns test dN/dS versus 1; P-trend is a separate binned logistic test.',
             x=0,y=.26,hjust=0,size=7.6)

pa <- lapply(seq_along(sp_order),function(i) {
  sp <- sp_order[i]
  plot_bin(subset(pooled,species==sp),sp,'All_unique',paste0(letters[i],'  ',sp_title[i]),purple)
})
pa[[8]] <- plot_overall(pooled_overall)
fig_a <- plot_grid(plot_grid(plotlist=pa,ncol=2,align='hv'),foot,ncol=1,rel_heights=c(1,.045))
ggsave(file.path(out,'Figure4A_pooled_CI_trend.pdf'),fig_a,width=9,height=12.7,device=cairo_pdf,bg='white')
ggsave(file.path(out,'Figure4A_pooled_CI_trend.png'),fig_a,width=9,height=12.7,dpi=220,bg='white')

# Place individual plots directly on the page. Nested plot_grid letterboxes
# 8:3 source-style panels, leaving excessive blank space between species.
fig_b <- ggdraw()
for(i in seq_along(sp_order)) {
  sp <- sp_order[i]
  row <- (i-1L)%/%2L
  left <- if(i%%2L==1L) 0 else .5
  bottom <- .06+(3-row)*.235
  p1 <- plot_bin(subset(grouped,species==sp & lineage=='T_B_NK'),sp,'T_B_NK',
                 paste0(letters[i],'  ',sp_title[i],' - T_B_NK'),purple)
  p2 <- plot_bin(subset(grouped,species==sp & lineage=='Mono'),sp,'Mono',
                 paste0(sp_title[i],' - Mono'),teal)
  fig_b <- fig_b + draw_plot(p1,x=left,y=bottom,width=.25,height=.235) +
    draw_plot(p2,x=left+.25,y=bottom,width=.25,height=.235)
}
fig_b <- fig_b + draw_plot(plot_overall(grouped_overall,grouped=TRUE),
                           x=.5,y=.06,width=.5,height=.235) +
  draw_label('95% CI: Clopper-Pearson N/(N+S), transformed using the original PAML N:S opportunity ratio.',
             x=0,y=.044,hjust=0,size=8.2) +
  draw_label('Upward arrows: CI exceeds plot range. Stars/ns: dN/dS versus 1. P-trend: separate binned logistic test.',
             x=0,y=.019,hjust=0,size=8.2)
ggsave(file.path(out,'Figure4B_by_lineage_CI_trend.pdf'),fig_b,width=18,height=13.4,device=cairo_pdf,bg='white')
ggsave(file.path(out,'Figure4B_by_lineage_CI_trend.png'),fig_b,width=18,height=13.4,dpi=150,bg='white')

log <- c(paste('Generated:',format(Sys.time(),'%Y-%m-%d %H:%M:%S')),
         'Definition: legacy_all_positive, original 4 bins, original PAML background',
         'CI: Clopper-Pearson then N/S opportunity-ratio transform',
         'Trend: glm(cbind(N,S) ~ VAF_midpoint, family=binomial), one analysis per species',
         'BH family: seven planned species tests within each of All_unique, T_B_NK, Mono',
         'All 28 pooled and 56 grouped N/S rows matched validated inputs; all dN/dS matched.',
         paste('Pooled unique sites:',nrow(site)),
         paste('Lineage-site records used in pooled validation:',nrow(ls)),
         paste('Trend statuses:',paste(capture.output(print(table(trends$analysis,trends$status))),collapse=' | ')))
writeLines(log,file.path(out,'run_log.txt'))
print(trends[,c('analysis','species','slope','se','p_trend','q_BH_across_7_species','n_effective_bins','n_sites','status')],row.names=FALSE)
cat('Output:',out,'\n')
