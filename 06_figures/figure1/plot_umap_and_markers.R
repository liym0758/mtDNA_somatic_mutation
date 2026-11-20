#!/usr/bin/env Rscript
# Plot-only revision. Reads existing annotation and embeddings; never saves inputs.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(ggplot2); library(Matrix); library(ggrepel)})
sp <- tolower(commandArgs(trailingOnly=TRUE)[1])
species <- c("human","mouse","rat","chicken","cobra","sturgeon","stingray")
stopifnot(length(sp)==1, sp %in% species)
root <- Sys.getenv("MITOEVO_ARCHR_ROOT")
if (!nzchar(root)) stop("Set MITOEVO_ARCHR_ROOT to the local annotation root")
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
if (length(script_arg)!=1L) stop("Run with Rscript")
code <- dirname(normalizePath(sub("^--file=", "", script_arg), winslash="/", mustWork=TRUE))
out <- Sys.getenv("MITOEVO_FIG1_OUTPUT_ROOT")
if (!nzchar(out)) stop("Set MITOEVO_FIG1_OUTPUT_ROOT to a new output directory")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
plan <- read.delim(file.path(code,"marker_plan_v3.tsv"),check.names=FALSE)
plan <- plan[plan$species==sp,,drop=FALSE]
if (!nrow(plan)) stop("No marker plan for ",sp)
title <- paste0(toupper(substr(sp,1,1)),substring(sp,2))
ct_order <- list(
 human=c("Naive B","Memory B","Plasma B","CD4+ T","CD8+ T","NK cell","CD14+ Monocyte","CD16+ Monocyte","pDC"),
 mouse=c("B cell","CD4+ T","CD8+ T","NKT","NK cell","Monocyte","Neutrophil"),
 rat=c("B cell","CD4+ T","CD8+ T","NK cell","CD14+ Monocyte","CD16+ Monocyte"),
 chicken=c("Naive B","Memory B","CD4+ T","CD8+ T","Delta T","NK cell","Monocyte","Platelet"),
 cobra=c("B cell","T/NK cell","Monocyte_1","Monocyte_2","Platelet","Erythrocyte"),
 sturgeon=c("B cell","T/NK cell","Monocyte","Neutrophil","Platelet"),
 stingray=c("B cell","T/NK cell","Monocyte","Neutrophil","Platelet","Erythrocyte"))
colors <- c("B cell"="#D490B9","Naive B"="#E5D2DD","Memory B"="#B3659F","Plasma B"="#853280",
 "CD4+ T"="#ACBDDE","CD4+ Naive T"="#CCE0F5","CD4+ Memory T"="#8291B0",
 "CD8+ T"="#AF95C4","CD8+ Naive T"="#CCC9E6","CD8+ Memory T"="#8C549C",
 "NK cell"="#A63376","T/NK cell"="#8E90B9","Delta T"="#4F418E",,
 "NKT"="#5F3D69","Monocyte"="#7EA795","Monocyte_1"="#7EA475","Monocyte_2"="#3B6E4F",
 "CD14+ Monocyte"="#C4DAB0","CD16+ Monocyte"="#7EA475","Neutrophil"="#274B2F",
 "pDC"="#C6902E","cDC"="#DEA091","Platelet"="#D1B56C","Erythrocyte"="#712820",
 "HSP"="#BF5D58","Myeloid Progenitor"="#D5B9B2","Early Baso"="#2C575E","Others"="#959492")
global_key <- read.delim(file.path(code,"global_umap_number_key.tsv"),check.names=FALSE)
colors <- setNames(global_key$color,global_key$sub_celltype)
norm <- function(x) toupper(sub("\\.[0-9]+$","",as.character(x)))
collapse_ct <- function(x) {
 x <- as.character(x)
 x[x %in% c("CD4+ Naive T","CD4+ Memory T")] <- "CD4+ T"
 x[x %in% c("CD8+ Naive T","CD8+ Memory T")] <- "CD8+ T"
 x
}
logmsg <- function(...) cat(format(Sys.time(),"%F %T")," ",sp,": ",paste0(...,collapse=""),"\n",sep="",file=stderr())

if (sp=="human") {
 suppressPackageStartupMessages({library(ArchR);library(SummarizedExperiment)})
 fullpath <- "/PATH/TO/pbmc.mt/data/human_bottleneck/output_pbmc/proj_res1_annotated/lym_anno-ArchR-Project.rds"
 pbmcpath <- "/PATH/TO/pbmc.mt/data/human_bottleneck/output_pbmc/proj_PBMC_only/Save-ArchR-Project.rds"
 logmsg("Reading PBMC-only ArchR project for saved UMAP")
 pbmc <- readRDS(pbmcpath)
 um <- as.data.frame(ArchR::getEmbedding(pbmc,embedding="UMAP"))
 colnames(um)[1:2] <- c("x","y")
 md <- as.data.frame(pbmc@cellColData)
 common <- intersect(rownames(um),rownames(md))
 if (!length(common)) stop("Human PBMC UMAP and metadata cell IDs do not overlap")
 um <- data.frame(x=um[common,"x"],y=um[common,"y"],celltype=collapse_ct(md[common,"sub_celltype"]),row.names=common)
 rm(pbmc,md); gc(verbose=FALSE)
 logmsg("Reading annotated ArchR project for existing GeneScoreMatrix")
 proj <- readRDS(fullpath)
 rec <- as.character(proj@sampleColData$ArrowFiles)
 if (!all(file.exists(rec))) {
   repl <- file.path(dirname(fullpath),"ArrowFiles",basename(rec))
   if (!all(file.exists(repl))) stop("Human ArrowFiles missing")
   proj@sampleColData$ArrowFiles <- normalizePath(repl,winslash="/",mustWork=TRUE)
 }
 meta <- as.data.frame(proj@cellColData)
 if (!all(c("Clusters2","sub_celltype") %in% names(meta))) stop("Human label fields missing")
 hm <- as.data.frame(table(Clusters2=as.character(meta$Clusters2),sub_celltype=as.character(meta$sub_celltype)),stringsAsFactors=FALSE)
 hm <- hm[hm$Freq>0,,drop=FALSE]
 names(hm)[3] <- "n_cells"
 write.table(hm,file.path(out,"human_observed_label_mapping.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
 logmsg("Loading GeneScoreMatrix")
 se <- ArchR::getMatrixFromProject(proj,useMatrix="GeneScoreMatrix")
 rd <- as.data.frame(SummarizedExperiment::rowData(se))
 syms <- if ("name" %in% names(rd)) as.character(rd$name) else if ("gene_name" %in% names(rd)) as.character(rd$gene_name) else rownames(se)
 feat <- rownames(se)
 if (is.null(feat)) feat <- paste0("r",seq_len(nrow(se)))
 ix <- which(norm(syms) %in% norm(plan$gene))
 allcells <- intersect(colnames(se),rownames(meta))
 allcells <- allcells[as.character(meta[allcells,"Sample"]) %in% c("PBMC_rep1","PBMC_rep2")]
 if (!length(allcells)) stop("No PBMC_rep1/2 cells in Human GeneScoreMatrix")
 groups <- intersect(ct_order[[sp]],unique(collapse_ct(meta[allcells,"sub_celltype"])))
 cells <- allcells[collapse_ct(meta[allcells,"sub_celltype"]) %in% groups]
 x <- SummarizedExperiment::assay(se[ix,cells,drop=FALSE])
 rownames(x) <- feat[ix]
 genes <- syms[ix]
 labels <- collapse_ct(meta[cells,"sub_celltype"])
 rm(se,proj,meta,rd); gc(verbose=FALSE)
} else {
 suppressPackageStartupMessages({library(Seurat);library(Signac)})
 path <- file.path(root,"RDS",paste0(sp,".ATAC.RDS"))
 logmsg("Reading ",path)
 obj <- readRDS(path)
 md <- obj@meta.data
 if (!"sub_celltype" %in% names(md)) stop("sub_celltype missing")
 em <- as.data.frame(SeuratObject::Embeddings(obj[["umap"]]))
 colnames(em)[1:2] <- c("x","y")
 common <- intersect(rownames(em),rownames(md))
 if (!length(common)) stop("UMAP and metadata cell IDs do not overlap")
 um <- data.frame(x=em[common,"x"],y=em[common,"y"],celltype=collapse_ct(md[common,"sub_celltype"]),row.names=common)
 labels_all <- collapse_ct(md$sub_celltype)
 groups <- intersect(ct_order[[sp]],unique(labels_all))
 groups <- c(groups,setdiff(sort(unique(labels_all)),c(groups,NA_character_)))
 cells <- rownames(md)[labels_all %in% groups]
 xa <- NULL; used_layer <- NULL
 for (ly in c("counts","data")) {
   z <- tryCatch(SeuratObject::GetAssayData(obj,assay="ACTIVITY",layer=ly),error=function(e) NULL)
   if (!is.null(z) && nrow(z)>0) {xa <- z; used_layer <- ly; break}
 }
 if (is.null(xa)) stop("ACTIVITY assay unavailable")
 cells <- intersect(cells,colnames(xa))
 ix <- which(norm(rownames(xa)) %in% norm(plan$gene))
 x <- xa[ix,cells,drop=FALSE]
 genes <- rownames(x)
 labels <- collapse_ct(md[cells,"sub_celltype"])
 rm(obj,xa,md); gc(verbose=FALSE)
}
if (length(groups)<2 || !nrow(x) || !ncol(x)) stop("Insufficient matrix/cell data")
logmsg("Summarizing ",nrow(x)," genes across ",ncol(x)," cells in ",length(groups)," subtypes")
g <- match(labels,groups)
H <- Matrix::sparseMatrix(i=seq_along(g),j=g,x=1,dims=c(length(g),length(groups)))
n <- tabulate(g,length(groups))
if (inherits(x,"sparseMatrix")) {
 y <- as(x,"dgCMatrix")
 if (sp=="human" || used_layer=="counts") y@x <- log1p(y@x)
} else y <- if (sp=="human" || used_layer=="counts") log1p(x) else x
avg <- sweep(as.matrix(y %*% H),2,n,"/")
pct <- sweep(as.matrix((x>0) %*% H),2,n,"/")
rownames(avg) <- rownames(pct) <- rownames(x)
colnames(avg) <- colnames(pct) <- groups
rm(x,y,H); gc(verbose=FALSE)

# Candidate order is biological priority from marker_plan_v3.tsv. Weak but plausible
# lineage markers are retained when two subtype-specific features do not exist.
found <- logical(nrow(plan))
res <- vector('list',nrow(plan))
for (i in seq_len(nrow(plan))) {
 target <- plan$celltype[i]; gene <- plan$gene[i]
 if (!target %in% groups) next
 ii <- which(norm(genes)==norm(gene) | norm(rownames(avg))==norm(gene))
 if (!length(ii)) next
 found[i] <- TRUE
 others <- setdiff(groups,target)
 score <- avg[ii,target]-rowMeans(avg[ii,others,drop=FALSE]) +
          0.35*(pct[ii,target]-rowMeans(pct[ii,others,drop=FALSE]))
 ii <- ii[which.max(score)]
 othermax <- max(avg[ii,others]); pctmax <- max(pct[ii,others])
 res[[i]] <- data.frame(species=sp,celltype=target,feature_id=rownames(avg)[ii],
  gene_symbol=if (sp=='human') genes[ii] else if (sp=='cobra' && gene=='ENSNNAG00000002567') 'TUBB1' else if (grepl('^(ENSNNAG|GWHGAOTN|LOC)',gene)) '' else gene,
  plot_label=if (sp=='cobra' && gene=='ENSNNAG00000002567') paste0('TUBB1 (',gene,')') else gene,
  biological_basis=plan$biological_basis[i],mapping_source=plan$mapping_source[i],
  assessment=plan$assessment[i],pct_in_type=100*pct[ii,target],
  pct_max_other=100*pctmax,mean_log1p_in_type=avg[ii,target],
  mean_log1p_max_other=othermax,delta_mean_max=avg[ii,target]-othermax,
  relative_score=score[which.max(score)],stringsAsFactors=FALSE)
}
availability <- cbind(plan,present_in_matrix=found)
write.table(availability,file.path(out,paste0(sp,'.candidate_availability.tsv')),
 sep='\t',quote=FALSE,row.names=FALSE)
tab <- do.call(rbind,Filter(Negate(is.null),res))
if (is.null(tab) || !nrow(tab)) stop('No candidate genes matched the matrix')
tab$signal_status <- ifelse(tab$pct_in_type<1 | tab$mean_log1p_in_type<=0,'no_usable_signal',
 ifelse(tab$relative_score<=0,'not_enriched',ifelse(tab$delta_mean_max<=0,'shared_or_weak','relative_enriched')))
tab$display <- tab$signal_status!='no_usable_signal'
tab$display[duplicated(paste(tab$celltype,tab$feature_id))] <- FALSE
cap <- ifelse(tab$celltype %in% c('CD4+ T','CD8+ T'),4L,2L)
tab$display <- tab$display & ave(as.integer(tab$display),tab$celltype,FUN=cumsum)<=cap
write.table(tab,file.path(out,paste0(sp,'.marker_evidence.tsv')),
 sep='\t',quote=FALSE,row.names=FALSE,na='')
coverage <- data.frame(species=sp,celltype=groups,n_cells=n,
 n_displayed=vapply(groups,function(ct) sum(tab$display & tab$celltype==ct),integer(1)),stringsAsFactors=FALSE)
coverage$interpretation <- ifelse(coverage$celltype=='Others','unassigned_category',
 ifelse(coverage$n_displayed<2,'fewer_than_two_usable_candidates','two_candidates_check_evidence'))
write.table(coverage,file.path(out,paste0(sp,'.subtype_coverage.tsv')),
 sep='\t',quote=FALSE,row.names=FALSE)
sel <- tab[tab$display,,drop=FALSE]
if (!nrow(sel)) stop('No displayed candidate')
bub <- do.call(rbind,lapply(seq_len(nrow(sel)),function(i) {
 ii <- match(sel$feature_id[i],rownames(avg))
 a <- avg[ii,groups]
 z <- if (sd(a)>0) as.numeric(scale(a)) else rep(0,length(a))
 data.frame(marker=paste0('m',i),marker_group=sel$celltype[i],celltype=groups,
  scaled=pmax(-1,pmin(2,z)),pct=100*pct[ii,groups])
}))
bub$marker <- factor(bub$marker,levels=paste0('m',seq_len(nrow(sel))))
bub$marker_group <- factor(bub$marker_group,levels=groups)
bub$celltype <- factor(bub$celltype,levels=groups)
marker_labs <- setNames(sel$plot_label,paste0('m',seq_len(nrow(sel))))
pale <- function(f) {v <- grDevices::col2rgb('#2C4A7F')[,1]/255; grDevices::rgb(v[1]+(1-v[1])*f,v[2]+(1-v[2])*f,v[3]+(1-v[3])*f)}
pct_cap <- min(100,max(10,ceiling(max(bub$pct)/10)*10))
pct_breaks <- pretty(c(0,pct_cap),n=3)
pct_breaks <- pct_breaks[pct_breaks>0 & pct_breaks<=pct_cap]
bp <- ggplot(bub,aes(celltype,marker)) + geom_point(aes(size=pct,color=scaled)) +
 facet_grid(marker_group~.,scales='free_y',space='free_y',drop=TRUE) +
 scale_y_discrete(labels=local({v <- marker_labs; function(z) unname(v[z])})) +
 scale_color_gradientn(colors=c(pale(.88),pale(.5),'#2C4A7F'),limits=c(-1,2),
  oob=scales::squish,name='Scaled activity') +
 scale_size_continuous(range=c(.4,3.7),limits=c(0,pct_cap),breaks=pct_breaks,
  name='Accessible cells (%)') +
 labs(title=paste('ATAC Bubble of',title,'PBMC'),x=NULL,y=NULL,
  caption=if(any(sel$assessment=='unverified_identity')) 'Unverified gene IDs: quantitative candidates only' else NULL) +
 theme_classic(base_size=11) + theme(plot.title=element_text(hjust=.5,face='bold',size=13),
  axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=8,color='black'),
  axis.text.y=element_text(size=7,color='black'),
  strip.background=element_rect(fill='grey96',color='grey80'),
  strip.text.y=element_text(size=8,face='bold',angle=0),
  panel.spacing.y=grid::unit(.05,'cm'),
  plot.margin=margin(8,8,7,8,'pt'),legend.position='bottom',
  legend.box='horizontal',legend.title=element_text(size=8),legend.text=element_text(size=8),
  plot.caption=element_text(size=7,hjust=0)) +
 guides(size=guide_legend(order=1,nrow=1),
        color=guide_colorbar(order=2,barwidth=grid::unit(2,'cm'),barheight=grid::unit(.25,'cm')))
invisible(ggplot_build(bp))
ggsave(file.path(out,paste0(sp,'.bubble.pdf')),bp,width=6,height=8,units='in',device='pdf',useDingbats=FALSE)
saveRDS(bp,file.path(out,paste0(sp,'.bubble.plot.rds')))

um <- um[!is.na(um$celltype) & nzchar(um$celltype),,drop=FALSE]
uorder <- global_key$sub_celltype[global_key$sub_celltype %in% unique(um$celltype)]
if (length(uorder)!=length(unique(um$celltype))) stop('UMAP labels missing from global legend')
um$celltype <- factor(um$celltype,levels=uorder)
key <- global_key[match(uorder,global_key$sub_celltype),,drop=FALSE]
key <- data.frame(species=sp,key,check.names=FALSE)
write.table(key,file.path(out,paste0(sp,'.umap_number_key.tsv')),
 sep='\t',quote=FALSE,row.names=FALSE)
anchors <- do.call(rbind,lapply(seq_along(uorder),function(i) {
 d <- um[um$celltype==uorder[i],c('x','y')]
 center <- c(median(d$x),median(d$y))
 dx <- (d$x-center[1])/diff(range(um$x)); dy <- (d$y-center[2])/diff(range(um$y))
 j <- which.min(dx^2+dy^2)
 data.frame(x=d$x[j],y=d$y[j],number=as.character(key$number[i]),celltype=uorder[i])
}))
up <- ggplot(um,aes(x,y,color=celltype)) + geom_point(size=.3,alpha=.85) +
 ggrepel::geom_text_repel(data=anchors,aes(x=x,y=y,label=number),inherit.aes=FALSE,
  color='black',size=5,fontface='bold',segment.color=NA,
  box.padding=.1,point.padding=.05,min.segment.length=0,
  max.overlaps=Inf,seed=128,force=.4,force_pull=3) +
 scale_color_manual(values=colors[uorder],breaks=uorder,drop=FALSE,guide='none') +
 labs(title=paste('ATAC UMAP of',title,'PBMC'),x='UMAP 1',y='UMAP 2') +
 coord_equal() + theme_classic(base_size=20) +
 theme(axis.line=element_line(size=.5),axis.ticks=element_line(size=.5),
  plot.title=element_text(face='bold'),axis.text=element_text(size=12),
  plot.margin=margin(10,10,10,10,'pt'))
invisible(ggplot_build(up))
ggsave(file.path(out,paste0(sp,'.umap_numbered.pdf')),up,width=6,height=8,units='in',device='pdf',useDingbats=FALSE)
saveRDS(up,file.path(out,paste0(sp,'.umap.plot.rds')))
logmsg('DONE: ',nrow(um),' UMAP cells; ',length(uorder),' subtypes; ',nrow(sel),' displayed candidates')
