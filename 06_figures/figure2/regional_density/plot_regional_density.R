#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

job <- Sys.getenv('MITOEVO_REGION_DENSITY_ROOT')
if(!nzchar(job) || !dir.exists(job)) stop('Set MITOEVO_REGION_DENSITY_ROOT to the regional density run')
infile <- file.path(job, 'Tables/regional_density_N_L.tsv')
out_pdf <- file.path(job, 'Figures/regional_mutation_density_five_celltypes_7x1_freeY.pdf')
out_png <- file.path(job, 'Figures/regional_mutation_density_five_celltypes_7x1_freeY.png')

# Plot-only script: read the finalized N/L/density table and do not recalculate it.
d <- fread(infile)
species_order <- c('human','mouse','rat','chicken','cobra','sturgeon','stingray')
region_order <- c('CDS','rRNA','tRNA','Noncoding')
stopifnot(nrow(d) == 28L, setequal(d$species, species_order), setequal(d$region, region_order))
d[, species := factor(species, levels=species_order,
  labels=c('Human','Mouse','Rat','Chicken','Cobra','Sturgeon','Stingray'))]
d[, region := factor(region, levels=region_order)]
d[, label := paste0('N=', N, '\nL=', L)]
d[, low_coverage := as.character(species)=='Sturgeon' & as.character(region)=='Noncoding' & L < 100L]
cols <- c(CDS='#625D9E', rRNA='#B3659F', tRNA='#C7B1D4', Noncoding='#D3D3D3')

p <- ggplot(d, aes(x=region, y=density_sites_per_kb, fill=region)) +
  geom_col(color='black', linewidth=.35, width=.76) +
  geom_text(aes(label=label), vjust=-.25, size=2.5, lineheight=.9) +
  # Place the coverage warning at the top-left of the Sturgeon facet, separate from its N/L label.
  geom_text(data=data.table(species=factor('Sturgeon', levels=levels(d$species)),
      region=factor('CDS', levels=region_order), y=Inf, note='Coverage limited'),
    aes(x=region, y=y, label=note), inherit.aes=FALSE, hjust=0, vjust=1.35,
    size=2.6, color='#8B1A1A', fontface='bold') +
  facet_wrap(~species, ncol=1, scales='free_y') +
  scale_fill_manual(values=cols, drop=FALSE) +
  scale_x_discrete(drop=FALSE) +
  scale_y_continuous(limits=c(0, NA), expand=expansion(mult=c(0,.30))) +
  labs(x=NULL, y='Regional mutation density (mutated positions/kb)',
    title='Regional mutation density',
    subtitle='N = unique positive-VAF sites in five cell types
L = unmasked positions with at least one nonzero coverage record',
    caption='Sturgeon Noncoding coverage limited; N and density retained.') +
  theme_bw(base_size=8.5) +
  theme(panel.background=element_rect(fill='white', colour=NA),
    plot.background=element_rect(fill='white', colour=NA),
    panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
    strip.background=element_rect(fill='white', colour='black'),
    strip.text=element_text(face='bold', size=9),
    axis.text.x=element_text(size=8), axis.text.y=element_text(size=7),
    axis.title.y=element_text(size=9), legend.position='none',
    plot.title=element_text(face='bold', hjust=.5, size=10),
    plot.subtitle=element_text(size=6.5), plot.caption=element_text(size=5.8, hjust=0),
    panel.spacing=grid::unit(.28,'lines'), plot.margin=margin(5,6,5,5))

ggsave(out_pdf, p, width=4, height=19, device=cairo_pdf, limitsize=FALSE)
ggsave(out_png, p, width=4, height=19, dpi=320, bg='white', limitsize=FALSE)
cat('Input rows:', nrow(d), '\n')
cat('Outputs:\n', out_pdf, '\n', out_png, '\n')


