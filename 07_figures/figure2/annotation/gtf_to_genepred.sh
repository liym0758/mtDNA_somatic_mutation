#!/bin/bash
#RefSeq does not have mitochondria gene definition. So ANNOVAR users need to use either UCSC Known Gene or Ensembl Gene.

#Human
cd /PATH/TO/biosoft/annovar_2020/humandb/
awk '$3 == "chrM"' hg38_ensGene.txt > hg38_chrM_ensGene.txt

#mouse
gtfToGenePred -genePredExt -geneNameAsName2 mm39_MT_genes.gtf mm39_MT_ensGene0.txt
nl mm39_MT_ensGene0.txt > mm39_MT_ensGene.txt

#Rat
cd /PATH/TO/biosoft/annovar_2020/rat_MT_db/
gtfToGenePred -genePredExt -geneNameAsName2 MT_genes.gtf rn7_MT_ensGene0.txt
nl rn7_MT_ensGene0.txt > rn7_MT_ensGene.txt

#Chicken
cd /PATH/TO/biosoft/annovar_2020/chicken_MT_db/
gtfToGenePred -genePredExt -geneNameAsName2 MT_genes.gtf gg1_MT_ensGene0.txt
nl gg1_MT_ensGene0.txt > gg1_MT_ensGene.txt

#Cobra
cd /PATH/TO/biosoft/annovar_2020/cobra_MT_db/
gtfToGenePred -genePredExt -geneNameAsName2 Naja_atra_mt.gtf na_MT_ensGene0.txt
nl na_MT_ensGene0.txt > na_MT_ensGene.txt

#Sturgeon
cd /PATH/TO/biosoft/annovar_2020/sturgeon_MT_db/
gtfToGenePred -genePredExt -geneNameAsName2 MT_genes.gtf as_MT_ensGene0.txt
nl as_MT_ensGene0.txt > as_MT_ensGene.txt

gtfToGenePred -genePredExt -geneNameAsName2 CNCBdb_fixed_Acipenser_sinensis_MT.gtf as_CNCB_MT_ensGene0.txt
nl as_CNCB_MT_ensGene0.txt > as_CNCB_MT_ensGene.txt

#Stingray
cd /PATH/TO/biosoft/annovar_2020/stingray_MT_db/
gtfToGenePred -genePredExt -geneNameAsName2 MT_genes.gtf pl_MT_ensGene0.txt
nl pl_MT_ensGene0.txt > pl_MT_ensGene.txt


