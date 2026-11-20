#!/bin/bash

#Mouse
cd /PATH/TO/biosoft/annovar_2020/mm39_MT_db/

perl /PATH/TO/biosoft/annovar_2020/retrieve_seq_from_fasta.pl \
--format ensGene \
--seqfile mm39_mito.fa mm39_MT_ensGene.txt \
--out mm39_MT_ensGeneMrna0.fa

sed 's/Warning: this coding transcript does not have correct ORF annotation\.//g' \
  mm39_MT_ensGeneMrna0.fa > mm39_MT_ensGeneMrna.fa


cd /PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/mouse/

perl /PATH/TO/biosoft/annovar_2020/table_annovar.pl \
  "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/mouse/mouse.AnnovarInput" \
  "/PATH/TO/biosoft/annovar_2020/mm39_MT_db/" \
  -buildver mm39_MT \
  -out result/mm39_anno \
  -protocol ensGene \
  -operation g \
  -remove -polish -nastring .



#Human
cd /PATH/TO/biosoft/annovar_2020/humandb/
mv hg38_chrM_ensGeneMrna.fa hg38_chrM_ensGeneMrna0.fa
sed 's/Warning: this coding transcript does not have correct ORF annotation\.//g' \
  hg38_chrM_ensGeneMrna0.fa > hg38_chrM_ensGeneMrna.fa

cd /PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/human/
perl /PATH/TO/biosoft/annovar_2020/table_annovar.pl \
  "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/human/human.AnnovarInput" \
  "/PATH/TO/biosoft/annovar_2020/humandb/" \
  -buildver hg38_chrM \
  -out result/hg38_anno \
  -protocol ensGene \
  -operation g \
  -remove -polish -nastring .


#Rat
#Mrna
cd /PATH/TO/biosoft/annovar_2020/rat_MT_db/

perl /PATH/TO/biosoft/annovar_2020/retrieve_seq_from_fasta.pl \
--format ensGene \
--seqfile rat.germline.rm.MT.genome.fa rn7_MT_ensGene.txt \
--out rn7_MT_ensGeneMrna0.fa 

sed 's/Warning: this coding transcript does not have correct ORF annotation\.//g' \
  rn7_MT_ensGeneMrna0.fa > rn7_MT_ensGeneMrna.fa

#Annotation
cd /PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/rat/

perl /PATH/TO/biosoft/annovar_2020/table_annovar.pl \
  "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/rat/rat.AnnovarInput" \
  "/PATH/TO/biosoft/annovar_2020/rat_MT_db/" \
  -buildver rn7_MT \
  -out result/rn7_MT_anno \
  -protocol ensGene \
  -operation g \
  -remove -polish -nastring .

#Chicken
#Mrna
cd /PATH/TO/biosoft/annovar_2020/chicken_MT_db/

perl /PATH/TO/biosoft/annovar_2020/retrieve_seq_from_fasta.pl \
--format ensGene \
--seqfile chicken.germline.rm.MT.genome.fa gg1_MT_ensGene.txt \
--out gg1_MT_ensGeneMrna0.fa 

sed 's/Warning: this coding transcript does not have correct ORF annotation\.//g' \
  gg1_MT_ensGeneMrna0.fa > gg1_MT_ensGeneMrna.fa

#Annotation
cd /PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/chicken/

perl /PATH/TO/biosoft/annovar_2020/table_annovar.pl \
  "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/chicken/chicken.AnnovarInput" \
  "/PATH/TO/biosoft/annovar_2020/chicken_MT_db/" \
  -buildver gg1_MT \
  -out result/gg1_MT_anno \
  -protocol ensGene \
  -operation g \
  -remove -polish -nastring .

#Cobra
#Mrna
cd /PATH/TO/biosoft/annovar_2020/cobra_MT_db/

perl /PATH/TO/biosoft/annovar_2020/retrieve_seq_from_fasta.pl \
--format ensGene \
--seqfile cobra.germline.rm.MT.genome.fa na_MT_ensGene.txt \
--out na_MT_ensGeneMrna0.fa 

sed 's/Warning: this coding transcript does not have correct ORF annotation\.//g' \
  na_MT_ensGeneMrna0.fa > na_MT_ensGeneMrna.fa

#Annotation
cd /PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/cobra/

perl /PATH/TO/biosoft/annovar_2020/table_annovar.pl \
  "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/cobra/cobra.AnnovarInput" \
  "/PATH/TO/biosoft/annovar_2020/cobra_MT_db/" \
  -buildver na_MT \
  -out result/na_MT_anno \
  -protocol ensGene \
  -operation g \
  -remove -polish -nastring .

#Sturgeon
#Mrna
cd /PATH/TO/biosoft/annovar_2020/sturgeon_MT_db/

perl /PATH/TO/biosoft/annovar_2020/retrieve_seq_from_fasta.pl \
--format ensGene \
--seqfile sturgeon.germline.rm.MT.genome.fa as_CNCB_MT_ensGene.txt \
--out as_CNCB_MT_ensGeneMrna0.fa 

sed 's/Warning: this coding transcript does not have correct ORF annotation\.//g' \
  as_CNCB_MT_ensGeneMrna0.fa > as_CNCB_MT_ensGeneMrna.fa

#Annotation
cd /PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/sturgeon/

perl /PATH/TO/biosoft/annovar_2020/table_annovar.pl \
  "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/sturgeon/sturgeon.AnnovarInput" \
  "/PATH/TO/biosoft/annovar_2020/sturgeon_MT_db/" \
  -buildver as_CNCB_MT \
  -out result/as_CNCB_MT_anno \
  -protocol ensGene \
  -operation g \
  -remove -polish -nastring .

#Stingray
#Mrna
cd /PATH/TO/biosoft/annovar_2020/stingray_MT_db/

perl /PATH/TO/biosoft/annovar_2020/retrieve_seq_from_fasta.pl \
--format ensGene \
--seqfile stingray.germline.rm.MT.genome.fa pl_MT_ensGene.txt \
--out pl_MT_ensGeneMrna0.fa 

sed 's/Warning: this coding transcript does not have correct ORF annotation\.//g' \
  pl_MT_ensGeneMrna0.fa > pl_MT_ensGeneMrna.fa

#Annotation
cd /PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/stingray/

perl /PATH/TO/biosoft/annovar_2020/table_annovar.pl \
  "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/stingray/stingray.AnnovarInput" \
  "/PATH/TO/biosoft/annovar_2020/stingray_MT_db/" \
  -buildver pl_MT \
  -out result/pl_MT_anno \
  -protocol ensGene \
  -operation g \
  -remove -polish -nastring .



