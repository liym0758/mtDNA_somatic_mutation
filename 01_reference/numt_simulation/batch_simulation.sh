#!/bin/bash

set -e
SPECIES=$1
if [[ -z "$SPECIES" ]]; then
  echo "Provide a species name as an argument"
  exit 1
fi

BASE_DIR="/PATH/TO/pbmc.mt/SNP"
READ_LENGTH=150
FRAGMENT_LENGTH=200
TARGET_COVERAGE=500

declare -A MITO_FA_MAP
declare -A NUMT_BED_MAP
declare -A NUC_GENOME_MAP

MITO_FA_MAP[human]="/PATH/TO/pbmc.mt/SNP/human/mgatk/downstream/rm.germline/human.germline.rm.MT.genome.fa"
NUMT_BED_MAP[human]="/PATH/TO/pbmc.mt/ref/human/mt.mask/hg38.full.blacklist.bed"
NUC_GENOME_MAP[human]="/PATH/TO/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa"

MITO_FA_MAP[rat]="/PATH/TO/pbmc.mt/SNP/rat/mgatk/downstream/rm.germline/rat.germline.rm.MT.genome.fa"
NUMT_BED_MAP[rat]="/PATH/TO/pbmc.mt/ref/rat/mt.mask/combinedBlacklist/mRatBN7.2.full.blacklist.bed"
NUC_GENOME_MAP[rat]="/PATH/TO/pbmc.mt/ref/rat/Rat_ref/mRatBN7.2_Ensembl/Rattus_norvegicus_rn7_toplevel/fasta/genome.fa"

MITO_FA_MAP[chicken]="/PATH/TO/pbmc.mt/SNP/chicken/mgatk/downstream/rm.germline/chicken.germline.rm.MT.genome.fa"
NUMT_BED_MAP[chicken]="/PATH/TO/pbmc.mt/ref/chicken/ensembl/mt.mask/combinedBlacklist/bGalGal1.full.blacklist.bed"
NUC_GENOME_MAP[chicken]="/PATH/TO/pbmc.mt/ref/chicken/ensembl/GalGal1_toplevel/fasta/genome.fa"

MITO_FA_MAP[cobra]="/PATH/TO/pbmc.mt/SNP/cobra/mgatk/downstream/rm.germline/cobra.germline.rm.MT.genome.fa"
NUMT_BED_MAP[cobra]="/PATH/TO/pbmc.mt/ref/cobra/Naja_atra/mt.mask/combinedBlacklist/Nana_v5.full.blacklist.bed"
NUC_GENOME_MAP[cobra]="/PATH/TO/pbmc.mt/ref/cobra/Naja_naja.Nana_v5.dna.toplevel.fa"

MITO_FA_MAP[sturgeon]="/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk/downstream/rm.germline/sturgeon.germline.rm.MT.genome.fa"
NUMT_BED_MAP[sturgeon]="/PATH/TO/pbmc.mt/ref/XunFish/NCBI.mt/mt.mask/combinedBlacklist/Acipenser_sinensis.full.blacklist.bed"
NUC_GENOME_MAP[sturgeon]="/PATH/TO/pbmc.mt/ref/XunFish/NCBI.mt/Acipenser_sinensis_toplevel_mt/fasta/genome.fa"

MITO_FA_MAP[stingray]="/PATH/TO/pbmc.mt/SNP/stingray/mgatk/downstream/rm.germline/stingray.germline.rm.MT.genome.fa"
NUMT_BED_MAP[stingray]="/PATH/TO/pbmc.mt/ref/HongFish/NCBI.mt/new/genome.mask/combinedBlacklist/Potamotrygon_leopoldi.full.blacklist.bed"
NUC_GENOME_MAP[stingray]="/PATH/TO/pbmc.mt/ref/HongFish/NCBI.mt/new/genome.mask/nuclear/genome.fa"

MITO_FA="${MITO_FA_MAP[$SPECIES]}"
NUMT_BED="${NUMT_BED_MAP[$SPECIES]}"
NUC_GENOME="${NUC_GENOME_MAP[$SPECIES]}"

SPECIES_DIR="$BASE_DIR/$SPECIES/mgatk/false.positive.test"
mkdir -p "$SPECIES_DIR"/{ref,simulate,align,results,logs}
LOG_FILE="$SPECIES_DIR/logs/${SPECIES}.log"

{
echo "===== Starting species: $SPECIES ====="
echo "Reference file: $MITO_FA"

echo "[$SPECIES] step1: Building NUMT reference sequence..."
bedtools getfasta -fi "$NUC_GENOME" -bed "$NUMT_BED" -fo "$SPECIES_DIR/ref/numt_raw.fa"
awk 'BEGIN {print ">fake_numt"} !/^>/ {printf "%s", $0} END {print ""}' \
  "$SPECIES_DIR/ref/numt_raw.fa" > "$SPECIES_DIR/ref/numt_merged.fa"

NUMT_LEN=$(grep -v ">" "$SPECIES_DIR/ref/numt_merged.fa" | tr -d '\n' | wc -c)
NREADS=$(echo "($TARGET_COVERAGE * $NUMT_LEN) / (2 * $READ_LENGTH)" | bc)
echo "[$SPECIES] NUMT length: $NUMT_LEN bp | target coverage: ${TARGET_COVERAGE}X | required reads: $NREADS"

echo "[$SPECIES] step2: Simulating reads..."
cd "$SPECIES_DIR/simulate"
art_illumina -p -m $FRAGMENT_LENGTH -s 10 -ss HSXn \
  -i "$SPECIES_DIR/ref/numt_merged.fa" \
  -c $NREADS \
  -l $READ_LENGTH \
  -o "$SPECIES.simulated"

echo "[$SPECIES] step3: Aligning to mitochondrial reference..."
bwa index -p mito_index "$MITO_FA" 2>/dev/null
bwa mem -t 8 mito_index \
  "$SPECIES_DIR/simulate/${SPECIES}.simulated1.fq" \
  "$SPECIES_DIR/simulate/${SPECIES}.simulated2.fq" |
  samtools sort -o "$SPECIES_DIR/align/${SPECIES}.numt_mapped.bam"

echo "[$SPECIES] step4: Repairing BAM header..."
samtools view -H "$SPECIES_DIR/align/${SPECIES}.numt_mapped.bam" |
  awk 'NR==1 || !/^@HD/' > "$SPECIES_DIR/align/fixed_header.sam"
samtools reheader "$SPECIES_DIR/align/fixed_header.sam" \
  "$SPECIES_DIR/align/${SPECIES}.numt_mapped.bam" > \
  "$SPECIES_DIR/align/${SPECIES}.fixed.bam"
samtools index "$SPECIES_DIR/align/${SPECIES}.fixed.bam"

echo "[$SPECIES] step5: Variant calling and filtering..."
freebayes -f "$MITO_FA" \
  "$SPECIES_DIR/align/${SPECIES}.fixed.bam" \
  --min-mapping-quality 59 \
  --min-alternate-count 2 \
  --min-coverage 10 > \
  "$SPECIES_DIR/results/${SPECIES}.raw.vcf"

bcftools filter -O v \
  -o "$SPECIES_DIR/results/${SPECIES}.filtered.vcf" \
  -i 'QUAL > 20 && INFO/DP > 5 && INFO/AO > 2 && SAF > 0 && SAR > 0' \
  "$SPECIES_DIR/results/${SPECIES}.raw.vcf"

echo "===== species $SPECIES processing complete! ====="
echo "Result file: $SPECIES_DIR/results/${SPECIES}.filtered.vcf"
} | tee "$LOG_FILE"

