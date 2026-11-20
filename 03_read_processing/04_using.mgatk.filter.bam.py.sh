#!/bin/bash

set -euo pipefail

SPECIES=(
#    "mouse"
#    "human"
    "rat"
    "chicken"
    "cobra"
    "sturgeon"
    "stingray"
)

REF_GENOMES=(
#    "/PATH/TO/pbmc.mt/ref/mouse/new/mm39_mito.fa"
#    "/PATH/TO/pbmc.mt/SNP/human/mgatk/downstream/rm.germline/human.germline.rm.MT.genome.fa"
    "/PATH/TO/pbmc.mt/SNP/rat/mgatk/downstream/rm.germline/rat.germline.rm.MT.genome.fa"
    "/PATH/TO/pbmc.mt/SNP/chicken/mgatk/downstream/rm.germline/chicken.germline.rm.MT.genome.fa"
    "/PATH/TO/pbmc.mt/SNP/cobra/mgatk/downstream/rm.germline/cobra.germline.rm.MT.genome.fa"
    "/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk/downstream/rm.germline/sturgeon.germline.rm.MT.genome.fa"
    "/PATH/TO/pbmc.mt/SNP/stingray/mgatk/downstream/rm.germline/stingray.germline.rm.MT.genome.fa"
)

if [ ${#SPECIES[@]} -ne ${#REF_GENOMES[@]} ]; then
    echo " !Error: species count(${#SPECIES[@]})and reference genome count(${#REF_GENOMES[@]})do not match！"
    exit 1
fi

CHR="MT"
NM_MAX=2
NH_VALUE=1
PROPER_PAIR="True"

FILTER_SCRIPT="/PATH/TO/pbmc.mt/SNP/pre.process/clean.bam/04_mgatk.filter.bam.py"

for i in "${!SPECIES[@]}"; do
    species="${SPECIES[$i]}"
    ref_fasta="${REF_GENOMES[$i]}"

    echo "=========================================="
    echo "▶ Processing Specie: $species"
    echo "▶ Reference: $ref_fasta"
    echo "=========================================="

    INPUT_BAM="/PATH/TO/pbmc.mt/SNP/$species/mgatk.CleanBam/bam/possorted_bam.chrM.MAPQ40.rmdup.filtered.realign.bam"
    OUTPUT_DIR="/PATH/TO/pbmc.mt/SNP/$species/mgatk.CleanBam/bam/"
    mkdir -p "$OUTPUT_DIR"

    OUTPUT_BAM="$OUTPUT_DIR/$(basename "$INPUT_BAM" .bam).mgatk.fil.bam"
    LOG_FILE="$OUTPUT_DIR/mgatk.bam.filter.log"

    python "$FILTER_SCRIPT" "$INPUT_BAM" "$LOG_FILE" "$CHR" "$PROPER_PAIR" "$NH_VALUE" "$NM_MAX"  > "$OUTPUT_BAM"
    samtools index "$OUTPUT_BAM"


    done

INPUT_BAM="/PATH/TO/pbmc.mt/SNP/human/mgatk.CleanBam/bam/possorted_bam.chrM.MAPQ40.rmdup.filtered.realign.bam"
OUTPUT_DIR="/PATH/TO/pbmc.mt/SNP/human/mgatk.CleanBam/bam/"
mkdir -p "$OUTPUT_DIR"

OUTPUT_BAM="$OUTPUT_DIR/$(basename "$INPUT_BAM" .bam).mgatk.fil.bam"
LOG_FILE="$OUTPUT_DIR/mgatk.bam.filter.log"
python "$FILTER_SCRIPT" "$INPUT_BAM" "$LOG_FILE" "chrM" "$PROPER_PAIR" "$NH_VALUE" "$NM_MAX"  > "$OUTPUT_BAM"
samtools index "$OUTPUT_BAM"



echo "=================All Species Processed===================="

