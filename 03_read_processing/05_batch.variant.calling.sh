#!/bin/bash

set -euo pipefail

SPECIES=(
    "mouse"
    #"human"
    #"rat"
    #"chicken"
    #"cobra"
    #"sturgeon"
    #"stingray"
)

REF_GENOMES=(
    "/PATH/TO/pbmc.mt/ref/mouse/new/mm39_mito.fa"
    #"/PATH/TO/pbmc.mt/SNP/human/mgatk/downstream/rm.germline/human.germline.rm.MT.genome.fa"
    #"/PATH/TO/pbmc.mt/SNP/rat/mgatk/downstream/rm.germline/rat.germline.rm.MT.genome.fa"
    #"/PATH/TO/pbmc.mt/SNP/chicken/mgatk/downstream/rm.germline/chicken.germline.rm.MT.genome.fa"
    #"/PATH/TO/pbmc.mt/SNP/cobra/mgatk/downstream/rm.germline/cobra.germline.rm.MT.genome.fa"
    #"/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk/downstream/rm.germline/sturgeon.germline.rm.MT.genome.fa"
    #"/PATH/TO/pbmc.mt/SNP/stingray/mgatk/downstream/rm.germline/stingray.germline.rm.MT.genome.fa"
)

BASE_QUAL=30
MAPQ=30
MPILEUP_DEPTH=0
MIN_DEPTH=10
MIN_SUPPORTING_READS=2


VARSCAN_JAR="/PATH/TO/biosoft/MitoSort/VarScan.v2.4.6.jar"

for i in "${!SPECIES[@]}"; do
    species="${SPECIES[$i]}"
    ref_fasta="${REF_GENOMES[$i]}"

    echo "=============================="
    echo "▶ Processing: $species"
    echo "▶ Reference:  $ref_fasta"
    echo "=============================="

    INPUT_BAM="/PATH/TO/pbmc.mt/SNP/$species/mgatk.bam.QC/bam/possorted_bam.chrM.MAPQ40.rmdup.filtered.realign.mgatk.fil.bam"
    OUTPUT_DIR="/PATH/TO/pbmc.mt/SNP/$species/mgatk.bam.QC/bam/bulk.var.calling"
    mkdir -p "$OUTPUT_DIR"

    MPILEUP_FILE="$OUTPUT_DIR/$species.mpileup.txt"
    VCF_FILE="$OUTPUT_DIR/$species.snp.vcf"

    echo "▶ Generating mpileup..."
    samtools mpileup -f "$ref_fasta" -Q "$BASE_QUAL" -q "$MAPQ" -d "$MPILEUP_DEPTH" "$INPUT_BAM" > "$MPILEUP_FILE"

    echo "▶ Calling SNPs with VarScan..."

    java -jar "$VARSCAN_JAR" pileup2snp --min-avg-qual "$BASE_QUAL" --min-coverage "$MIN_DEPTH"  --min-ads2 "$MIN_SUPPORTING_READS" --min-var-freq 0 --p-value 0.05 < "$MPILEUP_FILE" > "$VCF_FILE"

    echo "▶ Finished: $VCF_FILE"
done

echo "================= All Species Processed ================="
