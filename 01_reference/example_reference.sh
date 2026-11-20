#!/usr/bin/env bash
set -euo pipefail

# Example of the shared reference preparation steps. Set paths for one species.
: "${MT_FASTA:?}"
: "${MERGED_FASTA:?}"
: "${NUCLEAR_BOWTIE2_INDEX:?}"
: "${MT_CONTIG:?}"
: "${BLACKLIST_BED:?}"
: "${MASKED_FASTA:?}"
: "${CELLRANGER_ARC_CONFIG:?}"
: "${OUTPUT_DIR:?}"

mkdir -p "$OUTPUT_DIR"/{fastq,bam,peaks}
art_illumina -ss NS50 -c 10000000 -i "$MT_FASTA" -l 20 -o "$OUTPUT_DIR/fastq/mt_"
gzip -f "$OUTPUT_DIR/fastq/mt_.fq"
bowtie2 -p 4 -x "$NUCLEAR_BOWTIE2_INDEX" -U "$OUTPUT_DIR/fastq/mt_.fq.gz" |
  samtools view -bS - | samtools sort -@ 4 -o "$OUTPUT_DIR/bam/all.sorted.bam" -
samtools index "$OUTPUT_DIR/bam/all.sorted.bam"
# The source scripts excluded the mitochondrial contig.
samtools view -H "$OUTPUT_DIR/bam/all.sorted.bam" |
  awk -F '\t' -v mt="$MT_CONTIG" '$1=="@SQ" {sub(/^SN:/,"",$2); if($2!=mt) print $2}' > "$OUTPUT_DIR/nuclear_contigs.txt"
mapfile -t nuclear_contigs < "$OUTPUT_DIR/nuclear_contigs.txt"
samtools view -b -o "$OUTPUT_DIR/bam/nuclear.bam" "$OUTPUT_DIR/bam/all.sorted.bam" "${nuclear_contigs[@]}"
samtools index "$OUTPUT_DIR/bam/nuclear.bam"
macs2 callpeak -t "$OUTPUT_DIR/bam/nuclear.bam" --nomodel --nolambda --keep-dup all -n "$OUTPUT_DIR/peaks/numt"
cat "$BLACKLIST_BED" "$OUTPUT_DIR/peaks/numt_peaks.narrowPeak" |
  awk '{print $1"\t"$2"\t"$3}' | sortBed | mergeBed -i stdin -d 10 > "$OUTPUT_DIR/combined_blacklist.bed"
bedtools maskfasta -fi "$MERGED_FASTA" -bed "$OUTPUT_DIR/combined_blacklist.bed" -fo "$MASKED_FASTA"
bwa index "$MASKED_FASTA"
samtools faidx "$MASKED_FASTA"
cellranger-arc mkref --config "$CELLRANGER_ARC_CONFIG"
