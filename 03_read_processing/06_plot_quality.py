import pysam
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import argparse
import os

## Usage ==============================================
# python 06_plot_quality.py your_file.bam --out ./figs/quality.png
## ====================================================

def plot_base_quality(bam_file, output_png):
    print(f"\n📂 Processing BAM file: {bam_file}")
    
    bam = pysam.AlignmentFile(bam_file, "rb")
    qual_counter = Counter()
    total_reads = 0
    total_bases = 0
    unmapped_reads = 0

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            unmapped_reads += 1
            continue
        quals = read.query_qualities
        if quals is None:
            continue
        total_reads += 1
        total_bases += len(quals)
        for q in quals:
            qual_counter[q] += 1

    bam.close()

    print(f"Total reads processed: {total_reads}")
    print(f"Unmapped reads skipped: {unmapped_reads}")
    print(f"Total bases (with quality): {total_bases}")
    print("\n Quality counts:")
    for q in sorted(qual_counter):
        print(f"  Q{q}: {qual_counter[q]}")
    print()

    if not qual_counter:
        print("❌ No quality scores found! Check if the BAM file has mapped reads with qualities.")
        return

    qualities = sorted(qual_counter.keys())
    counts = [qual_counter[q] for q in qualities]

    sns.set(style="whitegrid")
    plt.figure(figsize=(12, 6))
    sns.barplot(x=qualities, y=counts, color='steelblue')

    plt.title("Base Quality Distribution")
    plt.xlabel("Phred Quality Score")
    plt.ylabel("Base Count")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(output_png)
    print(f"Plot saved to: {output_png}")

def main():
    parser = argparse.ArgumentParser(description="Plot base quality distribution from a BAM file.")
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("--out", default="base_quality_distribution.png", help="Output image file (default: base_quality_distribution.png)")
    args = parser.parse_args()

    out_dir = os.path.dirname(args.out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    plot_base_quality(args.bam, args.out)

if __name__ == "__main__":
    main()

