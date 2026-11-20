import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse

## USAGE: =======================================================================
# python visualize_mpileup.py \
#  --sample Sample1 --file /data/Sample1.mpileup.Statistics.tsv \
#  --sample Sample2 --file /data/Sample2.mpileup.Statistics.tsv \
#  --outdir ./my_output_plots
## =============================================================================

parser = argparse.ArgumentParser(description="Generate base statistics visualization from mpileup files.")
parser.add_argument('--sample', action='append', required=True, help='Sample name. Can be specified multiple times.')
parser.add_argument('--file', action='append', required=True, help='Path to corresponding TSV file. Order must match --sample.')
parser.add_argument('--outdir', required=True, help='Directory to save output plots.')

args = parser.parse_args()

if len(args.sample) != len(args.file):
    raise ValueError("Number of --sample and --file arguments must match!")

sample_to_file = dict(zip(args.sample, args.file))

output_dir = args.outdir
os.makedirs(output_dir, exist_ok=True)

def process_sample(sample, file_path):
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return
    
    df = pd.read_csv(file_path, sep='\t')
    valid_bases = ['A', 'T', 'C', 'G']
    df = df[df['base'].isin(valid_bases)]
    
    grouped = df.groupby('base').agg({
        'total_count': 'sum',
        'mean_base_quality': lambda x: np.average(x, weights=df.loc[x.index, 'total_count']),
        'pos_strand_count': 'sum',
        'neg_strand_count': 'sum'
    }).reset_index()
    
    plt.figure(figsize=(18, 5))
    plt.suptitle(f"Analysis of {sample}", fontsize=16)
    
    plt.subplot(131)
    plt.bar(grouped['base'], grouped['total_count'], color=['red', 'blue', 'green', 'purple'])
    plt.title('Total Base Counts')
    plt.ylabel('Count')
    plt.xlabel('Base')
    
    plt.subplot(132)
    plt.bar(grouped['base'], grouped['mean_base_quality'], color=['red', 'blue', 'green', 'purple'])
    plt.title('Mean Base Quality')
    plt.ylabel('Quality Score')
    plt.xlabel('Base')
    plt.ylim(0, max(grouped['mean_base_quality']) * 1.2)
    
    plt.subplot(133)
    positions = np.arange(len(grouped))
    plt.bar(positions, grouped['pos_strand_count'], width=0.6, label='Positive Strand')
    plt.bar(positions, grouped['neg_strand_count'], width=0.6, bottom=grouped['pos_strand_count'], label='Negative Strand')
    plt.title('Strand Distribution')
    plt.ylabel('Count')
    plt.xlabel('Base')
    plt.xticks(positions, grouped['base'])
    plt.legend()
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"{output_dir}/{sample}_analysis.png", dpi=150)
    plt.close()

for sample, file_path in sample_to_file.items():
    process_sample(sample, file_path)

def create_summary_comparison(sample_to_file):
    all_data = []
    for sample, file_path in sample_to_file.items():
        if not os.path.exists(file_path):
            continue
        df = pd.read_csv(file_path, sep='\t')
        df = df[df['base'].isin(['A', 'T', 'C', 'G'])]
        grouped = df.groupby('base').agg({
            'total_count': 'sum',
            'mean_base_quality': lambda x: np.average(x, weights=df.loc[x.index, 'total_count']),
            'pos_strand_count': 'sum',
            'neg_strand_count': 'sum'
        }).reset_index()
        grouped['sample'] = sample
        all_data.append(grouped)

    if not all_data:
        return
    combined = pd.concat(all_data)
    samples = list(sample_to_file.keys())
    valid_bases = ['A', 'T', 'C', 'G']

    plt.figure(figsize=(12, 6))
    for i, base in enumerate(valid_bases):
        base_data = combined[combined['base'] == base]
        positions = np.arange(len(samples)) + i * 0.2
        plt.bar(positions, base_data['total_count'], width=0.2, label=base)
    plt.title('Total Base Counts Comparison')
    plt.ylabel('Count')
    plt.xticks(np.arange(len(samples)) + 0.3, samples, rotation=45)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{output_dir}/all_samples_base_counts.png", dpi=150)
    plt.close()

    plt.figure(figsize=(12, 6))
    for base in valid_bases:
        base_data = combined[combined['base'] == base]
        plt.plot(base_data['sample'], base_data['mean_base_quality'], 'o-', label=base)
    plt.title('Mean Base Quality Comparison')
    plt.ylabel('Quality Score')
    plt.xlabel('Sample')
    plt.xticks(rotation=45)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/all_samples_quality.png", dpi=150)
    plt.close()

create_summary_comparison(sample_to_file)

print(f"All plots saved to: {output_dir}/")

