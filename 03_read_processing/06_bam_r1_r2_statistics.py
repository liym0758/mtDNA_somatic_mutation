import os
import pysam
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import numpy as np

samples = [
    "mouse", 
    "human", 
    "rat",
    "chicken", 
    "cobra", 
    "sturgeon", 
    "stingray"
]

base_dir = "/PATH/TO/pbmc.mt/SNP/"
output_dir = os.path.join(base_dir,  "pre.process/clean.bam/7species_flag_analysis_results")
os.makedirs(output_dir, exist_ok=True)
print("Output dir:", output_dir)

results = []
pair_analysis_results = []

for sample in samples:
    print("==================Processing: ",sample,"=====================")
    bam_path = os.path.join(base_dir, sample, "mgatk.CleanBam/bam/possorted_bam.chrM.MAPQ40.rmdup.filtered.realign.bam")
    
    print("BAM path:", bam_path)

    if not os.path.exists(bam_path):
        print(f"Warning: BAM file not found for sample {sample}. Skipping...")
        continue
    
    count_99 = 0
    count_147 = 0
    count_83 = 0
    count_163 = 0
    
    total_read1 = 0
    total_read2 = 0
    
    total_positive = 0
    total_negative = 0
    
    read_dict = defaultdict(lambda: {"read1": None, "read2": None})
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_read1:
                read_dict[read.query_name]["read1"] = read
            elif read.is_read2:
                read_dict[read.query_name]["read2"] = read
    
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            flag = read.flag
            
            if flag == 99:
                count_99 += 1
                total_positive += 1
            elif flag == 147:
                count_147 += 1
                total_negative += 1
            elif flag == 83:
                count_83 += 1
                total_negative += 1
            elif flag == 163:
                count_163 += 1
                total_positive += 1
            
            if read.is_read1:
                total_read1 += 1
            elif read.is_read2:
                total_read2 += 1
    
    pair_status = defaultdict(int)
    for qname, reads in read_dict.items():
        has_read1 = reads["read1"] is not None
        has_read2 = reads["read2"] is not None
        
        if has_read1 and has_read2:
            pair_status["complete_pairs"] += 1
        elif has_read1:
            pair_status["read1_only"] += 1
        elif has_read2:
            pair_status["read2_only"] += 1
        else:
            pair_status["no_reads"] += 1
    
    total_reads = total_read1 + total_read2
    read_imbalance = abs(total_read1 - total_read2)
    read_imbalance_pct = (read_imbalance / total_reads) * 100
    
    pct_positive = (total_positive / total_reads) * 100
    pct_negative = (total_negative / total_reads) * 100
    
    results.append({
        "Sample": sample,
        "Flag_99": count_99,
        "Flag_147": count_147,
        "Flag_83": count_83,
        "Flag_163": count_163,
        "Total_Reads": total_reads,
        "Total_Read1": total_read1,
        "Total_Read2": total_read2,
        "Total_Positive": total_positive,
        "Total_Negative": total_negative,
        "Pct_Positive": pct_positive,
        "Pct_Negative": pct_negative,
        "Read_Imbalance": read_imbalance,
        "Read_Imbalance_Pct": read_imbalance_pct
    })
    
    pair_analysis_results.append({
        "Sample": sample,
        "Complete_Pairs": pair_status["complete_pairs"],
        "Read1_Only": pair_status["read1_only"],
        "Read2_Only": pair_status["read2_only"],
        "Total_Fragments": pair_status["complete_pairs"] + 
                           pair_status["read1_only"] + 
                           pair_status["read2_only"]
    })
    
    print(f"Processed {sample}:")
    print(f"  Total reads: {total_reads:,}")
    print(f"  Read1: {total_read1:,} ({total_read1/total_reads*100:.2f}%)")
    print(f"  Read2: {total_read2:,} ({total_read2/total_reads*100:.2f}%)")
    print(f"  Positive strand: {total_positive:,} ({pct_positive:.2f}%)")
    print(f"  Negative strand: {total_negative:,} ({pct_negative:.2f}%)")
    print(f"  Imbalance: {read_imbalance:,} reads ({read_imbalance_pct:.4f}%)")
    print(f"  Complete pairs: {pair_status['complete_pairs']:,}")
    print(f"  Read1 without mate: {pair_status['read1_only']:,}")
    print(f"  Read2 without mate: {pair_status['read2_only']:,}")

df = pd.DataFrame(results)
pair_df = pd.DataFrame(pair_analysis_results)

output_csv = os.path.join(output_dir, "realigned_bam_flag_detailed_stats.csv")
df.to_csv(output_csv, index=False)
print(f"\nDetailed statistics saved to: {output_csv}")

pair_csv = os.path.join(output_dir, "realigned_pair_analysis_stats.csv")
pair_df.to_csv(pair_csv, index=False)
print(f"Pair statistics saved to: {pair_csv}")

plt.rcParams.update({'font.size': 12})
sns.set_style("whitegrid")
plt.figure(figsize=(18, 12))

ax1 = plt.subplot(2, 2, 1)
x = np.arange(len(samples))
width = 0.35

rects1 = ax1.bar(x - width/2, df["Total_Read1"], width, label='Read1', color='#1f77b4')
rects2 = ax1.bar(x + width/2, df["Total_Read2"], width, label='Read2', color='#ff7f0e')

ax1.set_ylabel('Read Count')
ax1.set_title('Read1 vs Read2 Count by Sample')
ax1.set_xticks(x)
ax1.set_xticklabels(samples, rotation=45)
ax1.legend()

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        ax1.annotate(f'{height/1e6:.2f}M',
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)

ax2 = plt.subplot(2, 2, 2)
width = 0.35

rects3 = ax2.bar(x - width/2, df["Total_Positive"], width, label='Positive Strand', color='#2ca02c')
rects4 = ax2.bar(x + width/2, df["Total_Negative"], width, label='Negative Strand', color='#d62728')

ax2.set_ylabel('Read Count')
ax2.set_title('Positive vs Negative Strand Distribution')
ax2.set_xticks(x)
ax2.set_xticklabels(samples, rotation=45)
ax2.legend()

for i, rect in enumerate(rects3):
    height = rect.get_height()
    ax2.text(rect.get_x() + rect.get_width()/2, height + height*0.01, 
             f"{df.loc[i, 'Pct_Positive']:.2f}%", ha='center', fontsize=9)

for i, rect in enumerate(rects4):
    height = rect.get_height()
    ax2.text(rect.get_x() + rect.get_width()/2, height + height*0.01, 
             f"{df.loc[i, 'Pct_Negative']:.2f}%", ha='center', fontsize=9)

ax3 = plt.subplot(2, 2, 3)
complete = pair_df["Complete_Pairs"]
read1_only = pair_df["Read1_Only"]
read2_only = pair_df["Read2_Only"]

p1 = ax3.bar(samples, complete, label='Complete Pairs', color='#1f77b4')
p2 = ax3.bar(samples, read1_only, bottom=complete, label='Read1 Only', color='#ff7f0e')
p3 = ax3.bar(samples, read2_only, bottom=complete+read1_only, label='Read2 Only', color='#2ca02c')

ax3.set_ylabel('Fragment Count')
ax3.set_title('Fragment Pairing Status by Sample')
ax3.legend()

for i, (c, r1, r2) in enumerate(zip(complete, read1_only, read2_only)):
    total = c + r1 + r2
    ax3.text(i, total + total*0.01, f"{total:,}", ha='center', va='bottom', fontsize=9)
    ax3.text(i, c/2, f"{c:,}", ha='center', va='center', color='white', fontsize=8)
    ax3.text(i, c + r1/2, f"{r1:,}", ha='center', va='center', fontsize=8)
    ax3.text(i, c + r1 + r2/2, f"{r2:,}", ha='center', va='center', fontsize=8)

ax4 = plt.subplot(2, 2, 4)
flag_types = ['99', '147', '83', '163']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

for i, sample in enumerate(samples):
    flags = [
        df.loc[df['Sample'] == sample, 'Flag_99'].values[0],
        df.loc[df['Sample'] == sample, 'Flag_147'].values[0],
        df.loc[df['Sample'] == sample, 'Flag_83'].values[0],
        df.loc[df['Sample'] == sample, 'Flag_163'].values[0]
    ]
    
    for j, (count, color) in enumerate(zip(flags, colors)):
        ax4.bar(i + j*0.2, count, width=0.2, color=color, alpha=0.8)
    
    total = sum(flags)
    ax4.text(i + 0.3, max(flags) + max(flags)*0.1, f"{total/1e6:.2f}M", 
             ha='center', fontsize=9)

ax4.set_xticks(np.arange(len(samples)))
ax4.set_xticklabels(samples, rotation=45)
ax4.set_ylabel('Read Count')
ax4.set_title('Flag Type Distribution by Sample')

legend_labels = [plt.Rectangle((0,0),1,1, color=color) for color in colors]
ax4.legend(legend_labels, flag_types, title='Flag Type')

plt.tight_layout()

output_png = os.path.join(output_dir, "realigned_bam_comprehensive_analysis.png")
plt.savefig(output_png, dpi=300)
print(f"Combined plot saved to: {output_png}")

plt.figure(figsize=(14, 8))
for i, sample in enumerate(samples):
    flags = [
        df.loc[df['Sample'] == sample, 'Flag_99'].values[0],
        df.loc[df['Sample'] == sample, 'Flag_147'].values[0],
        df.loc[df['Sample'] == sample, 'Flag_83'].values[0],
        df.loc[df['Sample'] == sample, 'Flag_163'].values[0]
    ]
    
    for j, (count, color) in enumerate(zip(flags, colors)):
        plt.bar(i + j*0.2, count, width=0.2, color=color, alpha=0.8)
    
    total = sum(flags)
    plt.text(i + 0.3, max(flags) + max(flags)*0.1, f"{total/1e6:.2f}M", 
             ha='center', fontsize=10)

plt.xticks(np.arange(len(samples)), samples, rotation=45)
plt.ylabel('Read Count')
plt.title('Detailed Flag Type Distribution')
legend_labels = [plt.Rectangle((0,0),1,1, color=color) for color in colors]
plt.legend(legend_labels, flag_types, title='Flag Type')
plt.tight_layout()

flag_png = os.path.join(output_dir, "realigned_bam_flag_distribution_detailed.png")
plt.savefig(flag_png, dpi=300)
print(f"Detailed flag distribution saved to: {flag_png}")

print("\nAnalysis summary:")
print(df[["Sample", "Total_Read1", "Total_Read2", "Total_Positive", "Total_Negative", "Read_Imbalance"]])
print("\nPair status summary:")
print(pair_df[["Sample", "Complete_Pairs", "Read1_Only", "Read2_Only"]])
