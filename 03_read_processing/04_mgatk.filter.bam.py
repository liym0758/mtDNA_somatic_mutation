#!/usr/bin/env python3

import sys
import pysam

###############################################
# Filter BAM by:
# - proper pairs (optional)
# - chromosome name
# - NM <= NMmax
# - NH == NHvalue
# - MAPQ >= mapq_min
###############################################
# USAGE:
# python 04_mgatk.filter.bam.py input.bam output.log chrM True 1 2 > filtered.bam
###############################################

bamfile = sys.argv[1]
logfile = sys.argv[2]
mtchr = sys.argv[3]
proper_pair = sys.argv[4]
NHmax = int(sys.argv[5])
NMmax = int(sys.argv[6])

bam = pysam.AlignmentFile(bamfile, "rb")

new_header = str(bam.header).replace("@HD\tSO:coordinate", "@HD\tVN:1.5\tSO:coordinate")
out = pysam.AlignmentFile("-", "wb", text=new_header)

keepCount = 0
filtCount = 0

def filterReadTags(tags):
    for tag, value in tags:
        if tag == "NH" and int(value) > NHmax:
            return False
        if (tag == "NM" or tag == "nM") and int(value) > NMmax:
            return False
    return True

def pairingCheck(read):
    if proper_pair != "True":
        return True
    else:
        return read.is_proper_pair

for read in bam:
    if filterReadTags(read.tags) and read.reference_name == mtchr and pairingCheck(read):
        out.write(read)
        keepCount += 1
    else:
        filtCount += 1

bam.close()

with open(logfile, "w") as logfile:
    logfile.write(f"Kept {keepCount}\nRemoved {filtCount}\n")
