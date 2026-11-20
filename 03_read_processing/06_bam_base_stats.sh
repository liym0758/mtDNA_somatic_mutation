#!/bin/bash

# Initialize variables
BAM_FILE=""
MPILEUP_FILE=""
OUTPUT_FILE=""
HELP_MESSAGE="Usage: $0 -b <BAM path> | -m <mpileup path> [-o <output path>]
Options:
  -b    Path to input BAM file (requires samtools)
  -m    Path to existing mpileup file (avoids recalculating)
  -o    Path to output file (optional, default: print to screen)
  
Note: -b and -m are mutually exclusive (use only one)"

while getopts ":b:m:o:h" opt; do
  case $opt in
    b) BAM_FILE="$OPTARG" ;;
    m) MPILEUP_FILE="$OPTARG" ;;
    o) OUTPUT_FILE="$OPTARG" ;;
    h) echo "$HELP_MESSAGE"; exit 0 ;;
    \?) echo "Invalid Argument: -$OPTARG" >&2; echo "$HELP_MESSAGE"; exit 1 ;;
    :) echo "Option -$OPTARG requires an argument." >&2; echo "$HELP_MESSAGE"; exit 1 ;;
  esac
done

if [ -z "$BAM_FILE" ] && [ -z "$MPILEUP_FILE" ]; then
  echo "ERROR: Must specify either -b or -m parameter"
  echo "$HELP_MESSAGE"
  exit 1
fi

if [ -n "$BAM_FILE" ] && [ -n "$MPILEUP_FILE" ]; then
  echo "ERROR: Cannot use both -b and -m parameters simultaneously"
  echo "$HELP_MESSAGE"
  exit 1
fi

if [ -n "$BAM_FILE" ]; then
  INPUT_TYPE="BAM"
  INPUT_SOURCE="$BAM_FILE"
  [ ! -f "$BAM_FILE" ] && echo "ERROR: BAM file not found: $BAM_FILE" && exit 1
  else
  INPUT_TYPE="mpileup"
  INPUT_SOURCE="$MPILEUP_FILE"
  [ ! -f "$MPILEUP_FILE" ] && echo "ERROR: mpileup file not found: $MPILEUP_FILE" && exit 1
fi

SCRIPT_DIR=$(dirname "$0")
TEMP_AWK="${SCRIPT_DIR}/strand_count_$$.awk"

cat > "$TEMP_AWK" << 'EOF'
function revcomp(base) {
    base = toupper(base)
    if (base == "A") return "T"
    else if (base == "T") return "A"
    else if (base == "G") return "C"
    else if (base == "C") return "G"
    else return "N"
}

function ord(c) {
    return index("\0\1\2\3\4\5\6\7\10\11\12\13\14\15\16\17\20\21\22\23\24\25\26\27\30\31 !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", c) - 1
}


function process_bases(bases, quals, ref,   c, i, qual_idx, len_str, len, base_char, strand, qual_val) {
    delete counts
    delete sum_qual
    delete strand_pos
    delete strand_neg

    i = 1
    qual_idx = 1
    while (i <= length(bases)) {
        c = substr(bases, i, 1)

        if (c == "^") { i += 2; continue }
        if (c == "$") { i++; continue }
        if (c == "*") { i++; qual_idx++; continue }

        if (c == "+" || c == "-") {
            i++
            len_str = ""
            while (substr(bases, i, 1) ~ /[0-9]/) {
                len_str = len_str substr(bases, i, 1)
                i++
            }
            len = len_str ? int(len_str) : 0
            i += len
            continue
        }

        if (qual_idx > length(quals)) {
            qual_val = 0
        } else {
            q_char = substr(quals, qual_idx, 1)
            qual_val = q_char != "" ? ord(q_char) - 33 : 0
        }
        qual_idx++

        strand = "+"
        base_char = toupper(c)

        if (c ~ /[.,]/) {
            base_char = ref
            strand = (c == ".") ? "+" : "-"
        }
        else if (c ~ /[a-z]/) {
            base_char = revcomp(toupper(c))
            strand = "-"
        }
        else if (c ~ /[A-Z]/) {
            strand = "+"
        }

        if (base_char ~ /[ATGC]/) {
            counts[base_char]++
            sum_qual[base_char] += qual_val
            if (strand == "+") strand_pos[base_char]++
            else strand_neg[base_char]++
        }
        i++
    }
}

BEGIN {
    OFS = "\t"
    print "chr", "position", "base", "total_count", "mean_base_quality", "pos_strand_count", "neg_strand_count"
}

{
    if (NF < 6) next
    chr = $1; pos = $2; ref = toupper($3); bases = $5; quals = $6
    process_bases(bases, quals, ref)
    for (base in counts) {
        avg_qual = sum_qual[base] / counts[base]
        pos_count = (base in strand_pos ? strand_pos[base] : 0)
        neg_count = (base in strand_neg ? strand_neg[base] : 0)
	print chr, pos, base, counts[base], sprintf("%.2f", avg_qual), pos_count, neg_count
    }
}
EOF

echo "================ START PROCESSING $INPUT_TYPE FILE ================="
echo "Input source: $INPUT_SOURCE"

if [ "$INPUT_TYPE" = "BAM" ]; then
  PROCESSING_PIPELINE="samtools mpileup -Q 0 -A $BAM_FILE"
else
  PROCESSING_PIPELINE="cat $MPILEUP_FILE"
fi

if [ -z "$OUTPUT_FILE" ]; then
  echo "Processing pipeline: $PROCESSING_PIPELINE | awk | sort | column"
  eval "$PROCESSING_PIPELINE" | awk -f "$TEMP_AWK" | sort -k1,1 -k2,2n -k3,3 | column -t -s $'\t'
  echo "================ ANALYSIS COMPLETED ================="
else
  echo "Writing output to: $OUTPUT_FILE"
  echo "Processing pipeline: $PROCESSING_PIPELINE | awk | sort > '$OUTPUT_FILE'"
  eval "$PROCESSING_PIPELINE" | awk -f "$TEMP_AWK" | sort -k1,1 -k2,2n -k3,3 > "$OUTPUT_FILE"
  line_count=$(wc -l < "$OUTPUT_FILE")
  echo "================ ANALYSIS COMPLETED ================="
  echo "Output file: $OUTPUT_FILE"
  echo "Total positions analyzed: $line_count"
fi

rm -f "$TEMP_AWK"
exit 0

