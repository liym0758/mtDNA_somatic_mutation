function revcomp(base) {
    switch(base) {
        case "A": return "T";
        case "T": return "A";
        case "G": return "C";
        case "C": return "G";
        default: return "N";
    }
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

function ord(c) {
    return sprintf("%d", c)
}

BEGIN {
    OFS = "\t";
    print "chr", "position", "base", "reads#", "mean_base_quality", "strand"
}

{
    chr = $1
    pos = $2
    ref = toupper($3)
    depth = $4
    bases = $5
    quals = $6
    
    process_bases(bases, quals, ref)
    
    for (base in counts) {
        avg_qual = sum_qual[base] / counts[base]
        strand_info = strand_pos[base] "+" strand_neg[base] "-"
        print chr, pos, base, counts[base], sprintf("%.2f", avg_qual), strand_info
    }
}
