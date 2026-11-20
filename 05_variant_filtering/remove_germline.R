setwd("/PATH/TO/pbmc.mt/SNP/")
library(Biostrings)
library(data.table)

Remove_Germline <- function(species, variant_file, ref_fa, output_dir = ".", vmr_cutoff = 0.01) {
  message("Processing species: ", species)
  
  ref <- readDNAStringSet(ref_fa)
  ref_seq <- ref[[1]]
  ref_name <- names(ref)
  
  variants <- fread(variant_file, sep = "\t", header = TRUE)
  # position, nucleotide, variant, vmr, ...
  if (!all(c("position", "nucleotide", "variant", "vmr") %in% colnames(variants))) {
    stop("variant_file format error: must contain 'position', 'nucleotide', 'variant', 'vmr' columns")
  }
  
  germline <- variants %>% dplyr::filter(vmr <= vmr_cutoff & n_cells_conf_detected > 3)
  if (nrow(germline) == 0) {
    warning("No germline mutations found for ", species, ". Reference genome unchanged.")
    output_fa <- file.path(output_dir, paste0(species, ".mt_genome.germline_corrected.fa"))
    writeXStringSet(ref, filepath = output_fa)
    return(invisible(NULL))
  }
  
  germline[, c("Ref", "Alt") := tstrsplit(nucleotide, ">")]
  
  seq_chars <- strsplit(as.character(ref_seq), "")[[1]]
  warnings <- character()
  
  for (i in seq_len(nrow(germline))) {
    pos <- germline$position[i]
    ref_base <- germline$Ref[i]
    alt_base <- germline$Alt[i]
    
    ref_base_genome <- seq_chars[pos]
    if (ref_base_genome != ref_base) {
      warnings <- c(warnings, paste0("Warning [", species, "]: Position ", pos, 
                                     " expected ", ref_base, " but found ", ref_base_genome))
    } else {
      seq_chars[pos] <- alt_base
    }
  }
  
  new_seq <- DNAString(paste(seq_chars, collapse = ""))
  new_ref <- DNAStringSet(new_seq)
  names(new_ref) <- ref_name
  output_fa <- file.path(output_dir, paste0(species, ".germline.rm.MT.genome.fa"))
  writeXStringSet(new_ref, filepath = output_fa)
  
  ref_seq_chars <- strsplit(as.character(ref_seq), "")[[1]]
  
  report <- germline[, .(position, Ref, Alt)]
  report[, ref_genome := ref_seq_chars[position]]
  
  report[, status := ifelse(ref_genome == Ref, "OK", "WARNING_CONFLICT")]
  
  report_file <- file.path(output_dir, paste0(species, ".germline_remove_report.tsv"))
  fwrite(report, report_file, sep = "\t")
  message("Report written to: ", report_file)
  
  
  
  if (any(duplicated(germline$position))) {
    conflict_file <- file.path(output_dir, paste0(species, ".conflict_positions.tsv"))
    fwrite(germline[duplicated(germline$position) | duplicated(germline$position, fromLast = TRUE)], conflict_file, sep = "\t")
    warning("Conflicts detected: multiple mutations at the same position. See ", conflict_file)
  }
  
  if (length(warnings) > 0) {
    cat(paste(warnings, collapse = "\n"), "\n")
  } else {
    message("All germline mutations successfully applied for ", species)
  }
  
  return(invisible(output_fa))
}


mgatk.out.path <- c(
  #"human" = "/PATH/TO/pbmc.mt/SNP/human/mgatk/final/",
  #"mouse" = "/PATH/TO/pbmc.mt/SNP/mouse/mgatk/final/",
  #"rat" = "/PATH/TO/pbmc.mt/SNP/rat/mgatk/final/",
  "chicken" = "/PATH/TO/pbmc.mt/SNP/chicken/mgatk/final/",
  #"cobra" = "/PATH/TO/pbmc.mt/SNP/cobra/mgatk/final/",
  "sturgeon" = "/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk/final/"#,
  #"stingray" = "/PATH/TO/pbmc.mt/SNP/stingray/mgatk/final/"
)

variant_files <- paste0(mgatk.out.path, names(mgatk.out.path), ".mgatk.variant_stats.tsv.gz")
ref_list <- c(
  #"human" = "/PATH/TO/pbmc.mt/ref/human/GRCh38.MT.fasta",
  #"mouse" = "/PATH/TO/pbmc.mt/ref/mouse/new/mm39_mito.fa",
  #"rat" = "/PATH/TO/pbmc.mt/ref/rat/Rat_ref/mRatBN7.2_Ensembl/mRatBN7_mt.fa",
  "chicken" = "/PATH/TO/pbmc.mt/ref/chicken/ensembl/Gallus_gallus.bGalGal1_mt.fa",
  #"cobra" = "/PATH/TO/pbmc.mt/ref/cobra/Naja_atra/Naja_atra.mitochondrion.complete.genome.fa",
  "sturgeon" = "/PATH/TO/pbmc.mt/ref/XunFish/NCBI.mt/Acipenser_sinensis.mt.complete.genome.fa"#,
  #"stingray" = "/PATH/TO/pbmc.mt/ref/HongFish/NCBI.mt/Potamotrygon_leopoldi.mt.genome.fa"
)
output_dir <- paste0(dirname(mgatk.out.path),"/downstream/rm.germline/")

if (!dir.exists(output_dir)) dir.create(output_dir)

for (i in seq_along(species_list)) {
  Remove_Germline(
    species = names(mgatk.out.path)[i],
    variant_file = variant_files[i],
    ref_fa = ref_list[i],
    output_dir = output_dir[i]
  )
}

paste0(output_dir,names(mgatk.out.path), ".germline.rm.MT.genome.fa")
