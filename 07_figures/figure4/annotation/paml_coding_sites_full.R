
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)
library(purrr)
library(IRanges)



read_genepred <- function(path) {
  message("Reading GenePred file: ", path)
  
  gp <- fread(path, header = FALSE, sep = "\t")
  
  message("Detected column count: ", ncol(gp))
  
  if (ncol(gp) >= 12) {
    colnames(gp)[1:12] <- c("name", "chrom", "strand", "txStart", "txEnd", 
                            "cdsStart", "cdsEnd", "exonCount", "exonStarts", 
                            "exonEnds", "score", "name2")
    
    if (ncol(gp) >= 15) {
      colnames(gp)[13:15] <- c("cdsStartStat", "cdsEndStat", "exonFrames")
    }
    if (ncol(gp) >= 16) {
      colnames(gp)[1] <- "bin"  # First column may be a bin column
      colnames(gp)[2:16] <- c("name", "chrom", "strand", "txStart", "txEnd", 
                              "cdsStart", "cdsEnd", "exonCount", "exonStarts", 
                              "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")
    }
  }
  
  message("Column names used: ", paste(colnames(gp), collapse = ", "))
  
  gp <- gp %>%
    filter(toupper(chrom) %in% c("MT", "CHRM", "CHRMT", "M") | 
             grepl("MT|MITO", chrom, ignore.case = TRUE))
  
  numeric_cols <- c("txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "score")
  for (col in numeric_cols) {
    if (col %in% names(gp)) {
      gp[[col]] <- as.numeric(gp[[col]])
    }
  }
  
  gp <- gp %>%
    mutate(
      gene_name = ifelse(!is.na(name2) & name2 != "" & name2 != "none", name2, name),
      gene_length = txEnd - txStart + 1
    )
  
  gp$gene_type <- character(nrow(gp))
  for (i in 1:nrow(gp)) {
    gp$gene_type[i] <- classify_gene_type_genepred(gp$gene_name[i], gp$gene_length[i])
  }
  
  gp_genes <- gp
  unique_genes <- unique(gp_genes$gene_name)
  result_list <- list()
  
  for (gene in unique_genes) {
    gene_subset <- gp_genes[gp_genes$gene_name == gene, ]
    max_length_idx <- which.max(gene_subset$gene_length)
    if (length(max_length_idx) > 0) {
      result_list[[gene]] <- gene_subset[max_length_idx[1], ]
    }
  }
  
  gp_genes <- do.call(rbind, result_list)
  rownames(gp_genes) <- NULL
  
  gp_genes <- gp_genes[, c("gene_name", "chrom", "strand", "txStart", "txEnd", 
                           "exonStarts", "exonEnds", "exonCount", "gene_length", "gene_type")]
  
  message("Identified in GenePred file ", nrow(gp_genes), " mitochondrial genes")
  message("Gene-type summary:")
  print(table(gp_genes$gene_type))
  
  if (nrow(gp_genes) > 0) {
    message("First five gene examples:")
    print(head(gp_genes[, c("gene_name", "gene_type", "gene_length")], 5))
  }
  
  return(gp_genes)
}

classify_gene_type_genepred <- function(gene_name, gene_length = NA) {
  gene_name <- as.character(gene_name)
  if (is.na(gene_name)) return("ncRNA")
  
  gene_upper <- toupper(gene_name)
  
  if (grepl("^MT-T|^T", gene_upper)) return("tRNA")
  
  if (grepl("^MT-R|^R|16S|12S", gene_upper)) return("rRNA")
  
  if (grepl("^MT-", gene_upper)) {
    if (grepl("^MT-ND|^MT-CO|^MT-ATP|^MT-CY", gene_upper)) {
      return("CDS")
    }
  }
  
  cds_pattern <- "ND[1-6]|CO[I-III1-3]|ATP[6-8]|CYB|CYTB|NAD[1-6]|COX[1-3]"
  if (grepl(cds_pattern, gene_upper)) return("CDS")
  
  if (!is.na(gene_length)) {
    if (gene_length >= 50 && gene_length <= 100) return("tRNA")      # Typical tRNA length
    if (gene_length >= 500 && gene_length <= 3000) return("rRNA")    # Typical rRNA length
    if (gene_length >= 100 && gene_length <= 2000) return("CDS")     # Typical protein-coding gene length
  }
  
  "ncRNA"
}


extract_cds_dna_simple <- function(mt_dna, gp_row) {
  tryCatch({
    s <- as.numeric(gp_row$txStart) + 1  # Convert to 1-based
    e <- as.numeric(gp_row$txEnd)
    
    if (s < 1 || e > length(mt_dna) || s > e) {
      stop(sprintf("Invalid coordinates: start=%d, end=%d, sequence length=%d", s, e, length(mt_dna)))
    }
    
    cds <- subseq(mt_dna, start = s, end = e)
    
    if (length(cds) == 0) {
      stop("Extracted CDS length is zero")
    }
    
    if (as.character(gp_row$strand) == "-") {
      cds <- reverseComplement(cds)
    }
    
    L <- length(cds)
    if (L < 3) {
      stop("CDS length is under three bases: ", gp_row$gene_name, " (length=", L, ")")
    }
    
    if (L %% 3 != 0) {
      drop <- L %% 3
      cds <- cds[1:(L - drop)]
      message("Truncating: ", gp_row$gene_name, " | original length=", L, " | new length=", length(cds))
    }
    
    seq_vec <- strsplit(as.character(cds), "")[[1]]
    valid_bases <- c("A", "C", "G", "T")
    seq_vec[!seq_vec %in% valid_bases] <- "N"
    cds_clean <- DNAString(paste(seq_vec, collapse = ""))
    
    cds_clean
  }, error = function(e) {
    stop("CDS extraction failed (", gp_row$gene_name, "): ", conditionMessage(e))
  })
}

write_codeml_files <- function(outdir, basename, dna_cds, icode = 1) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  seqfile <- file.path(outdir, paste0(basename, "_same_codon.txt"))
  treefile <- file.path(outdir, paste0(basename, "_same.tree"))
  ctlfile <- file.path(outdir, paste0(basename, "_same.ctl"))
  
  L <- length(dna_cds)
  if (L %% 3 != 0) {
    drop <- L %% 3
    dna_cds <- dna_cds[1:(L - drop)]
    message("Truncating: ", basename, " | original length=", L, " | new length=", length(dna_cds))
  }
  
  cds_string <- as.character(dna_cds)
  nsite <- nchar(cds_string)
  
  seq_content <- paste0("2 ", nsite, "\n",
                        "sp1       ", cds_string, "\n",
                        "sp2       ", cds_string)
  
  writeLines(seq_content, con = seqfile)
  
  tree_content <- "(sp1, sp2);"
  writeLines(tree_content, con = treefile)
  
  ctl_lines <- c(
    paste0("seqfile = ", basename, "_same_codon.txt"),
    paste0("treefile = ", basename, "_same.tree"),
    "outfile = mlc",
    "",
    "noisy = 9   * 0,1,2,3,9: how much rubbish on the screen",
    "verbose = 1   * 1: detailed output, 0: concise output", 
    "runmode = 0   * 0: user tree;  1: semi-automatic;  2: automatic",
    "            * 3: StepwiseAddition; (4,5):PerturbationNNI", 
    "",
    "seqtype = 1   * 1:codons; 2:AAs; 3:codons-->AAs",
    "CodonFreq = 2   * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table",
    "clock = 0   * 0: no clock, unrooted tree, 1: clock, rooted tree", 
    "model = 0",
    "            * models for codons:",
    "                * 0:one, 1:b, 2:2 or more dN/dS ratios for branches",
    "",
    "NSsites = 0   * dN/dS among sites. 0:no variation, 1:neutral, 2:positive", 
    paste0("icode = ", icode, "   * 0:standard genetic code; 1:mammalian mt; 2-10:see below"),
    "",
    "fix_kappa = 0   * 1: kappa fixed, 0: kappa to be estimated",
    "kappa = 2   * initial or fixed kappa", 
    "fix_omega = 0   * 1: omega or omega_1 fixed, 0: estimate", 
    "omega = 2   * initial or fixed omega, for codons or codon-transltd AAs",
    "",
    "fix_alpha = 1   * 0: estimate gamma shape parameter; 1: fix it at alpha",
    "alpha = .0  * initial or fixed alpha, 0:infinity (constant rate)",
    "Malpha = 0   * different alphas for genes",
    "ncatG = 4   * # of categories in the dG or AdG models of rates",
    "",
    "getSE = 0   * 0: don't want them, 1: want S.E.s of estimates",
    "RateAncestor = 0   * (1/0): rates (alpha>0) or ancestral states (alpha=0)",
    "method = 0   * 0: simultaneous; 1: one branch at a time",
    "fix_blength = 0  * 0: ignore, -1: random, 1: initial, 2: fixed, 3: proportional"
  )
  
  writeLines(ctl_lines, con = ctlfile)
  
  list(seqfile = seqfile, treefile = treefile, ctlfile = ctlfile)
}

run_codeml_and_parse_SN <- function(workdir, ctlfile, codeml = "codeml") {
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  setwd(workdir)
  
  if (file.exists("mlc")) {
    file.remove("mlc")
  }
  
  ctl_basename <- basename(ctlfile)
  result <- system2(codeml, args = ctl_basename, stdout = TRUE, stderr = TRUE)
  
  if (!file.exists("mlc")) {
    warning("codeml did not produce an mlc file")
    writeLines(result, "codeml_output.log")
    return(c(S = NA, N = NA))
  }
  
  mlc_content <- readLines("mlc", warn = FALSE)
  
  if (length(mlc_content) == 0) {
    warning("mlc file is empty")
    return(c(S = NA, N = NA))
  }
  
  S <- NA_real_
  N <- NA_real_
  
  for (i in rev(seq_along(mlc_content))) {
    line <- mlc_content[i]
    
    if (grepl("^\\s*3\\.\\.", line)) {
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      
      if (length(parts) >= 6) {
        N <- as.numeric(parts[3])
        S <- as.numeric(parts[4])
        
        if (!is.na(S) && !is.na(N) && S > 0 && N > 0) {
          message("Parsed successfully: S = ", round(S, 1), ", N = ", round(N, 1))
          break
        } else {
          N <- NA_real_
          S <- NA_real_
        }
      }
    }
  }
  
  if (is.na(S) || is.na(N)) {
    warning("Cannot parse S and N values")
    tail_lines <- tail(mlc_content, 20)
    writeLines(tail_lines, "debug_mlc_tail.txt")
    message("Last 20 lines of mlc saved to debug_mlc_tail.txt")
  }
  
  return(c(S = S, N = N))
}


compute_background_with_paml <- function(fasta_path,
                                         genepred_path,
                                         out_root,
                                         codeml = "codeml",
                                         icode = 1) {
  
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  
  log_file <- file.path(out_root, "analysis.log")
  log_con <- file(log_file, "w")
  
  log_message <- function(...) {
    msg <- paste(...)
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    full_msg <- paste0("[", timestamp, "] ", msg)
    message(full_msg)
    writeLines(full_msg, log_con)
    flush(log_con)
  }
  
  tryCatch({
    log_message("Starting species: ", basename(out_root))
    
    fa <- readDNAStringSet(fasta_path)
    mt_idx <- which(grepl("MT|mito|chrM", names(fa), ignore.case = TRUE))
    if (length(mt_idx) == 0) stop("Mitochondrial sequence not found")
    mt <- fa[[mt_idx[1]]]
    mt_len <- length(mt)
    log_message("Mitochondrial sequence length: ", mt_len)
    
    log_message("Reading GenePred file...")
    gp <- read_genepred(genepred_path)
    
    if (nrow(gp) == 0) {
      stop("No mitochondrial genes found in GenePred file")
    }
    
    gene_intervals <- IRanges::IRanges(start = gp$txStart + 1, end = gp$txEnd)
    gene_union <- IRanges::reduce(gene_intervals)
    total_gene_len <- sum(IRanges::width(gene_union))
    
    
    gp_sorted <- gp[order(gp$txStart), ]
    if (nrow(gp_sorted) >= 2) {
      starts_sorted <- gp_sorted$txStart + 1
      ends_sorted   <- gp_sorted$txEnd
      gap_vec <- pmax(0, starts_sorted[-1] - ends_sorted[-length(ends_sorted)] - 1)
      len_intergenic <- sum(gap_vec)
    } else {
      len_intergenic <- 0
    }
    
    len_dloop <- max(0, mt_len - max(gp_sorted$txEnd) + min(gp_sorted$txStart))
    
    len_tRNA <- if (sum(gp$gene_type == "tRNA") > 0) {
      trna_intervals <- IRanges::IRanges(start = gp$txStart[gp$gene_type == "tRNA"] + 1, 
                                         end = gp$txEnd[gp$gene_type == "tRNA"])
      sum(IRanges::width(IRanges::reduce(trna_intervals)))
    } else { 0 }
    
    len_rRNA <- if (sum(gp$gene_type == "rRNA") > 0) {
      rrna_intervals <- IRanges::IRanges(start = gp$txStart[gp$gene_type == "rRNA"] + 1, 
                                         end = gp$txEnd[gp$gene_type == "rRNA"])
      sum(IRanges::width(IRanges::reduce(rrna_intervals)))
    } else { 0 }
    
    len_ncRNA <- if (sum(gp$gene_type == "ncRNA") > 0) {
      ncrna_intervals <- IRanges::IRanges(start = gp$txStart[gp$gene_type == "ncRNA"] + 1, 
                                          end = gp$txEnd[gp$gene_type == "ncRNA"])
      sum(IRanges::width(IRanges::reduce(ncrna_intervals)))
    } else { 0 }
    
    gp_cds <- gp[gp$gene_type == "CDS", ]
    
    if (nrow(gp_cds) == 0) {
      log_message("No CDS genes found, Trying pattern matching...")
      possible_cds <- gp[grepl("ND[1-6]|CO[I-III]|COX[1-3]|ATP[6-8]|CYB", gp$gene_name, ignore.case = TRUE), ]
      if (nrow(possible_cds) > 0) {
        log_message("Found possible CDS genes: ", paste(possible_cds$gene_name, collapse = ", "))
        gp_cds <- possible_cds
        gp_cds$gene_type <- "CDS"
      } else {
        stop("No CDS genes found")
      }
    }
    
    outdir <- file.path(out_root, "paml_runs")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    log_message("CDS gene count: ", nrow(gp_cds))
    
    SN_tbl <- data.frame(gene = character(), S = numeric(), N = numeric(), stringsAsFactors = FALSE)
    
    for (i in 1:nrow(gp_cds)) {
      row <- gp_cds[i, , drop = FALSE]  # Keep as a data frame
      safe_gene <- gsub("[^A-Za-z0-9]", "_", row$gene_name)
      gdir <- file.path(outdir, safe_gene)
      dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
      
      log_message("Processing gene: ", row$gene_name)
      
      tryCatch({
        dna_cds <- extract_cds_dna_simple(mt, row)
        log_message("  CDS length: ", length(dna_cds), " bp")
        
        files <- write_codeml_files(gdir, safe_gene, dna_cds, icode)
        SN <- run_codeml_and_parse_SN(gdir, files$ctlfile, codeml)
        
        SN_tbl <- rbind(SN_tbl, data.frame(
          gene = row$gene_name, 
          S = SN["S"], 
          N = SN["N"],
          stringsAsFactors = FALSE
        ))
        
        if (!is.na(SN["S"]) && !is.na(SN["N"])) {
          log_message("  Success: S = ", round(SN["S"], 1), ", N = ", round(SN["N"], 1))
        } else {
          log_message("  failed: Cannot obtain S/N values")
        }
        
      }, error = function(e) {
        log_message("  Gene processing failed: ", conditionMessage(e))
        SN_tbl <- rbind(SN_tbl, data.frame(
          gene = row$gene_name, 
          S = NA_real_, 
          N = NA_real_,
          stringsAsFactors = FALSE
        ))
      })
    }
    
    total_S <- sum(SN_tbl$S, na.rm = TRUE)
    total_N <- sum(SN_tbl$N, na.rm = TRUE)
    
    log_message("PAML analysis complete")
    log_message("Total synonymous sites (S): ", round(total_S, 1))
    log_message("Total nonsynonymous sites (N): ", round(total_N, 1))
    
    write.table(SN_tbl, file.path(out_root, "SN_by_gene.txt"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    bg_sites <- c(
      "non-synonymous" = total_N,
      "synonymous" = total_S,
      "tRNA" = len_tRNA,
      "rRNA" = len_rRNA,
      "ncRNA" = len_ncRNA,
      "D-loop" = len_dloop,
      "intergenic" = len_intergenic
    )
    
    bg_df <- data.frame(
      region = names(bg_sites),
      sites = bg_sites,
      stringsAsFactors = FALSE
    )
    write.table(bg_df, file.path(out_root, "background_sites.txt"), 
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    list(
      SN_by_gene = SN_tbl,
      background_sites = bg_sites,
      mt_length = mt_len,
      gene_annotations = gp
    )
    
  }, error = function(e) {
    log_message("Analysis failed: ", conditionMessage(e))
    stop(e)
  }, finally = {
    if (exists("log_con")) close(log_con)
  })
}


species <- c('human', 'mouse', 'rat', 'chicken', 'cobra', 'sturgeon', 'stingray')

mt_fasta <- c(
  "human" = "/PATH/TO/pbmc.mt/SNP/human/mgatk/downstream/rm.germline/human.germline.rm.MT.genome.fa",
  "mouse" = "/PATH/TO/pbmc.mt/ref/mouse/new/mm39_mito.fa",
  "rat" = "/PATH/TO/pbmc.mt/SNP/rat/mgatk/downstream/rm.germline/rat.germline.rm.MT.genome.fa",
  "chicken" = "/PATH/TO/pbmc.mt/SNP/chicken/mgatk/downstream/rm.germline/chicken.germline.rm.MT.genome.fa",
  "cobra" = "/PATH/TO/pbmc.mt/SNP/cobra/mgatk/downstream/rm.germline/cobra.germline.rm.MT.genome.fa",
  "sturgeon" = "/PATH/TO/pbmc.mt/SNP/sturgeon/mgatk/downstream/rm.germline/sturgeon.germline.rm.MT.genome.fa",
  "stingray" = "/PATH/TO/pbmc.mt/SNP/stingray/mgatk/downstream/rm.germline/stingray.germline.rm.MT.genome.fa"
)

gene_pred_file_list <- c(
  "human" = "/PATH/TO/annovar/humandb/hg38_chrM_ensGene.txt",
  "mouse" = "/PATH/TO/annovar/mm39_MT_db/mm39_MT_ensGene.txt",
  "rat" = "/PATH/TO/annovar/rat_MT_db/rn7_MT_ensGene.txt",
  "chicken" = "/PATH/TO/annovar/chicken_MT_db/gg1_MT_ensGene.txt",
  "cobra" = "/PATH/TO/annovar/cobra_MT_db/na_MT_ensGene.txt",
  "sturgeon" = "/PATH/TO/annovar/sturgeon_MT_db/as_CNCB_MT_ensGene.txt",
  "stingray" = "/PATH/TO/annovar/stingray_MT_db/pl_MT_ensGene.txt"
)

codeml_bin <- "/PATH/TO/codeml"
out_root_base <- "/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/PAML"
icode_val <- 1  # mammalian mitochondrialGenetic code table


dir.create(out_root_base, recursive = TRUE, showWarnings = FALSE)

all_results <- list()

for (sp in species) {
  message("\n", strrep("=", 50))
  message("Starting species: ", sp)
  message(strrep("=", 50))
  
  fasta_path <- unname(mt_fasta[sp])
  genepred_path <- unname(gene_pred_file_list[sp])
  
  if (!file.exists(fasta_path)) {
    warning("FASTA file does not exist: ", fasta_path, ", Skipping species ", sp)
    next
  }
  if (!file.exists(genepred_path)) {
    warning("GenePred file does not exist: ", genepred_path, ", Skipping species ", sp)
    next
  }
  
  out_root <- file.path(out_root_base, sp)
  
  res <- tryCatch({
    compute_background_with_paml(
      fasta_path = fasta_path,
      genepred_path = genepred_path,
      out_root = out_root,
      codeml = codeml_bin,
      icode = icode_val
    )
  }, error = function(e) {
    warning("species ", sp, " Analysis failed: ", conditionMessage(e))
    NULL
  })
  
  if (!is.null(res)) {
    all_results[[sp]] <- res
    message("species ", sp, " analysis complete, Results saved to: ", out_root)
  }
}


message("\nGenerating cross-species summary...")

summary_results <- data.frame()
for (sp in names(all_results)) {
  res <- all_results[[sp]]
  summary_results <- rbind(summary_results, data.frame(
    species = sp,
    mt_length = res$mt_length,
    non_synonymous_sites = res$background_sites["non-synonymous"],
    synonymous_sites = res$background_sites["synonymous"],
    tRNA_sites = res$background_sites["tRNA"],
    rRNA_sites = res$background_sites["rRNA"],
    ncRNA_sites = res$background_sites["ncRNA"],
    dloop_sites = res$background_sites["D-loop"],
    intergenic_sites = res$background_sites["intergenic"],
    total_coding_sites = res$background_sites["non-synonymous"] + res$background_sites["synonymous"],
    total_non_coding_sites = res$background_sites["tRNA"] + res$background_sites["rRNA"] + 
      res$background_sites["ncRNA"] + res$background_sites["D-loop"] + res$background_sites["intergenic"],
    stringsAsFactors = FALSE
  ))
}

summary_file <- file.path(out_root_base, "background_sites_summary.csv")
write.csv(summary_results, summary_file, row.names = FALSE)
message("Cross-species summary saved to: ", summary_file)

all_sn_results <- data.frame()
for (sp in names(all_results)) {
  res <- all_results[[sp]]
  sn_df <- res$SN_by_gene
  sn_df$species <- sp
  all_sn_results <- rbind(all_sn_results, sn_df)
}

sn_summary_file <- file.path(out_root_base, "all_species_SN_results.csv")
write.csv(all_sn_results, sn_summary_file, row.names = FALSE)
message("allspeciesSNresultsaveto: ", sn_summary_file)

message("\nMulti-species mitochondrial analysis complete！")
message("Root output directory: ", out_root_base)
message("\nSpecies summary:")
print(summary_results)

completion_file <- file.path(out_root_base, "ANALYSIS_COMPLETE.txt")
writeLines(paste("analysis completeTime:", Sys.time()), completion_file)
message("Analysis completion marker file: ", completion_file)


message("\nOutput directory structure:")
message(out_root_base)
message("├── background_sites_summary.csv    # Cross-species background-site summary")
message("├── all_species_SN_results.csv      # allspeciesSNResults summary")
message("├── ANALYSIS_COMPLETE.txt           # analysis completion marker")
message("└── [speciesname]/")
message("    ├── analysis.log                 # Analysis log")
message("    ├── analysis_summary.txt         # analysis summary")
message("    ├── SN_by_gene.txt              # geneSNresult")
message("    ├── background_sites.txt        # Background-site summary")
message("    ├── gene_annotations.txt        # Gene annotations")
message("    └── paml_runs/                  # PAMLRun files")
message("        └── [Gene name]/")
message("            ├── [Gene name]_same_codon.txt  # sequence file")
message("            ├── [Gene name]_same.tree       # Tree file")
message("            ├── [Gene name]_same.ctl        # codemlControl file")
message("            └── mlc                         # codemlOutput")


message("\nCheckingPAMLRun results...")
for (sp in names(all_results)) {
  res <- all_results[[sp]]
  valid_results <- sum(!is.na(res$SN_by_gene$S) & !is.na(res$SN_by_gene$N))
  total_genes <- nrow(res$SN_by_gene)
  message(sprintf("species %s: %d/%d genes obtained successfullyS/Nvalue", sp, valid_results, total_genes))
  
  if (valid_results > 0) {
    total_S <- sum(res$SN_by_gene$S, na.rm = TRUE)
    total_N <- sum(res$SN_by_gene$N, na.rm = TRUE)
    message(sprintf("  Total S value: %.1f, Total N value: %.1f, N/Sratio: %.3f", 
                    total_S, total_N, total_N/total_S))
  }
}
