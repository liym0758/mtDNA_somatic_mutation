library(utils)

create_annovar_input <- function(input_file, 
                                 final_variants, 
                                 chr_name, 
                                 output_file) {
  variant_stats_filter <- input_file %>% 
    dplyr::filter(variant %in% final_variants)
  
  nuc_split <- strsplit(variant_stats_filter$nucleotide, ">")
  ref <- sapply(nuc_split, `[`, 1)
  alt <- sapply(nuc_split, `[`, 2)
  
  AnnovarInput <- data.frame(
    chr = rep(chr_name, nrow(variant_stats_filter)),
    start = variant_stats_filter$position, 
    end   = variant_stats_filter$position, 
    ref   = ref, 
    alt   = alt
  )
  
  write.table(AnnovarInput, 
              file = output_file, 
              sep = "\t",
              quote = FALSE,
              row.names = FALSE, 
              col.names = FALSE)
}




base_path <- "/PATH/TO/pbmc.mt/SNP/"
chr_name <- c("human" = "chrM",
              "rat" = "MT",
              "mouse" = "MT",
              "chicken" = "MT",
              "cobra" = "MT",
              "sturgeon" = "MT",
              "stingray" = "MT")
species_list <- c("human","mouse","rat","chicken","cobra","sturgeon","stingray")
variant.list <- readRDS("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Frequency/VAF_histogram/fixedDepth_VAF/new_human/Final_Variant_List_7species_260128.RDS")

for (i in 1:7){
  create_annovar_input(
    input_file =  read_tsv(paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Filter_pct_conf_detected_cell/final_filter/",species_list[i],"/variant_stats_final_filtered.tsv")),
    final_variants = variant.list[[i]],
    chr_name = chr_name[i],
    output_file = paste0("/PATH/TO/pbmc.mt/SNP/mgatk.downstream/Variant_Type_Anno/new_AnnoVar_Result/",species_list[i],"/",species_list[i],".AnnovarInput")
  )
}






