# mtDNA_somatic_mutation
Somatic mtDNA Mutations and Selection across Vertebrate Peripheral Blood Mononuclear Cells
## Analysis order

| Step | Code | Main output |
| --- | --- | --- |
| 1. Prepare reference and NUMT mask | `scripts/01_reference/example_reference.sh`, `numt_simulation/` | Masked reference and NUMT candidates |
| 2. scATAC analysis | `scripts/02_scatac/example_cellranger_atac.sh`, `example_signac_qc_clustering.R`, `example_rna_atac_label_transfer.R` | Cell Ranger BAM and annotated cells |
| 3. mtDNA read processing and read QC | `scripts/03_read_processing/` | Filtered and realigned mtDNA BAM, read statistics |
| 4. mgatk variant calling | `scripts/04_variant_calling/example_mgatk.sh`, `call_mutations_mgatk.R` | Allele matrices and variant statistics |
| 5. Variant filtering | `scripts/05_variant_filtering/variant_filtering.R` | Final variant list |
| 6. Figure analyses and plots | `scripts/07_figures/` | Figure-specific calculations, annotation, and main/supplementary plots |
