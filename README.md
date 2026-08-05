# Supplementary analysis scripts for the EOSCZ study

This repository contains analysis scripts accompanying a study of early-onset schizophrenia (EOSCZ), including genome assembly, variant analysis, transcriptomics, DNA methylation, and genetic association analyses.


## Contents

| Directory | Analysis |
|---|---|
| `genome_assembly` | HiFi read processing and genome assembly |
| `assembly_quality` | Assembly quality assessment |
| `read_mapping_sv_calling` | Read mapping and structural-variant calling |
| `pangenie_genotyping` | Pangenome-based genotyping and quality control |
| `paragraph_genotyping` | Paragraph structural-variant genotyping |
| `snv_sv_gwas` | SNV/SV association, LD, annotation, and locus analyses |
| `copy_number_variants` | Copy-number analyses and figure generation |
| `mhc_hla` | MHC/HLA analyses |
| `complement_c4` | Complement C4 copy-number analyses |
| `tandem_repeats` | Tandem-repeat analyses |
| `5mc_calling` | CpG methylation calling and filtering |
| `bulk_rna_seq` | Bulk RNA-seq processing |
| `single_cell_rna_seq` | Single-cell RNA-seq analysis |
| `eQTL` | cis-eQTL analysis and fine-mapping |
| `heritability` | Heritability analyses |
| `polygenic_risk_scores` | Polygenic-risk and rare-variant burden analyses |
| `reference_population_analysis` | Reference-population utilities and comparisons |

## Notes

- Analysis methods and parameter choices should be interpreted together with the manuscript and its supplementary methods.
- Individual-level genomic, clinical, and molecular data are not included because they are subject to controlled-access and privacy requirements.
- Scripts were developed for the study's Linux/HPC environment. Required programs and packages are listed in script headers or module-level README files where applicable.
- Intermediate files and study results are not included.

## Citation

Please cite the associated article when using these scripts. Full citation details will be added after publication.

## License

Code in this repository is available under the [MIT License](LICENSE).
