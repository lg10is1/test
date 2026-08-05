# EOSCZ analysis code

## Overview

Analysis workflows for early-onset schizophrenia whole-genome, transcriptomic, methylation, and genetic-association studies. Scripts are organized by analysis module and use command-line arguments, environment variables, or example configuration files for local paths and computing resources.

## Modules

| Directory | Analysis |
|---|---|
| `5mc_calling` | pb-CpG-tools 5mC calling, chromosome-level matrix merging, and CpG filtering |
| `assembly_quality` | Meryl and Merqury assembly quality assessment |
| `bulk_rna_seq` | Read alignment, transcript assembly, quantification, and count matrices |
| `complement_c4` | C4 haplotype and sample copy-number analysis |
| `copy_number_variants` | CNV association, burden, pathway, PCA, and figure workflows |
| `eQTL` | RNA phenotype preparation, nominal cis-eQTL mapping, and SuSiEx fine-mapping |
| `genome_assembly` | BAM-to-FASTQ conversion, FASTQ summaries, and hifiasm assembly |
| `heritability` | GCTA-based heritability preparation, task generation, and result summaries |
| `mhc_hla` | HLA subtype, haplotype, association, and visualization workflows |
| `pangenie_genotyping` | Graph VCF preparation, PanGenie genotyping, QC, and GWAS |
| `paragraph_genotyping` | Paragraph multi-sample genotyping |
| `polygenic_risk_scores` | Common-variant PRS, rare-variant burden, and combined analyses |
| `read_mapping_sv_calling` | PacBio read mapping and pbsv/Sniffles structural-variant calling |
| `reference_population_analysis` | VCF utilities, inversion analysis, IGSR manifests, and CHM13 remapping |
| `single_cell_rna_seq` | Marker-based annotation and Seurat analysis |
| `snv_sv_gwas` | SNV/SV GWAS, clumping, annotation, LD, locus plots, and meQTL analysis |
| `tandem_repeats` | Tandem-repeat logistic association analysis |

Supporting files:

- `configs/config.example.yaml`: shared path and resource configuration example.
- `docs/workflow.md`: workflow relationships and execution order.
- `docs/input_formats.md`: expected input schemas.
- `docs/configuration.md`: configuration variables.
- `examples/synthetic/`: synthetic inputs for format and command checks.

## Workflow

1. Configure reference genomes, annotations, authorized input paths, output paths, and compute resources.
2. Assemble genomes or map reads and call structural variants.
3. Run PanGenie or Paragraph genotyping and module-specific QC.
4. Generate 5mC, RNA-seq, single-cell, CNV, C4, HLA, or tandem-repeat features.
5. Run GWAS, cis-eQTL, heritability, PRS, burden, and fine-mapping analyses.
6. Generate summary tables and figures.

See [docs/workflow.md](docs/workflow.md) for module connections.

## Environment

Python packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

R packages and command-line tools are listed in [software_versions.md](software_versions.md). Core workflows use Bash, Slurm, Python, R, bcftools, samtools, PLINK/PLINK2, pb-CpG-tools, QTLtools, PanGenie, Paragraph, pbsv, Sniffles, hifiasm, Meryl, Merqury, and related reference indexes.

Record the exact versions and checksums used for every analysis run.

## Configuration

Create a local configuration from the example:

```bash
cp configs/config.example.yaml configs/config.yaml
```

Keep `configs/config.yaml`, sample manifests, credentials, cluster-specific paths, and controlled-data locations outside version control. Shell modules also document required arguments and environment variables in their module README files or script headers.

## Inputs

Main input types include FASTQ, BAM/CRAM, FASTA, GTF/GFF, VCF/BCF, BED, PLINK datasets, count matrices, covariate tables, and de-identified sample manifests. Detailed fields are listed in [docs/input_formats.md](docs/input_formats.md).

The repository contains only code, configuration examples, public reference manifests, and synthetic examples. Store individual-level genomic, molecular, clinical, and linkage data in approved controlled storage.

## Usage examples

Run the synthetic transcript count-matrix example from the repository root:

```bash
mkdir -p demo_output
python3 bulk_rna_seq/prepare_count_matrices.py \
  -i examples/synthetic/gtf_files.tsv \
  -g demo_output/gene_count_matrix.csv \
  -t demo_output/transcript_count_matrix.csv
```

Filter a VCF by exact chromosome and position:

```bash
python3 reference_population_analysis/filter_vcf_by_positions.py \
  --position-file /path/to/positions.tsv \
  --vcf-file /path/to/input.vcf.gz \
  --output-vcf /path/to/filtered.vcf.gz
```

Validate the retained IGSR ONT manifest and print the download plan:

```bash
bash reference_population_analysis/1kgp/ont/download_igsr_ont.sh \
  --workdir /path/to/approved/storage
```

Module-specific commands are documented in:

- [5mc_calling/README.md](5mc_calling/README.md)
- [assembly_quality/README.md](assembly_quality/README.md)
- [bulk_rna_seq/README.md](bulk_rna_seq/README.md)
- [complement_c4/README.md](complement_c4/README.md)
- [copy_number_variants/README.md](copy_number_variants/README.md)
- [eQTL/README.md](eQTL/README.md)
- [genome_assembly/README.md](genome_assembly/README.md)
- [heritability/README.md](heritability/README.md)
- [mhc_hla/README.md](mhc_hla/README.md)
- [pangenie_genotyping/README.md](pangenie_genotyping/README.md)
- [paragraph_genotyping/README.md](paragraph_genotyping/README.md)
- [polygenic_risk_scores/README.md](polygenic_risk_scores/README.md)
- [read_mapping_sv_calling/README.md](read_mapping_sv_calling/README.md)
- [reference_population_analysis/README.md](reference_population_analysis/README.md)
- [single_cell_rna_seq/README.md](single_cell_rna_seq/README.md)
- [snv_sv_gwas/README.md](snv_sv_gwas/README.md)
- [tandem_repeats/README.md](tandem_repeats/README.md)

## Outputs

Outputs include assemblies, alignment files, variant calls, genotype matrices, methylation matrices, expression phenotypes, association statistics, burden results, and figures. Direct outputs and scheduler logs to project-specific output directories.

## License
Code in this repository is released under the MIT License.

## Contact
For questions or technical issues, please contact the corresponding author.
