# Input formats

This document describes formats inferred from retained code. Validate every schema against the author-approved methods before analysis.

| Input | Typical modules | Minimum inferred content | Privacy |
|---|---|---|---|
| FASTQ / BAM | assembly, mapping, SV calling | Sequencing reads and valid headers/read groups | Controlled individual-level data; not included |
| FASTA plus indexes | assembly QC, mapping, genotyping, annotation | Module-appropriate T2T-CHM13 or GRCh38 reference | Supply locally and record checksums |
| GTF/GFF | bulk RNA-seq | Standard tab-delimited records; count helper expects `gene_id`, `transcript_id`, and `cov` attributes plus StringTie metadata | Synthetic example included |
| VCF/BCF | genotyping, GWAS, LD, annotation | Standard variant fields; some scripts expect named INFO or genotype fields | Often individual-level controlled data |
| BED | tandem-repeat and interval filters | Chromosome, start, end; module-specific extra columns may be required | Review before sharing |
| Sample manifest/list | Slurm arrays, Paragraph, mapping | One de-identified sample identifier per row or tool-specific manifest columns | Keep untracked; no real manifest included |
| PLINK files | GWAS, heritability, PRS | Tool-standard BIM/BED/FAM or PGEN/PVAR/PSAM sets | Individual-level files are controlled |
| CSV/TSV/XLSX tables | CNV, C4, MHC/HLA, PRS, downstream plots | Script-specific named columns | Study tables are absent; inspect the consuming script |
| Matrix Market plus metadata | single-cell annotation | `matrix.mtx`, `features.tsv`, `barcodes.tsv`, `meta.csv`, and `matrix_log.csv` | Individual/cell metadata may be controlled |
| Marker YAML | single-cell annotation | Nested cell types with `marker_genes` and optional `subtypes` | Supply locally |
| GWAS Catalog zip | local annotation | Zip containing an associations TSV with mapped-gene and trait columns | Public source may be used, but the file is not downloaded automatically |
| Exact-position TSV | reference-population VCF filtering | Tab-separated `CHROM` and one-based `POS` in the first two columns | Do not include study linkage fields |
| Single-sample inversion VCF | inversion clustering | Filename `<sample>_filtered.vcf[.gz]`; standard VCF columns, one sample genotype column, and `END` or `SVLEN` in INFO | Controlled individual-level data; not included |
| Case/control sample lists | inversion clustering | One unique de-identified sample ID per line; lists must not overlap | Controlled linkage data; keep untracked |
| Public IGSR URL/MD5 manifests | reference-population download | One approved HTTPS URL per row and one `md5sum` row with the matching basename | Included public external-data metadata |
| CRAM path list | CHM13 remapping template | One authorized local CRAM path per line; Slurm array index is one-based unless `OFFSET` is set | Local paths and sequencing data are controlled; keep untracked |

## Synthetic GTF example

`examples/synthetic/gtf_files.tsv` has two whitespace-delimited columns: a synthetic sample ID and a repository-relative GTF path. The GTF files contain invented coordinates, gene IDs, transcript IDs, and coverage values. They are designed only for the count-matrix helper.
