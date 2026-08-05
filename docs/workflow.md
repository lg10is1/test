# Workflow guide

## Scope

The source collection contains related analysis modules, not one fully automated workflow. The sequence below reflects code relationships and file semantics; it is not a substitute for an author-approved protocol.

## Upstream paths

### Long-read assembly path

1. `genome_assembly/convert_bam_to_fastq.slurm` converts authorized PacBio BAM inputs to FASTQ.
2. `genome_assembly/summarize_fastq.slurm` summarizes FASTQ files with SeqKit.
3. `genome_assembly/assemble_hifiasm_array.sh` generates and submits one hifiasm job per de-identified FASTQ.
4. `assembly_quality/link_assembly_fastas.sh` prepares local assembly links.
5. `assembly_quality/evaluate_assembly_qv.slurm` runs Meryl and Merqury quality assessment.

### Read-mapping and structural-variant path

1. Copy `read_mapping_sv_calling/config.example.sh` to the ignored `config.sh`.
2. Supply a de-identified sample list and authorized input paths.
3. `read_mapping_sv_calling/call_structural_variants.sh` converts reads as needed, aligns with pbmm2, discovers/calls with pbsv, and calls with Sniffles.

### Genotyping paths

- `pangenie_genotyping` contains graph VCF preparation, PanGenie execution, single-sample QC, PLINK conversion/QC, and a GWAS wrapper.
- `paragraph_genotyping/genotype_sample.slurm` invokes a local Paragraph `multigrmpy.py` installation for a supplied manifest and VCF.

### Reference-population and inversion path

1. `reference_population_analysis/1kgp/ont/download_igsr_ont.sh` validates the retained public IGSR manifest and, only with `--execute`, downloads and verifies the listed ONT CRAM/CRAI files.
2. `reference_population_analysis/1kgp/ngs/remap_cram_to_chm13.sh` is a resumable Slurm-array template that converts authorized GRCh38-referenced CRAMs to paired FASTQ and realigns them to an author-approved CHM13 reference.
3. `reference_population_analysis/filter_vcf_by_positions.py` retains exact CHROM/POS records from a VCF, and `list_vcf_samples.sh` extracts VCF sample names.
4. `reference_population_analysis/combine_large_inversions.py` clusters single-sample inversion VCFs and compares case/control carrier rates.

These utilities are not automatically connected to the EOSCZ production workflow. Confirm whether and where the external reference data and inversion analysis belong in the manuscript-specific analysis graph.

## Association and downstream analyses

The numbered scripts in `snv_sv_gwas/scripts` prepare variant sets, run or merge GWAS results, clump signals, generate Manhattan/QQ and locus plots, add local GWAS Catalog annotations, evaluate LD, and run selected molecular-QTL analyses. `snv_sv_gwas/scripts/run_pipeline.sh` documents an internal ordering but still requires local data and author-confirmed configuration.

CNV, C4, MHC/HLA, tandem-repeat, heritability, PRS, bulk RNA-seq, and single-cell modules consume module-specific derived tables. These modules can be run independently after their inputs have been generated and reviewed.

## Dependency boundaries

- Assemblies feed assembly quality and can support downstream variant discovery.
- Mapped reads feed pbsv and Sniffles calling.
- Variant call sets feed PanGenie/Paragraph genotyping, QC, GWAS, LD, and annotation.
- CNV and HLA analyses expect precomputed per-sample or per-haplotype tables.
- Transcriptomic modules are not wired into the genomic workflow by a repository-level launcher.
- Public 1000 Genomes manifests can support external reference-data preparation, but the files they describe are not repository inputs and are never downloaded automatically.
