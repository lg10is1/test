# Assembly quality assessment

Meryl k-mer counting and Merqury QV evaluation of hifiasm assemblies.

## Requirements

- Meryl and Merqury (Conda environment, default name `merqury`)
- Slurm

## Scripts

| Script | Purpose |
|---|---|
| `link_assembly_fastas.sh` | Creates symlinks of assembly FASTA files with normalized names (`ASSEMBLY_SOURCE_DIR`, `ASSEMBLY_LINK_DIR`). |
| `evaluate_assembly_qv.slurm` | Runs Meryl counting and Merqury QV for one sample. |

## Usage

```bash
ASSEMBLY_SOURCE_DIR=/path/to/assemblies \
ASSEMBLY_LINK_DIR=/path/to/links \
bash link_assembly_fastas.sh

MERYL_INPUT_DIR=/path/to/reads \
MERYL_OUTPUT_DIR=/path/to/meryl \
ASSEMBLY_FASTA_DIR=/path/to/links \
MERQURY_SH=/path/to/merqury.sh \
sbatch evaluate_assembly_qv.slurm SAMPLE_ID [ASSEMBLY_ID]
```

Optional variables: `THREADS` (20), `MEMORY_GB` (80), `KMER_SIZE` (21), `CONDA_ENV` (`merqury`), `MERYL_ANALYSIS_DIR`.
