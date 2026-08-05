# Single-cell RNA-seq

Marker-based cell-type annotation (snapseed) and Seurat downstream analysis.

## Requirements

- Python 3 with `scanpy`, `anndata`, `pandas`, `scipy`, `snapseed`, and `hnoca`
- R with `Seurat`, `scCustomize`, `FLASHMM`, `data.table`, `dplyr`, `ComplexHeatmap`, `clusterProfiler`, `org.Hs.eg.db`

## Scripts

| Script | Purpose |
|---|---|
| `snapseed_annotation.py` | Loads a 10x-style matrix (`matrix.mtx`, `features.tsv`, `barcodes.tsv`, `meta.csv`, `matrix_log.mtx`) and annotates cell types from a marker YAML file. |
| `seurat_analysis.R` | Integration, clustering, marker scoring, and downstream analyses. Sample phenotype mapping is read from a de-identified TSV with columns `sample_id`, `time`, `mutation_group`, `disease_group`, `trio` (see `../examples/synthetic/single_cell_sample_metadata.example.tsv`). |

## Usage

```bash
EOSCZ_SCRNA_INPUT_DIR=/path/to/single_cell_matrix \
EOSCZ_SNAPSEED_MARKERS=/path/to/marker_definitions.yaml \
python3 snapseed_annotation.py

EOSCZ_SCRNA_INPUT_DIR=/path/to/seurat_objects \
EOSCZ_SCRNA_SAMPLE_METADATA=/path/to/sample_metadata.tsv \
EOSCZ_SCRNA_OUTPUT_DIR=seurat_output \
Rscript seurat_analysis.R
```
