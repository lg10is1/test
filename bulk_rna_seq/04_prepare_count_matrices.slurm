#!/usr/bin/env bash
#SBATCH --job-name=bulk_rna_counts
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --output=task_file/%x_%j.out
#SBATCH --error=task_file/%x_%j.err

set -euo pipefail
: "${BULK_RNA_WORK_DIR:?Set BULK_RNA_WORK_DIR}"
: "${EOSCZ_REFERENCE_DIR:?Set EOSCZ_REFERENCE_DIR}"

module load miniconda3
source activate "${CONDA_ENV:-rnaseq}"
cd "$BULK_RNA_WORK_DIR"
mkdir -p count/post_rna_seq_analysis task_file

output_file="gtf_files_list.txt"
rm -f -- "$output_file"
find ballgown -type f -name '*.gtf' -print0 | while IFS= read -r -d '' file; do
  filename="${file##*/}"
  filename="${filename%.gtf}"
  printf '%s\t%s\n' "$filename" "$file" >> "$output_file"
done
[[ -s "$output_file" ]] || { echo "No GTF files found under ballgown" >&2; exit 1; }

python3 prepare_count_matrices.py -i "$output_file" \
  -g count/gene_count_matrix.csv -t count/transcript_count_matrix.csv
cd count
grep 'protein_coding' "$EOSCZ_REFERENCE_DIR/T2T-CHM13V2.0_UCSC.sorted.gtf" \
  | awk -F'\t' '$3 == "gene"' | awk -F';' '{print $1}' | cut -f9 \
  | sed 's/"//g' | sed 's/gene_id //g' > protein_coding_genes.txt
grep 'gene_id' gene_count_matrix.csv > gene_count_matrix_protein_coding.csv
grep -Ff protein_coding_genes.txt gene_count_matrix.csv >> gene_count_matrix_protein_coding.csv
cp gene_count_matrix_protein_coding.csv post_rna_seq_analysis/
grep 'CAMTA1' "$EOSCZ_REFERENCE_DIR/T2T-CHM13V2.0_UCSC.sorted.gtf" > camta1_reference_records.txt
