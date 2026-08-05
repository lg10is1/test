#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_BASE="${PROJECT_BASE:-${REPO_DIR}/data}"
FIGURE7_BASE="${FIGURE7_BASE:-${PROJECT_BASE}/Figure7_heritability}"
TR_GRM_PREFIX="${TR_GRM_PREFIX:-${FIGURE7_BASE}/tr_grm/ngstr}"
REMOVE_FILE="${REMOVE_FILE:-${FIGURE7_BASE}/tr_grm/remove_sample}"
OUT_DIR="${OUT_DIR:-${REPO_DIR}/results/heritability/tr_pca}"
PHENO_FILE="${PHENO_FILE:-${PROJECT_BASE}/TGS_callset/Pangenie_v3/06.gwas/SCZ_pheno.txt}"
BATCH_FILE="${BATCH_FILE:-${PROJECT_BASE}/TGS_callset/Pangenie_v3/06.gwas/SCZ_batch.txt}"
GCTA_BIN="${GCTA_BIN:-gcta64}"
RSCRIPT_BIN="${RSCRIPT_BIN:-Rscript}"
THREADS="${THREADS:-28}"
PCA_COUNT="${PCA_COUNT:-20}"

usage() {
  cat <<'EOF'
Usage:
  bash run_tr_pca_remove_samples.sh [options]

Options:
  --tr-grm-prefix PREFIX  Existing TR GRM prefix
  --remove-file FILE      Samples to remove; one-column IID or two-column FID IID
  --out-dir DIR           Output directory
  --pheno FILE            Phenotype file used only for plot annotation
  --batch FILE            Batch file used only for plot annotation
  --threads N             GCTA threads (default: 28)
  --pca-count N           Number of PCs (default: 20; must be >=4)
  --gcta PATH             GCTA executable (default: gcta64)
  --rscript PATH          Rscript executable (default: Rscript)
  -h, --help              Show help
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tr-grm-prefix) TR_GRM_PREFIX="$2"; shift 2 ;;
    --remove-file) REMOVE_FILE="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --pheno) PHENO_FILE="$2"; shift 2 ;;
    --batch) BATCH_FILE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --pca-count) PCA_COUNT="$2"; shift 2 ;;
    --gcta) GCTA_BIN="$2"; shift 2 ;;
    --rscript) RSCRIPT_BIN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if ! [[ "$PCA_COUNT" =~ ^[0-9]+$ ]] || (( PCA_COUNT < 4 )); then
  echo "[ERROR] --pca-count must be an integer >=4" >&2
  exit 2
fi

for file in "${TR_GRM_PREFIX}.grm.bin" "${TR_GRM_PREFIX}.grm.N.bin" \
  "${TR_GRM_PREFIX}.grm.id" "$REMOVE_FILE" "$PHENO_FILE" "$BATCH_FILE"; do
  [[ -s "$file" ]] || { echo "[ERROR] Missing or empty file: $file" >&2; exit 1; }
done

mkdir -p "$OUT_DIR"
KEEP_FILE="${OUT_DIR}/tr_pca.keep"
MATCHED_REMOVE="${OUT_DIR}/remove_sample.matched.tsv"
UNMATCHED_REMOVE="${OUT_DIR}/remove_sample.unmatched.tsv"
PCA_PREFIX="${OUT_DIR}/tr_after_remove"
LOG_FILE="${OUT_DIR}/tr_after_remove.pca.log"
: > "$MATCHED_REMOVE"
: > "$UNMATCHED_REMOVE"

# Accept one-column IID or two-column FID/IID input. Headers, blank lines,
# comments, commas, spaces, tabs, and Windows CRLF are tolerated.
awk -v matched="$MATCHED_REMOVE" '
  BEGIN { FS="[ ,\t]+"; OFS="\t" }
  NR==FNR {
    gsub(/\r/, "")
    if ($0 ~ /^[[:space:]]*$/ || $0 ~ /^[[:space:]]*#/) next
    u1=toupper($1); u2=toupper($2)
    if (u1=="FID" || u1=="IID" || u2=="IID") next
    if (NF >= 2) remove_pair[$1 SUBSEP $2]=1
    else remove_iid[$1]=1
    next
  }
  {
    fid=$1; iid=$2
    if ((fid SUBSEP iid) in remove_pair || iid in remove_iid) print fid, iid > matched
    else print fid, iid
  }
' "$REMOVE_FILE" "${TR_GRM_PREFIX}.grm.id" > "$KEEP_FILE"

total_n="$(awk 'NF >= 2 {n++} END {print n+0}' "${TR_GRM_PREFIX}.grm.id")"
keep_n="$(awk 'NF >= 2 {n++} END {print n+0}' "$KEEP_FILE")"
removed_n="$(awk 'NF >= 2 {n++} END {print n+0}' "$MATCHED_REMOVE" 2>/dev/null || true)"
removed_n="${removed_n:-0}"

if (( removed_n == 0 )); then
  echo "[ERROR] No remove_sample entry matched ${TR_GRM_PREFIX}.grm.id" >&2
  exit 1
fi
if (( keep_n <= PCA_COUNT )); then
  echo "[ERROR] Only ${keep_n} samples remain; cannot compute ${PCA_COUNT} PCs" >&2
  exit 1
fi

# Record requested removal entries that did not match any GRM IID/pair.
awk '
  BEGIN { FS="[ ,\t]+"; OFS="\t" }
  NR==FNR { gsub(/\r/, ""); if (NF>=2) pair[$1 SUBSEP $2]=1; iid[$2]=1; next }
  {
    gsub(/\r/, "")
    if ($0 ~ /^[[:space:]]*$/ || $0 ~ /^[[:space:]]*#/) next
    u1=toupper($1); u2=toupper($2)
    if (u1=="FID" || u1=="IID" || u2=="IID") next
    if (NF>=2) { if (!(($1 SUBSEP $2) in pair)) print $1, $2 }
    else { if (!($1 in iid)) print $1 }
  }
' "${TR_GRM_PREFIX}.grm.id" "$REMOVE_FILE" > "$UNMATCHED_REMOVE"

echo "[INFO] TR GRM samples: ${total_n}"
echo "[INFO] Matched samples removed: ${removed_n}"
echo "[INFO] Samples retained for PCA: ${keep_n}"

{
  printf '[CMD] '
  printf '%q ' "$GCTA_BIN" --grm "$TR_GRM_PREFIX" --keep "$KEEP_FILE" \
    --pca "$PCA_COUNT" --thread-num "$THREADS" --out "$PCA_PREFIX"
  printf '\n'
} > "$LOG_FILE"

"$GCTA_BIN" --grm "$TR_GRM_PREFIX" --keep "$KEEP_FILE" \
  --pca "$PCA_COUNT" --thread-num "$THREADS" --out "$PCA_PREFIX" \
  >> "$LOG_FILE" 2>&1

[[ -s "${PCA_PREFIX}.eigenvec" ]] || { echo "[ERROR] Missing ${PCA_PREFIX}.eigenvec; see $LOG_FILE" >&2; exit 1; }
[[ -s "${PCA_PREFIX}.eigenval" ]] || { echo "[ERROR] Missing ${PCA_PREFIX}.eigenval; see $LOG_FILE" >&2; exit 1; }

"$RSCRIPT_BIN" - "${PCA_PREFIX}.eigenvec" "${PCA_PREFIX}.eigenval" \
  "$PHENO_FILE" "$BATCH_FILE" "$OUT_DIR" "$PCA_COUNT" "$total_n" "$removed_n" "$keep_n" <<'RSCRIPT'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) stop("Internal argument error")
eigenvec_file <- args[[1]]
eigenval_file <- args[[2]]
pheno_file <- args[[3]]
batch_file <- args[[4]]
out_dir <- args[[5]]
pca_count <- as.integer(args[[6]])
total_n <- as.integer(args[[7]])
removed_n <- as.integer(args[[8]])
keep_n <- as.integer(args[[9]])

read_annotation <- function(file, value_name) {
  x <- fread(file, header = "auto")
  header_ok <- ncol(x) >= 2L && all(toupper(names(x)[1:2]) %in% c("FID", "IID", "#FID", "#IID"))
  if (!header_ok) x <- fread(file, header = FALSE)
  if (ncol(x) < 3L) stop("Expected at least 3 columns: ", file)
  x <- x[, 1:3]
  setnames(x, c("FID", "IID", value_name))
  x[, `:=`(FID = as.character(FID), IID = as.character(IID))]
  if (anyDuplicated(x[, .(FID, IID)])) stop("Duplicated IDs: ", file)
  x
}

pc <- fread(eigenvec_file, header = FALSE)
if (ncol(pc) < pca_count + 2L) stop("Eigenvec has fewer columns than requested PCs")
pc <- pc[, seq_len(pca_count + 2L), with = FALSE]
setnames(pc, c("FID", "IID", paste0("PC", seq_len(pca_count))))
pc[, `:=`(FID = as.character(FID), IID = as.character(IID))]

pheno <- read_annotation(pheno_file, "phenotype")
batch <- read_annotation(batch_file, "batch")
pc <- merge(pc, pheno, by = c("FID", "IID"), all.x = TRUE)
pc <- merge(pc, batch, by = c("FID", "IID"), all.x = TRUE)
pc[, phenotype_label := fifelse(as.character(phenotype) %in% c("2", "case", "CASE"), "Case",
  fifelse(as.character(phenotype) %in% c("1", "control", "CONTROL"), "Control", as.character(phenotype)))]
pc[, batch_label := as.factor(batch)]

vals <- as.numeric(fread(eigenval_file, header = FALSE)[[1]])
positive_vals <- pmax(vals, 0)
denom <- sum(positive_vals)
eigen <- data.table(
  PC = paste0("PC", seq_along(vals)),
  eigenvalue = vals,
  percent = if (denom > 0) 100 * positive_vals / denom else NA_real_
)

fwrite(pc, file.path(out_dir, "tr_pc_values.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(pc, file.path(out_dir, "tr_pc_values.csv"), quote = TRUE, na = "NA")
fwrite(eigen, file.path(out_dir, "tr_pc_eigenvalues.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(eigen, file.path(out_dir, "tr_pc_eigenvalues.csv"), quote = TRUE, na = "NA")
summary <- data.table(total_grm_samples = total_n, removed_samples = removed_n,
  retained_pca_samples = keep_n, pca_count = pca_count,
  n_with_phenotype = sum(!is.na(pc$phenotype)), n_with_batch = sum(!is.na(pc$batch)))
fwrite(summary, file.path(out_dir, "tr_pc_summary.tsv"), sep = "\t", quote = FALSE)
fwrite(summary, file.path(out_dir, "tr_pc_summary.csv"), quote = TRUE)

pct <- setNames(eigen$percent, eigen$PC)
axis_lab <- function(nm) if (nm %in% names(pct) && is.finite(pct[[nm]])) sprintf("%s (%.2f%%)", nm, pct[[nm]]) else nm
make_plot <- function(x, y, colour_col, colour_title) {
  ggplot(pc, aes(x = .data[[x]], y = .data[[y]], colour = .data[[colour_col]])) +
    geom_point(size = 1.8, alpha = 0.8, na.rm = TRUE) +
    labs(title = "TR PCA after removing outlier samples",
      subtitle = sprintf("N=%d retained; %d removed", keep_n, removed_n),
      x = axis_lab(x), y = axis_lab(y), colour = colour_title) +
    theme_bw(base_size = 11)
}

max_pair_pc <- min(20L, pca_count)
pc_pairs <- data.table(
  x = paste0("PC", seq(1L, max_pair_pc - 1L, by = 2L)),
  y = paste0("PC", seq(2L, max_pair_pc, by = 2L))
)
plots <- list()
for (i in seq_len(nrow(pc_pairs))) {
  x <- pc_pairs$x[[i]]
  y <- pc_pairs$y[[i]]
  pair_name <- paste0(x, "_", y)
  plots[[paste0(pair_name, ".case_control")]] <- make_plot(x, y, "phenotype_label", "Phenotype")
  plots[[paste0(pair_name, ".batch")]] <- make_plot(x, y, "batch_label", "Batch")
}
scree <- eigen[seq_len(min(20L, .N))]
plots$scree <- ggplot(scree, aes(x = factor(PC, levels = PC), y = percent)) +
  geom_col(fill = "#356AA0") +
  labs(title = "TR PCA scree plot after sample removal", x = "PC", y = "Variance explained (%)") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

for (nm in names(plots)) {
  ggsave(file.path(out_dir, paste0("TR_", nm, ".png")), plots[[nm]], width = 8.2, height = 6.2, dpi = 180)
  ggsave(file.path(out_dir, paste0("TR_", nm, ".pdf")), plots[[nm]], width = 8.2, height = 6.2)
}
pdf(file.path(out_dir, "TR_all_pc_plots.pdf"), width = 8.2, height = 6.2, onefile = TRUE)
for (p in plots) print(p)
dev.off()
message("[DONE] TR PC values and plots: ", out_dir)
RSCRIPT

echo "[DONE] Outputs: ${OUT_DIR}"
echo "[DONE] PC values: ${OUT_DIR}/tr_pc_values.tsv and .csv"
echo "[DONE] Matched removals: ${MATCHED_REMOVE}"
