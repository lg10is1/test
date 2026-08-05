import os

import scanpy as sc
import pandas as pd
from pathlib import Path
from scipy.io import mmread
import anndata as ad
indir = Path(os.environ.get("EOSCZ_SCRNA_INPUT_DIR", "/path/to/single_cell_matrix"))
X = mmread(indir / "matrix.mtx").tocsr().T           # cells × genes
feat = pd.read_csv(indir / "features.tsv", sep="\t", header=None)
bar  = pd.read_csv(indir / "barcodes.tsv", sep="\t", header=None)

adata = ad.AnnData(X=X)
gene_ids   = feat.iloc[:, 0].astype(str).values
gene_names = (feat.iloc[:, 1] if feat.shape[1] > 1 else feat.iloc[:, 0]).astype(str).values
adata.var["gene_id"] = gene_ids
adata.var["gene_name"] = gene_names
adata.var_names = gene_names  # Use gene_ids instead when required by the matrix metadata.
adata.obs_names = bar.iloc[:, 0].astype(str).values
meta = pd.read_csv(f"{indir}/meta.csv").set_index("cell").loc[adata.obs_names]
adata.obs = pd.concat([adata.obs, meta], axis=1)
logX = mmread(Path(indir) / "matrix_log.mtx").tocsr().T
assert logX.shape == adata.shape
adata.layers["lognorm"] = logX



import numpy as np
import pandas as pd
import snapseed
from hnoca.snapseed.utils import read_yaml
import hnoca.snapseed as snap
marker = read_yaml(os.environ.get("EOSCZ_SNAPSEED_MARKERS", "/path/to/hnoca_markers.yaml"))

import numpy as np
import pandas as pd
import snap

# ============================================================
# Configuration
# ============================================================
cluster_key = "seurat_clusters"
layer       = "lognorm"
score_thr   = 0.7
marker_dict = marker  
# ============================================================


def navigate_to_node(marker_dict, path_str):
    """
    Navigate to a marker node using a path such as 'neuron > excitatory_neuron'.
    Return the complete node dictionary, including markers and subtypes, or None.
    """
    parts = [p.strip() for p in path_str.split(">")]
    current = marker_dict
    for part in parts:
        # Resolve the first level directly from the marker dictionary.
        if current is marker_dict and part in current:
            current = current[part]
        # Resolve deeper levels through nested subtypes.
        elif isinstance(current, dict) and "subtypes" in current and part in current["subtypes"]:
            current = current["subtypes"][part]
        else:
            return None
    return current


def get_marker_genes(node):
    """Return marker_genes from a node, or an empty list when absent."""
    if isinstance(node, dict) and "marker_genes" in node:
        return node["marker_genes"]
    return []


def get_subtypes(node):
    """Return the subtype dictionary from a node, or an empty dictionary when absent."""
    if isinstance(node, dict) and "subtypes" in node:
        return node["subtypes"]
    return {}


def extract_L1_markers(marker_dict):
    """
    Collect markers for every first-level cell type.
    Return a mapping from cell-type name to gene list.
    """
    result = {}
    for name, node in marker_dict.items():
        genes = get_marker_genes(node)
        if genes:
            result[name] = genes
    return result


def get_sub_markers_for_parent(marker_dict, parent_label):
    """
    Collect subtype markers for a parent label.
    Nested label paths use the ' > ' separator.
    Return a subtype-to-gene mapping, or an empty dictionary.
    """
    node = navigate_to_node(marker_dict, parent_label)
    if node is None:
        return {}
    
    subtypes = get_subtypes(node)
    result = {}
    for name, subnode in subtypes.items():
        genes = get_marker_genes(subnode)
        if genes:
            result[name] = genes
    return result


def _get_max_depth(marker_dict):
    """Return the maximum depth of the nested marker structure."""
    def _depth(node, d=1):
        subtypes = get_subtypes(node)
        if not subtypes:
            return d
        return max(_depth(v, d + 1) for v in subtypes.values())
    
    return max(_depth(v) for v in marker_dict.values())


# ============================================================
# Iterative annotation workflow
# ============================================================

adata.obs["L1"] = pd.Series(["NA"] * adata.n_obs, index=adata.obs.index, dtype="object")
adata.obs["score1"] = np.nan

all_clusters = adata.obs[cluster_key].unique()
print(f"Detected {len(all_clusters)} clusters")

current_level = 1
max_marker_depth = _get_max_depth(marker_dict)
print(f"Maximum marker depth: {max_marker_depth}")

while current_level <= max_marker_depth:
    print(f"\n{'='*60}")
    print(f"  Annotating level {current_level}")
    print(f"{'='*60}")
    
    # ----------------------------------------------------------
    # Level 1: annotate all clusters.
    # ----------------------------------------------------------
    if current_level == 1:
        L1_markers = extract_L1_markers(marker_dict)
        if not L1_markers:
            print("  Level 1 marker set is empty; annotation cannot continue")
            break
        
        print(f"  Level 1 types: {list(L1_markers.keys())}")
        print(f"  Level 1 marker counts: { {k: len(v) for k, v in L1_markers.items()} }")
        
        result = snap.annotate(adata, L1_markers, group_name=cluster_key, layer=layer)
        
        for _, row in result.iterrows():
            clust = row[cluster_key]
            mask = adata.obs[cluster_key] == clust
            adata.obs.loc[mask, "L1"] = row["class"]
            adata.obs.loc[mask, "score1"] = row["score"]
        
        print("  Level 1 annotation complete:")
        for _, row in result.iterrows():
            print(f"    Cluster {row[cluster_key]} -> {row['class']} (score={row['score']:.3f})")
        
        current_level = 2
        continue
    
    # ----------------------------------------------------------
    # Level 2 and deeper: refine annotations.
    # ----------------------------------------------------------
    prev_col = f"L{current_level - 1}"
    score_col = f"score{current_level - 1}"
    cur_col = f"L{current_level}"
    cur_score_col = f"score{current_level}"
    
    adata.obs[cur_col] = pd.Series(["NA"] * adata.n_obs, index=adata.obs.index, dtype="object")
    adata.obs[cur_score_col] = np.nan
    
    # Group clusters by parent label.
    prev_annotations = adata.obs.groupby(cluster_key)[prev_col].first()
    prev_scores = adata.obs.groupby(cluster_key)[score_col].first()
    
    # Collect groups that require refinement.
    label_to_clusters = {}
    for clust in all_clusters:
        lbl = prev_annotations.get(clust, "NA")
        scr = prev_scores.get(clust, 0.0)
        if lbl == "NA" or pd.isna(lbl):
            continue
        if scr < score_thr:
            continue
        label_to_clusters.setdefault(lbl, []).append(clust)
    
    if not label_to_clusters:
        print("  No clusters require refinement; stopping")
        break
    
    any_annotated = False
    
    for parent_label, cluster_list in label_to_clusters.items():
        n_clusters = len(cluster_list)
        
        # Check whether subtype markers are available.
        sub_markers = get_sub_markers_for_parent(marker_dict, parent_label)
        
        if n_clusters == 1:
            # Determine whether a unique label requires deeper annotation.
            if not sub_markers:
                print(f"  '{parent_label}' has one cluster and no subtypes; retaining the parent label")
                for clust in cluster_list:
                    mask = adata.obs[cluster_key] == clust
                    adata.obs.loc[mask, cur_col] = parent_label
                    adata.obs.loc[mask, cur_score_col] = prev_scores[clust]
            else:
                # A unique label with subtypes is refined even when it has one cluster.
                print(f"  '{parent_label}' has one cluster and subtype markers; refining")
                # Continue into the refinement logic below.
                pass  # Continue with the refinement logic below.
        
        if n_clusters > 1 or (n_clusters == 1 and sub_markers):
            if not sub_markers:
                print(f"  '{parent_label}' has {n_clusters} clusters but no subtype markers; retaining the parent label")
                for clust in cluster_list:
                    mask = adata.obs[cluster_key] == clust
                    adata.obs.loc[mask, cur_col] = parent_label
                    adata.obs.loc[mask, cur_score_col] = prev_scores[clust]
                continue
            
            print(f"  Refining '{parent_label}' ({n_clusters} clusters: {cluster_list})")
            print(f"    Subtypes: {list(sub_markers.keys())}")
            print(f"    Subtype marker counts: { {k: len(v) for k, v in sub_markers.items()} }")
            
            # Subset the AnnData object.
            sub_mask = adata.obs[cluster_key].isin(cluster_list)
            sub_adata = adata[sub_mask].copy()
            
            try:
                result = snap.annotate(sub_adata, sub_markers, group_name=cluster_key, layer=layer)
                
                for _, row in result.iterrows():
                    clust = row[cluster_key]
                    mask = adata.obs[cluster_key] == clust
                    combined_label = f"{parent_label} > {row['class']}"
                    adata.obs.loc[mask, cur_col] = combined_label
                    adata.obs.loc[mask, cur_score_col] = row["score"]
                    any_annotated = True
                    print(f"    Cluster {clust} -> {combined_label} (score={row['score']:.3f})")
            
            except Exception as e:
                print(f"    Refinement failed for '{parent_label}': {e}")
                for clust in cluster_list:
                    mask = adata.obs[cluster_key] == clust
                    adata.obs.loc[mask, cur_col] = parent_label
                    adata.obs.loc[mask, cur_score_col] = prev_scores[clust]
    
    # Cells not covered at this level retain their previous label.
    na_mask = adata.obs[cur_col] == "NA"
    if na_mask.any():
        print(f"  Retained previous labels for {na_mask.sum()} cells")
        adata.obs.loc[na_mask, cur_col] = adata.obs.loc[na_mask, prev_col]
        adata.obs.loc[na_mask, cur_score_col] = adata.obs.loc[na_mask, score_col]
    
    if not any_annotated:
        print("  No new refinements were produced; stopping")
        break
    
    # Check label uniqueness and annotation scores.
    current_annotations = adata.obs.groupby(cluster_key)[cur_col].first()
    current_scores = adata.obs.groupby(cluster_key)[cur_score_col].first()
    label_counts = current_annotations.value_counts()
    shared_labels = label_counts[label_counts > 1]
    
    print(f"\n  Current-level annotation results:")
    for clust in all_clusters:
        lbl = current_annotations.get(clust, "NA")
        scr = current_scores.get(clust, np.nan)
        print(f"    Cluster {clust} -> {lbl} (score={scr:.3f})" if not pd.isna(scr) else f"    Cluster {clust} -> {lbl}")
    
    if len(shared_labels) > 0:
        print(f"  Shared labels that still require refinement: {dict(shared_labels)}")
    else:
        all_scores_ok = all(current_scores.dropna() >= score_thr)
        if all_scores_ok:
            print("  All clusters have unique labels and sufficient scores; annotation is complete")
            break
        else:
            print(f"  All clusters are unique, but some scores are below {score_thr}; stopping")
            # Lower-scoring clusters are not refined further under the configured rule.
            # Scores below the threshold terminate refinement.
            break
    
    current_level += 1


# ============================================================
# Build final annotations.
# ============================================================
level_cols = [f"L{i}" for i in range(1, current_level + 1) if f"L{i}" in adata.obs.columns]
score_cols = [f"score{i}" for i in range(1, current_level + 1) if f"score{i}" in adata.obs.columns]

print(f"\nAnnotation-level columns: {level_cols}")

final_labels = []
final_scores = []
final_depths = []

for idx in adata.obs.index:
    label = "NA"
    score = np.nan
    depth = 0
    for i, (lc, sc) in enumerate(zip(level_cols, score_cols)):
        val = adata.obs.loc[idx, lc]
        if val != "NA" and not pd.isna(val):
            label = val
            score = adata.obs.loc[idx, sc]
            depth = i + 1
    final_labels.append(label)
    final_scores.append(score)
    final_depths.append(depth)

adata.obs["final_annotation"] = pd.Categorical(final_labels)
adata.obs["final_score"] = final_scores
adata.obs["annotation_depth"] = final_depths

# ============================================================
# Report results.
# ============================================================
print("\n" + "=" * 60)
print("  Final annotation results")
print("=" * 60)

summary = (
    adata.obs
    .groupby(cluster_key)
    .agg(
        final_label=("final_annotation", "first"),
        final_score=("final_score", "first"),
        depth=("annotation_depth", "first"),
        n_cells=("final_annotation", "count"),
    )
    .sort_values("final_label")
)
print(summary.to_string())