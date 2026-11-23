"""
Step 02 - turn one Xenium section into spot-level data for the decomposition.

For one sample this script:
  1. loads the H&E-aligned transcripts and keeps good-quality gene transcripts;
  2. thins them: each transcript goes at random to part A (50%), B1 (25%) or B2 (25%);
  3. assigns every transcript to the 100 um spot square it falls in;
  4. builds a cell x gene matrix from part A only, clusters cells (Leiden) and
     names each cluster with marker genes -> one cell type per cell;
  5. places each cell in a spot (by centroid) -> spot x cell-type fractions;
  6. saves, per spot: B1 counts, B2 counts, cell-type fractions, number of cells.

Only spots that have an H&E patch (HEST tissue patches) are kept.

Usage: python 02_build_spots.py TENX199
"""
import os
import sys
import json
import time

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pyarrow.compute as pc
import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.spatial import cKDTree

import config


def log(message):
    print(time.strftime("%H:%M:%S"), message, flush=True)


# ---------------------------------------------------------------- loading
def load_spot_grid(sample_id):
    """Spot barcodes + centre positions (H&E pixels) for spots that have an H&E patch."""
    st = ad.read_h5ad(os.path.join(config.HEST_DIR, "st", f"{sample_id}.h5ad"))
    centres = pd.DataFrame(st.obsm["spatial"], index=st.obs_names, columns=["x", "y"])
    with h5py.File(os.path.join(config.HEST_DIR, "patches", f"{sample_id}.h5"), "r") as h5:
        patch_barcodes = [b[0].decode() for b in h5["barcode"][:]]
    centres = centres.loc[patch_barcodes]          # same order as the patch file
    return centres


def load_pixel_size(sample_id):
    with open(os.path.join(config.HEST_DIR, "metadata", f"{sample_id}.json")) as f:
        return float(json.load(f)["pixel_size_um_estimated"])


def load_transcripts(sample_id):
    """Good-quality gene transcripts: gene code, cell code, H&E position."""
    path = os.path.join(config.HEST_DIR, "transcripts", f"{sample_id}_transcripts.parquet")
    table = pq.read_table(path, columns=["feature_name", "cell_id", "qv", "is_gene", "he_x", "he_y"])
    keep = pc.and_(pc.greater_equal(table["qv"], config.MIN_QV), table["is_gene"])
    table = table.filter(keep)

    # Integer codes instead of strings (much faster and lighter)
    gene_codes = table["feature_name"].combine_chunks().dictionary_encode()
    cell_codes = table["cell_id"].combine_chunks().dictionary_encode()
    gene_names = np.array(gene_codes.dictionary.to_pylist())
    cell_names = np.array(cell_codes.dictionary.to_pylist())

    transcripts = {
        "gene": gene_codes.indices.to_numpy(),
        "cell": cell_codes.indices.to_numpy(),
        "x": table["he_x"].to_numpy(),
        "y": table["he_y"].to_numpy(),
    }
    return transcripts, gene_names, cell_names


# ---------------------------------------------------------------- spots
def assign_to_spots(x, y, centres, half_width_px):
    """Index of the spot square containing each point, or -1 if none."""
    tree = cKDTree(centres[["x", "y"]].to_numpy())
    # p=inf gives the square (Chebyshev) distance, i.e. "inside the square"
    distance, index = tree.query(np.column_stack([x, y]), k=1, p=np.inf,
                                 distance_upper_bound=half_width_px, workers=-1)
    index[~np.isfinite(distance)] = -1
    return index


def count_matrix(row_index, gene_index, n_rows, n_genes):
    """Dense row x gene count matrix from one entry per transcript."""
    ok = row_index >= 0
    flat = row_index[ok].astype(np.int64) * n_genes + gene_index[ok]
    counts = np.bincount(flat, minlength=n_rows * n_genes)
    return counts.reshape(n_rows, n_genes).astype(np.int32)


# ---------------------------------------------------------------- cell types
def call_cell_types(cell_by_gene, gene_names):
    """Leiden clusters on part-A counts, each cluster named by its best marker set."""
    cells = ad.AnnData(cell_by_gene.astype(np.float32))
    cells.var_names = gene_names
    sc.pp.normalize_total(cells, target_sum=100)
    sc.pp.log1p(cells)
    sc.pp.scale(cells, max_value=10)
    sc.tl.pca(cells, n_comps=30, random_state=config.SEED)
    sc.pp.neighbors(cells, n_neighbors=15, random_state=config.SEED, transformer="sklearn")  # exact kNN; pynndescent path breaks with umap-learn 0.6
    sc.tl.leiden(cells, resolution=config.LEIDEN_RESOLUTION, flavor="igraph",
                 n_iterations=2, random_state=config.SEED)

    # Mean scaled expression per cluster, then z-score each gene across clusters
    clusters = cells.obs["leiden"].astype(int).to_numpy()
    scaled = pd.DataFrame(cells.X, columns=gene_names)
    cluster_means = scaled.groupby(clusters).mean()
    cluster_z = (cluster_means - cluster_means.mean()) / (cluster_means.std() + 1e-8)

    # Score each cluster for each cell type = mean z of that type's markers on the panel
    scores = pd.DataFrame(index=cluster_z.index)
    for cell_type, markers in config.CELL_TYPE_MARKERS.items():
        present = [m for m in markers if m in cluster_z.columns]
        scores[cell_type] = cluster_z[present].mean(axis=1) if present else -np.inf

    best = scores.idxmax(axis=1)
    best[scores.max(axis=1) < 0.3] = "Other"        # no clear marker signal
    summary = pd.DataFrame({"cluster": scores.index, "cell_type": best.values,
                            "n_cells": np.bincount(clusters)[scores.index],
                            "best_score": scores.max(axis=1).round(2).values})
    return best.loc[clusters].to_numpy(), summary


# ---------------------------------------------------------------- main
def main(sample_id):
    os.makedirs(config.SPOT_DIR, exist_ok=True)
    rng = np.random.default_rng(config.SEED)

    log(f"{sample_id}: loading")
    centres = load_spot_grid(sample_id)
    pixel_size = load_pixel_size(sample_id)
    half_width_px = config.SPOT_SIZE_UM / 2 / pixel_size
    tx, gene_names, cell_names = load_transcripts(sample_id)
    n_tx, n_genes, n_spots = len(tx["gene"]), len(gene_names), len(centres)
    log(f"{sample_id}: {n_tx:,} transcripts, {n_genes} genes, {n_spots} spots, "
        f"spot = {2 * half_width_px:.0f} px")

    # 1. Data thinning: part 0 = A, 1 = B1, 2 = B2
    u = rng.random(n_tx)
    part = np.where(u < config.THIN_A, 0, np.where(u < config.THIN_A + config.THIN_B1, 1, 2))

    # 2. Transcripts -> spots, then B1 / B2 spot counts (all transcripts, in cells or not)
    tx_spot = assign_to_spots(tx["x"], tx["y"], centres, half_width_px)
    log(f"{sample_id}: {np.mean(tx_spot >= 0):.1%} of transcripts fall inside a patch spot")
    counts_b1 = count_matrix(np.where(part == 1, tx_spot, -1), tx["gene"], n_spots, n_genes)
    counts_b2 = count_matrix(np.where(part == 2, tx_spot, -1), tx["gene"], n_spots, n_genes)

    # 3. Cells: part-A counts and centroids (centroid uses every transcript of the cell)
    unassigned = np.flatnonzero(cell_names == "UNASSIGNED")
    in_cell = ~np.isin(tx["cell"], unassigned)
    cell_code = tx["cell"][in_cell]
    n_cell_codes = len(cell_names)
    cx = np.bincount(cell_code, weights=tx["x"][in_cell], minlength=n_cell_codes)
    cy = np.bincount(cell_code, weights=tx["y"][in_cell], minlength=n_cell_codes)
    n_per_cell = np.bincount(cell_code, minlength=n_cell_codes)

    is_a = part[in_cell] == 0
    cell_by_gene = sparse.csr_matrix(
        (np.ones(is_a.sum(), dtype=np.float32), (cell_code[is_a], tx["gene"][in_cell][is_a])),
        shape=(n_cell_codes, n_genes))
    a_per_cell = np.asarray(cell_by_gene.sum(axis=1)).ravel()
    good_cells = np.flatnonzero(a_per_cell >= config.MIN_TRANSCRIPTS_PER_CELL)
    log(f"{sample_id}: {len(good_cells):,} cells with >= {config.MIN_TRANSCRIPTS_PER_CELL} part-A transcripts")

    # 4. Cell types from part A
    cell_types, cluster_summary = call_cell_types(cell_by_gene[good_cells], gene_names)
    cluster_summary.to_csv(os.path.join(config.SPOT_DIR, f"{sample_id}_clusters.csv"), index=False)
    log(f"{sample_id}: cell types\n{pd.Series(cell_types).value_counts().to_string()}")

    # 5. Cells -> spots -> composition fractions
    centroid_x = cx[good_cells] / n_per_cell[good_cells]
    centroid_y = cy[good_cells] / n_per_cell[good_cells]
    cell_spot = assign_to_spots(centroid_x, centroid_y, centres, half_width_px)
    type_names = list(config.CELL_TYPE_MARKERS) + ["Other"]
    type_index = pd.Categorical(cell_types, categories=type_names).codes
    type_counts = count_matrix(cell_spot, type_index, n_spots, len(type_names))
    n_cells = type_counts.sum(axis=1)
    fractions = type_counts / np.maximum(n_cells, 1)[:, None]

    # 6. Save
    out = os.path.join(config.SPOT_DIR, f"{sample_id}.npz")
    np.savez_compressed(out, barcodes=centres.index.to_numpy(), centres=centres.to_numpy(),
                        genes=gene_names, counts_b1=counts_b1, counts_b2=counts_b2,
                        type_names=np.array(type_names), fractions=fractions, n_cells=n_cells)
    log(f"{sample_id}: saved {out}  (spots with >= {config.MIN_CELLS_PER_SPOT} cells: "
        f"{np.sum(n_cells >= config.MIN_CELLS_PER_SPOT)})")


if __name__ == "__main__":
    main(sys.argv[1])
