"""
Step 02b (robustness) - cell types from ONE joint clustering of all samples.

Same part-A counts as step 02 (its functions are imported, same seed and thinning),
but instead of a Leiden run per sample, all cells are pooled:
  subsample <= 100k cells per sample -> normalise, log1p, scale -> PCA (30)
  -> Harmony by sample (all cells) -> Leiden 1.0 on the subsample
  -> same marker rule as 02 (best type, "Other" if score < 0.3)
  -> every cell labelled by kNN (k=15) transfer in the harmonised PCA space.
Writes data/spots_joint/ID.npz in the same format as data/spots/ID.npz;
only `fractions` change (checked: barcodes, centres, genes, B counts, n_cells identical).

The per-sample labels of step 02 are recomputed (same function, same seed) so the
old and new label of every cell can be compared; we check they reproduce data/spots.

Usage: python 02b_joint_celltypes.py
"""
import os
import sys
import time
import importlib.util

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse
from sklearn.neighbors import KNeighborsClassifier

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import config

spec = importlib.util.spec_from_file_location("build_spots", os.path.join(HERE, "02_build_spots.py"))
s02 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s02)
log = s02.log

SAMPLES = ["TENX191", "TENX195", "TENX199", "TENX202"]
MAX_CELLS_CLUSTER = 100_000
OUT_DIR = os.path.join(config.PROJECT, "data", "spots_joint")
CACHE_DIR = os.path.join(OUT_DIR, "_cache")
TYPE_NAMES = list(config.CELL_TYPE_MARKERS) + ["Other"]


def part_a_cells(sample_id):
    """Exactly steps 1-3 of 02.main: part-A cell x gene counts, B counts, cell -> spot."""
    rng = np.random.default_rng(config.SEED)
    centres = s02.load_spot_grid(sample_id)
    half_width_px = config.SPOT_SIZE_UM / 2 / s02.load_pixel_size(sample_id)
    tx, gene_names, cell_names = s02.load_transcripts(sample_id)
    n_tx, n_genes, n_spots = len(tx["gene"]), len(gene_names), len(centres)

    u = rng.random(n_tx)
    part = np.where(u < config.THIN_A, 0, np.where(u < config.THIN_A + config.THIN_B1, 1, 2))
    tx_spot = s02.assign_to_spots(tx["x"], tx["y"], centres, half_width_px)
    counts_b1 = s02.count_matrix(np.where(part == 1, tx_spot, -1), tx["gene"], n_spots, n_genes)
    counts_b2 = s02.count_matrix(np.where(part == 2, tx_spot, -1), tx["gene"], n_spots, n_genes)

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
    good = np.flatnonzero(a_per_cell >= config.MIN_TRANSCRIPTS_PER_CELL)
    cell_spot = s02.assign_to_spots(cx[good] / n_per_cell[good], cy[good] / n_per_cell[good],
                                    centres, half_width_px)
    return dict(centres=centres, genes=gene_names, counts_b1=counts_b1, counts_b2=counts_b2,
                X=cell_by_gene[good], cell_spot=cell_spot)


def fractions_from_labels(cell_spot, labels, n_spots):
    type_index = pd.Categorical(labels, categories=TYPE_NAMES).codes
    type_counts = s02.count_matrix(cell_spot, type_index, n_spots, len(TYPE_NAMES))
    n_cells = type_counts.sum(axis=1)
    return type_counts / np.maximum(n_cells, 1)[:, None], n_cells


def stage1(sample_id):
    """Part-A cells + the original per-sample labels (cached)."""
    path = os.path.join(CACHE_DIR, f"{sample_id}.npz")
    if os.path.exists(path):
        z = np.load(path, allow_pickle=True)
        return dict(X=sparse.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])),
                    genes=z["genes"], cell_spot=z["cell_spot"], old=z["old"])
    log(f"{sample_id}: rebuilding part A")
    d = part_a_cells(sample_id)
    ref = np.load(os.path.join(config.SPOT_DIR, f"{sample_id}.npz"), allow_pickle=True)
    assert np.array_equal(ref["counts_b1"], d["counts_b1"]) and np.array_equal(ref["counts_b2"], d["counts_b2"]), \
        f"{sample_id}: B counts differ from data/spots"
    assert np.array_equal(ref["barcodes"], d["centres"].index.to_numpy())
    log(f"{sample_id}: {d['X'].shape[0]:,} cells; B counts identical to data/spots; per-sample Leiden")
    old, _ = s02.call_cell_types(d["X"], d["genes"])
    frac_old, n_cells = fractions_from_labels(d["cell_spot"], old, len(d["centres"]))
    ok_frac = np.allclose(frac_old, ref["fractions"]) and np.array_equal(n_cells, ref["n_cells"])
    log(f"{sample_id}: recomputed per-sample labels reproduce data/spots fractions: {ok_frac}")
    X = d["X"].tocsr()
    np.savez(path, data=X.data, indices=X.indices, indptr=X.indptr, shape=np.array(X.shape),
             genes=d["genes"], cell_spot=d["cell_spot"], old=old)
    return dict(X=X, genes=d["genes"], cell_spot=d["cell_spot"], old=old)


def harmonise(pcs, batch):
    try:
        import harmonypy
    except ImportError:
        log("harmonypy not available -> plain joint PCA")
        return pcs, "pca"
    meta = pd.DataFrame({"sample": batch})
    ho = harmonypy.run_harmony(pcs, meta, ["sample"], random_state=config.SEED)
    z = np.asarray(ho.Z_corr)
    if z.shape[0] != pcs.shape[0]:
        z = z.T
    return z, "harmony"


def main():
    os.makedirs(CACHE_DIR, exist_ok=True)
    data = {s: stage1(s) for s in SAMPLES}

    genes = list(data[SAMPLES[0]]["genes"])
    for s in SAMPLES:   # same panel, possibly different order -> reorder columns by name
        assert set(data[s]["genes"]) == set(genes), "gene panels differ"
        col = {g: i for i, g in enumerate(data[s]["genes"])}
        data[s]["X"] = data[s]["X"][:, [col[g] for g in genes]]
    X = sparse.vstack([data[s]["X"] for s in SAMPLES]).tocsr()
    batch = np.concatenate([[s] * data[s]["X"].shape[0] for s in SAMPLES])
    old = np.concatenate([data[s]["old"] for s in SAMPLES])
    n = X.shape[0]
    log(f"pooled {n:,} cells x {len(genes)} genes")

    # normalise + log1p (per cell), subsample for scaling / PCA / clustering
    lib = np.asarray(X.sum(1)).ravel()
    L = np.log1p(np.asarray(X.multiply(100.0 / lib[:, None]).todense(), dtype=np.float32))
    rng = np.random.default_rng(config.SEED)
    sub = np.concatenate([rng.choice(np.flatnonzero(batch == s),
                                     min(MAX_CELLS_CLUSTER, (batch == s).sum()), replace=False)
                          for s in SAMPLES])
    sub.sort()
    mu, sd = L[sub].mean(0), L[sub].std(0, ddof=1)
    sd[sd == 0] = 1
    S = np.clip((L - mu) / sd, None, 10).astype(np.float32)
    del L

    cells = ad.AnnData(S[sub])
    cells.var_names = genes
    sc.tl.pca(cells, n_comps=30, random_state=config.SEED)
    pcs = (S - S[sub].mean(0)) @ cells.varm["PCs"]           # all cells, same loadings
    Z, method = harmonise(pcs, batch)
    log(f"PCA + {method} done; clustering {len(sub):,} subsampled cells")

    cells.obsm["X_emb"] = Z[sub]
    sc.pp.neighbors(cells, n_neighbors=15, use_rep="X_emb", random_state=config.SEED, transformer="sklearn")
    sc.tl.leiden(cells, resolution=config.LEIDEN_RESOLUTION, flavor="igraph",
                 n_iterations=2, random_state=config.SEED)
    clusters = cells.obs["leiden"].astype(int).to_numpy()

    # same marker rule as 02.call_cell_types
    scaled = pd.DataFrame(S[sub], columns=genes)
    cluster_means = scaled.groupby(clusters).mean()
    cluster_z = (cluster_means - cluster_means.mean()) / (cluster_means.std() + 1e-8)
    scores = pd.DataFrame(index=cluster_z.index)
    for cell_type, markers in config.CELL_TYPE_MARKERS.items():
        present = [m for m in markers if m in cluster_z.columns]
        scores[cell_type] = cluster_z[present].mean(axis=1) if present else -np.inf
    best = scores.idxmax(axis=1)
    best[scores.max(axis=1) < 0.3] = "Other"
    comp = pd.crosstab(clusters, batch[sub], normalize="index").round(2)
    summary = pd.DataFrame({"cluster": scores.index, "cell_type": best.values,
                            "n_cells": np.bincount(clusters)[scores.index],
                            "best_score": scores.max(axis=1).round(2).values}).join(comp, on="cluster")
    summary.to_csv(os.path.join(OUT_DIR, "joint_clusters.csv"), index=False)
    log(f"{len(scores)} joint clusters\n{summary.to_string(index=False)}")

    # kNN transfer to every cell (subsampled cells keep their own cluster label)
    sub_labels = best.loc[clusters].to_numpy()
    new = np.empty(n, dtype=object)
    new[sub] = sub_labels
    rest = np.setdiff1d(np.arange(n), sub)
    knn = KNeighborsClassifier(n_neighbors=15, n_jobs=-1).fit(Z[sub], sub_labels)
    for i in range(0, len(rest), 50_000):
        new[rest[i:i + 50_000]] = knn.predict(Z[rest[i:i + 50_000]])
    log("kNN transfer done")

    # per-sample outputs + checks + old-vs-new table
    rows, offset = [], 0
    for s in SAMPLES:
        m = slice(offset, offset + data[s]["X"].shape[0])
        offset = m.stop
        ref = np.load(os.path.join(config.SPOT_DIR, f"{s}.npz"), allow_pickle=True)
        frac, n_cells = fractions_from_labels(data[s]["cell_spot"], new[m], len(ref["barcodes"]))
        assert np.array_equal(n_cells, ref["n_cells"]), f"{s}: n_cells differ"
        out = os.path.join(OUT_DIR, f"{s}.npz")
        np.savez_compressed(out, barcodes=ref["barcodes"], centres=ref["centres"], genes=ref["genes"],
                            counts_b1=ref["counts_b1"], counts_b2=ref["counts_b2"],
                            type_names=np.array(TYPE_NAMES), fractions=frac, n_cells=n_cells)
        chk = np.load(out, allow_pickle=True)
        same = all(np.array_equal(chk[k], ref[k]) for k in
                   ["barcodes", "centres", "genes", "counts_b1", "counts_b2", "type_names", "n_cells"])
        same_dtypes = all(chk[k].dtype == ref[k].dtype and chk[k].shape == ref[k].shape for k in ref.files)
        log(f"{s}: wrote {out}; identical non-fraction arrays: {same}; same dtypes/shapes: {same_dtypes}; "
            f"Other {np.mean(data[s]['old'] == 'Other'):.1%} -> {np.mean(new[m] == 'Other'):.1%}")
        rows.append(pd.DataFrame({"sample": s, "old": data[s]["old"], "new": new[m]}))
    tab = pd.concat(rows)
    pd.crosstab([tab["sample"], tab["old"]], tab["new"]).to_csv(os.path.join(OUT_DIR, "old_vs_new.csv"))
    for s, o in [("TENX202", "Mast"), ("TENX191", "Other")]:
        v = tab[(tab["sample"] == s) & (tab["old"] == o)]["new"].value_counts(normalize=True).head(3)
        log(f"{s} old '{o}' -> " + ", ".join(f"{k} {p:.0%}" for k, p in v.items()))


if __name__ == "__main__":
    main()
