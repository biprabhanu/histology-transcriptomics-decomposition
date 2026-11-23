"""
Step 04: per-gene decomposition of an H&E -> expression model's score into
noise, cell-type composition and morphology.

Usage:
    python 04_decompose.py resnet50            # all samples in config that have spots + features
    python 04_decompose.py phikon2 --samples TENX191 TENX199

Folds are leave-one-patient-out over the patients that are present. Ridge alpha
is chosen per target by an inner leave-one-patient-out on the training patients
(fallback when only one training patient: leave-one-sample-out, else 5-fold spots).

Output: results/core/decomposition_<model>.csv and results/core/summary_<model>.txt
"""
import os
import argparse
import numpy as np
import pandas as pd

import config


# ---------------------------------------------------------------- loading

def load_sample(sample_id, model):
    """Load one sample's spots and image features, joined by barcode and filtered."""
    spots = np.load(os.path.join(config.SPOT_DIR, f"{sample_id}.npz"), allow_pickle=True)
    feats = np.load(os.path.join(config.FEATURE_DIR, model, f"{sample_id}.npz"), allow_pickle=True)

    # join spots and feature rows by barcode
    spot_bc = spots["barcodes"].astype(str)
    feat_row = {bc: i for i, bc in enumerate(feats["barcodes"].astype(str))}
    in_both = np.array([bc in feat_row for bc in spot_bc])
    feat_idx = np.array([feat_row[bc] for bc in spot_bc[in_both]], dtype=int)

    b1 = spots["counts_b1"][in_both].astype(np.float64)
    b2 = spots["counts_b2"][in_both].astype(np.float64)
    n_cells = spots["n_cells"][in_both]

    # keep spots with tissue and enough target counts
    keep = (n_cells >= config.MIN_CELLS_PER_SPOT) & ((b1 + b2).sum(1) >= config.MIN_COUNTS_PER_SPOT)
    return {
        "sample": sample_id,
        "patient": config.SAMPLE_TO_PATIENT[sample_id],
        "genes": spots["genes"].astype(str),
        "b1": b1[keep], "b2": b2[keep],
        "fractions": spots["fractions"][in_both][keep],
        "n_cells": n_cells[keep],
        "image": feats["features"][feat_idx][keep].astype(np.float64),
        "n_joined": int(in_both.sum()), "n_spots": int(spot_bc.size),
    }


def zscore_columns(x):
    """z-score each column; constant columns become 0."""
    sd = x.std(0)
    sd[sd == 0] = 1.0
    return (x - x.mean(0)) / sd


def log_normalise(counts):
    """log1p(C / total per spot * 100)."""
    total = np.maximum(counts.sum(1, keepdims=True), 1.0)
    return np.log1p(counts / total * 100.0)


def build_targets(counts, genes, gene_idx):
    """Gene targets (log-normalised) + programme scores, for one sample and one count matrix."""
    expr = log_normalise(counts)[:, gene_idx]          # library size uses all shared genes
    kept_genes = list(genes[gene_idx])
    z = zscore_columns(expr)                           # within-sample z-score for programmes
    prog_cols, prog_names = [], []
    for name, members in config.PROGRAMMES.items():
        present = [kept_genes.index(g) for g in members if g in kept_genes]
        if len(present) >= config.MIN_GENES_PER_PROGRAMME:
            prog_cols.append(z[:, present].mean(1))
            prog_names.append(f"PROG_{name}")
    if prog_cols:
        expr = np.hstack([expr, np.column_stack(prog_cols)])
    return expr, kept_genes + prog_names


# ---------------------------------------------------------------- ridge

def ridge_predict_all_alphas(x_train, y_train, x_test, alphas):
    """SVD ridge (no intercept; data are centred). Returns predictions (n_alpha, n_test, n_targets)."""
    u, s, vt = np.linalg.svd(x_train, full_matrices=False)
    uty = u.T @ y_train                    # (k, t)
    xv = x_test @ vt.T                     # (m, k)
    preds = np.empty((len(alphas), x_test.shape[0], y_train.shape[1]))
    for a, alpha in enumerate(alphas):
        d = s / (s ** 2 + alpha)
        preds[a] = xv @ (d[:, None] * uty)
    return preds


def ridge_fit_predict(x_train, y_train, x_test, alpha_per_target, alphas):
    """Fit ridge with one alpha per target; return train fit and test prediction."""
    u, s, vt = np.linalg.svd(x_train, full_matrices=False)
    uty = u.T @ y_train
    fit = np.empty_like(y_train)
    pred = np.empty((x_test.shape[0], y_train.shape[1]))
    xv_train, xv_test = x_train @ vt.T, x_test @ vt.T
    for a, alpha in enumerate(alphas):
        cols = np.where(alpha_per_target == a)[0]
        if cols.size == 0:
            continue
        w = (s / (s ** 2 + alpha))[:, None] * uty[:, cols]     # coefficients in SVD basis
        fit[:, cols] = xv_train @ w
        pred[:, cols] = xv_test @ w
    return fit, pred


def inner_groups(patients, samples, seed=0):
    """Groups for inner alpha selection: patients, else samples, else 5 random spot folds."""
    if len(np.unique(patients)) >= 2:
        return patients
    if len(np.unique(samples)) >= 2:
        print("  WARNING: 1 training patient -> inner leave-one-sample-out", flush=True)
        return samples
    print("  WARNING: 1 training sample -> inner 5-fold spot split", flush=True)
    return np.random.default_rng(seed).integers(0, 5, size=len(samples)).astype(str)


def choose_alpha(x, y, groups, alphas):
    """Index of the alpha with the lowest pooled inner-CV SSE, per target."""
    sse = np.zeros((len(alphas), y.shape[1]))
    for g in np.unique(groups):
        test = groups == g
        preds = ridge_predict_all_alphas(x[~test], y[~test], x[test], alphas)
        sse += ((preds - y[test][None]) ** 2).sum(1)
    return sse.argmin(0)


def fit_ridge(x_tr, y_tr, x_te, groups, alphas):
    """Inner-CV alpha per target, then fit on all training rows."""
    best = choose_alpha(x_tr, y_tr, groups, alphas)
    return ridge_fit_predict(x_tr, y_tr, x_te, best, alphas)


# ---------------------------------------------------------------- metrics

def r2(y, pred):
    """R^2 per column; y is centred within sample so SST = sum(y^2) around the sample mean."""
    sst = ((y - y.mean(0)) ** 2).sum(0)
    return 1 - ((y - pred) ** 2).sum(0) / np.where(sst > 0, sst, np.nan)


def pearson(a, b):
    a = a - a.mean(0)
    b = b - b.mean(0)
    den = np.sqrt((a ** 2).sum(0) * (b ** 2).sum(0))
    return (a * b).sum(0) / np.where(den > 0, den, np.nan)


# ---------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model", help="feature folder name, e.g. resnet50 or phikon2")
    parser.add_argument("--samples", nargs="*", default=None,
                        help="sample ids (default: all in config with spots + features)")
    parser.add_argument("--save-predictions", default=None, metavar="DIR",
                        help="opt-in: save per held-out sample predictions/metrics as DIR/<sample>.npz")
    args = parser.parse_args()
    model, alphas = args.model, config.RIDGE_ALPHAS

    # ---- 1. which samples exist, load and join
    candidates = args.samples or list(config.SAMPLE_TO_PATIENT)
    sample_ids = [s for s in candidates
                  if os.path.exists(os.path.join(config.SPOT_DIR, f"{s}.npz"))
                  and os.path.exists(os.path.join(config.FEATURE_DIR, model, f"{s}.npz"))]
    data = [load_sample(s, model) for s in sample_ids]
    for d in data:
        print(f"{d['sample']} ({d['patient']}): {d['n_joined']}/{d['n_spots']} spots joined, "
              f"{len(d['b1'])} kept", flush=True)
    patients = sorted({d["patient"] for d in data})
    print(f"{len(data)} samples, {len(patients)} patients: {patients}", flush=True)
    if len(patients) < 2:
        raise SystemExit("need at least 2 patients")

    # ---- 2. genes: shared across samples, then detected (B > 0) in >= MIN_SPOT_DETECTION of spots
    shared = [g for g in data[0]["genes"] if all(g in set(d["genes"]) for d in data)]
    shared = [g for g in shared if not g.startswith(config.CONTROL_PREFIXES)]
    for d in data:   # reorder every sample to the shared gene list
        col = {g: i for i, g in enumerate(d["genes"])}
        idx = [col[g] for g in shared]
        d["b1"], d["b2"], d["genes"] = d["b1"][:, idx], d["b2"][:, idx], np.array(shared)
    b_all = np.vstack([d["b1"] + d["b2"] for d in data])
    detected = (b_all > 0).mean(0) >= config.MIN_SPOT_DETECTION
    gene_idx = np.where(detected)[0]
    print(f"genes: {len(shared)} shared, {len(gene_idx)} detected in >= "
          f"{config.MIN_SPOT_DETECTION:.0%} of spots", flush=True)

    # ---- 3-5. targets, ceiling inputs and features per sample
    for d in data:
        y_b, names = build_targets(d["b1"] + d["b2"], d["genes"], gene_idx)
        y_b1, _ = build_targets(d["b1"], d["genes"], gene_idx)
        y_b2, _ = build_targets(d["b2"], d["genes"], gene_idx)
        # noise ceiling: Spearman-Brown of corr(B1, B2), clipped to [0, 1]
        r = pearson(y_b1, y_b2)
        d["ceiling"] = np.clip(np.nan_to_num(2 * r / (1 + r)), 0, 1)
        d["y"] = y_b - y_b.mean(0)                                   # centre B per sample
        comp = np.column_stack([d["fractions"], np.log1p(d["n_cells"])])
        d["comp"] = zscore_columns(comp)
        d["img"] = zscore_columns(d["image"])
    target_names = names
    n_genes = len(gene_idx)

    # mean raw B count per spot (programmes: mean over their member genes)
    gene_mean = b_all[:, gene_idx].mean(0)
    kept = list(np.array(shared)[gene_idx])
    mean_count = list(gene_mean)
    for name in target_names[n_genes:]:
        members = [kept.index(g) for g in config.PROGRAMMES[name[5:]] if g in kept]
        mean_count.append(gene_mean[members].mean())

    # stack everything
    y = np.vstack([d["y"] for d in data])
    x_comp = np.vstack([d["comp"] for d in data])
    x_img = np.vstack([d["img"] for d in data])
    pat = np.concatenate([[d["patient"]] * len(d["y"]) for d in data])
    smp = np.concatenate([[d["sample"]] * len(d["y"]) for d in data])
    run_shuffle = model == "resnet50"   # control run once, resnet50 only

    # ---- 6-7. leave-one-patient-out folds
    rng = np.random.default_rng(config.SEED)
    per_sample = []   # metric dicts per held-out sample
    for held in patients:
        print(f"fold: hold out {held}", flush=True)
        tr, te = pat != held, pat == held
        groups = inner_groups(pat[tr], smp[tr])
        y_tr = y[tr]

        # comp: ridge on composition
        comp_fit, comp_pred = fit_ridge(x_comp[tr], y_tr, x_comp[te], groups, alphas)
        # img: ridge on image features
        _, img_pred = fit_ridge(x_img[tr], y_tr, x_img[te], groups, alphas)
        # full: image ridge on the composition residuals, added to the comp prediction
        _, res_pred = fit_ridge(x_img[tr], y_tr - comp_fit, x_img[te], groups, alphas)
        full_pred = comp_pred + res_pred
        # img_shuffled: image rows permuted within each training sample (control)
        if run_shuffle:
            x_shuf = x_img[tr].copy()
            for s in np.unique(smp[tr]):
                rows = np.where(smp[tr] == s)[0]
                x_shuf[rows] = x_shuf[rng.permutation(rows)]
            _, shuf_pred = fit_ridge(x_shuf, y_tr, x_img[te], groups, alphas)

        # ---- 8. metrics per held-out sample
        te_samples = smp[te]
        for s in np.unique(te_samples):
            m = te_samples == s
            yt = y[te][m]
            sse_comp = ((yt - comp_pred[m]) ** 2).sum(0)
            sse_full = ((yt - full_pred[m]) ** 2).sum(0)
            per_sample.append({
                "r2_comp": r2(yt, comp_pred[m]),
                "r2_img": r2(yt, img_pred[m]),
                "r2_full": r2(yt, full_pred[m]),
                "naive_r": pearson(yt, img_pred[m]),
                "partial_r2": (sse_comp - sse_full) / np.where(sse_comp > 0, sse_comp, np.nan),
                "r2_img_shuffled": r2(yt, shuf_pred[m]) if run_shuffle else np.full(yt.shape[1], np.nan),
                "ceiling": next(d["ceiling"] for d in data if d["sample"] == s),
            })
            if args.save_predictions:
                os.makedirs(args.save_predictions, exist_ok=True)
                p = per_sample[-1]
                np.savez(os.path.join(args.save_predictions, f"{s}.npz"),
                         sample=s, patient=held, target_names=np.array(target_names),
                         y_true=yt, comp_pred=comp_pred[m], img_pred=img_pred[m], full_pred=full_pred[m],
                         r2_comp=p["r2_comp"], r2_img=p["r2_img"], r2_full=p["r2_full"],
                         ceiling=p["ceiling"])

    # average over held-out samples
    keys = per_sample[0].keys()
    avg = {k: np.nanmean(np.vstack([p[k] for p in per_sample]), 0) for k in keys}

    # ---- 9. results table
    df = pd.DataFrame({
        "target": target_names,
        "is_programme": [t.startswith("PROG_") for t in target_names],
        "mean_count": mean_count,
        "ceiling": avg["ceiling"],
        "naive_r": avg["naive_r"],
        "r2_img": avg["r2_img"],
        "r2_comp": avg["r2_comp"],
        "r2_full": avg["r2_full"],
        "delta_morph": avg["r2_full"] - avg["r2_comp"],
        "partial_r2": avg["partial_r2"],
        "r2_img_shuffled": avg["r2_img_shuffled"],
    })
    ceil_safe = df["ceiling"].where(df["ceiling"] >= 0.2)   # shares only where ceiling >= 0.2 (rest = noise_limited)
    df["comp_share"] = df["r2_comp"] / ceil_safe
    df["morph_share"] = df["delta_morph"] / ceil_safe
    df["category"] = np.select(
        [df["ceiling"] < 0.2,
         df["morph_share"] >= 0.10,
         df["r2_img"] >= 0.10 * df["ceiling"]],
        ["noise_limited", "morphology_visible", "composition_explained"],
        default="not_predicted")

    os.makedirs(config.RESULT_DIR, exist_ok=True)
    csv_path = os.path.join(config.RESULT_DIR, f"decomposition_{model}.csv")
    df.to_csv(csv_path, index=False)

    # ---- summary
    genes_df = df[~df["is_programme"]]
    success = genes_df[genes_df["r2_img"] > 0.05]
    lines = [
        f"model: {model}   samples: {len(data)}   patients: {len(patients)}   held-out samples: {len(per_sample)}",
        f"genes: {len(genes_df)}",
        f"median ceiling: {genes_df['ceiling'].median():.3f}",
        "category counts (genes): " + ", ".join(f"{k}={v}" for k, v in genes_df["category"].value_counts().items()),
        f"median r2_comp / r2_img among genes with r2_img > 0.05 (n={len(success)}): "
        f"{(success['r2_comp'] / success['r2_img']).median():.3f}",
        "",
        "top 15 genes by morph_share:",
        genes_df.nlargest(15, "morph_share")[["target", "ceiling", "r2_comp", "r2_img", "r2_full", "morph_share",
                                               "category"]].to_string(index=False, float_format="%.3f"),
        "",
        "programmes:",
        df[df["is_programme"]][["target", "ceiling", "r2_comp", "r2_img", "r2_full", "delta_morph",
                                "comp_share", "morph_share", "category"]].to_string(index=False, float_format="%.3f"),
    ]
    with open(os.path.join(config.RESULT_DIR, f"summary_{model}.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines[:5]), flush=True)
    print(f"wrote {csv_path}", flush=True)


if __name__ == "__main__":
    main()
