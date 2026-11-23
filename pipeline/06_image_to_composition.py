"""
Step 06 (robustness): can the H&E image encoder predict per-spot cell-type composition?

Target  = per-spot cell-type fractions (10 types + Other), centred within sample.
Feature = image features (resnet50 / phikon2), z-scored within sample.
Folds   = leave-one-patient-out; ridge alpha picked per target by inner LOPO (04's functions).
Metric  = R^2 per cell type per held-out sample (sample == patient here), then the mean.
A type that is constant in a held-out sample gives NaN and is kept as NaN
(the mean ignores NaN). Pre-set reading: R^2 >= 0.3 means 'image encodes this type'.

Usage: python 06_image_to_composition.py [resnet50 phikon2]
Output: results/robustness/image_to_composition.csv
"""
import os, sys, importlib.util
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config

# import 04 (filename starts with a digit)
spec = importlib.util.spec_from_file_location("decompose04", os.path.join(os.path.dirname(os.path.abspath(__file__)), "04_decompose.py"))
d04 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d04)

SAMPLES = ["TENX191", "TENX195", "TENX199", "TENX202"]
TYPES = list(config.CELL_TYPE_MARKERS) + ["Other"]   # same order as type_names in the spot npz


def run(model):
    data = [d04.load_sample(s, model) for s in SAMPLES]   # same spot filters as 04
    for d in data:
        d["y"] = d["fractions"] - d["fractions"].mean(0)      # centre within sample
        d["img"] = d04.zscore_columns(d["image"])
    y = np.vstack([d["y"] for d in data]); x = np.vstack([d["img"] for d in data])
    pat = np.concatenate([[d["patient"]] * len(d["y"]) for d in data])
    smp = np.concatenate([[d["sample"]] * len(d["y"]) for d in data])
    raw = np.vstack([d["fractions"] for d in data])
    assert y.shape[1] == len(TYPES)

    per_patient = {}
    for held in sorted(set(pat)):
        tr, te = pat != held, pat == held
        groups = d04.inner_groups(pat[tr], smp[tr])
        _, pred = d04.fit_ridge(x[tr], y[tr], x[te], groups, config.RIDGE_ALPHAS)
        per_patient[held] = d04.r2(y[te], pred)        # NaN where target is constant
    rows = []
    for j, t in enumerate(TYPES):
        r = {p: per_patient[p][j] for p in per_patient}
        rows.append({"encoder": model, "cell_type": t,
                     "r2_mean": np.nanmean(list(r.values())) if not all(np.isnan(list(r.values()))) else np.nan,
                     **{f"r2_{p}": v for p, v in r.items()}, "mean_fraction": raw[:, j].mean()})
    return rows


if __name__ == "__main__":
    models = sys.argv[1:] or ["resnet50", "phikon2"]
    df = pd.DataFrame([r for m in models for r in run(m)])
    out = os.path.join(config.PROJECT, "results", "robustness")
    os.makedirs(out, exist_ok=True)
    df.to_csv(os.path.join(out, "image_to_composition.csv"), index=False)
    print(df.pivot(index="cell_type", columns="encoder", values="r2_mean").round(3).to_string())
