"""
Step 06b: per-held-out-patient breakdown of the step-04 decomposition, plus a
patient-level bootstrap of the headline numbers.

Input : results/robustness/preds/<model>/<sample>.npz  (from 04 --save-predictions)
Output: results/robustness/per_patient.csv, results/robustness/per_patient_summary.txt

Aggregation mirrors 04: per-sample metric vectors are nanmean-ed over held-out
samples, then per-gene shares / categories / medians are computed on the means.
"""
import os
import glob
import numpy as np
import pandas as pd

import config

ROB = os.path.join(config.PROJECT, "results", "robustness")
MODELS = ["resnet50", "phikon2"]
FOCUS = ["ADIPOQ", "COLEC12", "NPY1R"]
N_BOOT, SEED = 1000, 0
KEYS = ["r2_comp", "r2_img", "r2_full", "ceiling"]


def load(model):
    out = []
    for f in sorted(glob.glob(os.path.join(ROB, "preds", model, "*.npz"))):
        z = np.load(f, allow_pickle=True)
        out.append({"sample": str(z["sample"]), "patient": str(z["patient"]),
                    "names": z["target_names"].astype(str), **{k: z[k] for k in KEYS}})
    return out


def table(samples):
    """Same as 04: mean over samples, then shares and categories. Genes only."""
    avg = {k: np.nanmean(np.vstack([s[k] for s in samples]), 0) for k in KEYS}
    df = pd.DataFrame({"target": samples[0]["names"], **avg})
    df["delta_morph"] = df["r2_full"] - df["r2_comp"]
    ceil_safe = df["ceiling"].where(df["ceiling"] >= 0.2)
    df["morph_share"] = df["delta_morph"] / ceil_safe
    df["category"] = np.select(
        [df["ceiling"] < 0.2, df["morph_share"] >= 0.10, df["r2_img"] >= 0.10 * df["ceiling"]],
        ["noise_limited", "morphology_visible", "composition_explained"], default="not_predicted")
    return df[~df["target"].str.startswith("PROG_")].reset_index(drop=True)


def headline(g):
    succ = g[g["r2_img"] > 0.05]
    return {
        "median_ceiling": g["ceiling"].median(),
        "noise_limited_share": (g["category"] == "noise_limited").mean(),
        "median_comp_over_img": (succ["r2_comp"] / succ["r2_img"]).median(),
        "n_img_gt_005": len(succ),
        "median_delta_morph": g["delta_morph"].median(),
        "n_morph_share_ge_010": int((g["morph_share"] >= 0.10).sum()),
        "n_morph_share_ge_010_r2full_pos": int(((g["morph_share"] >= 0.10) & (g["r2_full"] > 0)).sum()),
        "n_composition_explained": int((g["category"] == "composition_explained").sum()),
    }


def main():
    rows, lines = [], []
    boot_keys = ["median_ceiling", "noise_limited_share", "median_comp_over_img",
                 "median_delta_morph", "n_morph_share_ge_010_r2full_pos"]
    for model in MODELS:
        samples = load(model)
        assert len(samples) == 4, f"{model}: expected 4 held-out samples, found {len(samples)}"
        pooled = headline(table(samples))
        for s in samples:
            g = table([s])
            h = headline(g)
            dm = g.set_index("target")["delta_morph"]
            rows.append({"model": model, "patient": s["patient"], "sample": s["sample"], **h,
                         **{f"delta_morph_{x}": dm.get(x, np.nan) for x in FOCUS},
                         "headline_holds": bool(h["median_comp_over_img"] >= 0.9
                                                and h["n_morph_share_ge_010"] <= 0.05 * len(g))})
        rows.append({"model": model, "patient": "pooled", "sample": "all", **pooled,
                     **{f"delta_morph_{x}": table(samples).set_index("target")["delta_morph"].get(x, np.nan)
                        for x in FOCUS}})

        # patient bootstrap
        rng = np.random.default_rng(SEED)
        boot = []
        for _ in range(N_BOOT):
            idx = rng.integers(0, len(samples), size=len(samples))
            boot.append(headline(table([samples[i] for i in idx])))
        boot = pd.DataFrame(boot)
        lines.append(f"== {model}: patient bootstrap ({N_BOOT} resamples of 4 patients, seed {SEED}) ==")
        for k in boot_keys:
            lo, hi = np.nanpercentile(boot[k], [2.5, 97.5])
            lines.append(f"  {k:34s} point={pooled[k]:.3f}  95% CI [{lo:.3f}, {hi:.3f}]")
        lines.append("")

    df = pd.DataFrame(rows)
    os.makedirs(ROB, exist_ok=True)
    df.to_csv(os.path.join(ROB, "per_patient.csv"), index=False)
    cols = ["model", "patient", "sample", "median_ceiling", "noise_limited_share", "median_comp_over_img",
            "n_img_gt_005", "median_delta_morph", "n_morph_share_ge_010", "n_morph_share_ge_010_r2full_pos"] \
        + [f"delta_morph_{x}" for x in FOCUS] + ["headline_holds"]
    per = ["== per held-out patient (headline_holds: median r2_comp/r2_img >= 0.9 and "
           "genes with morph_share >= 0.10 <= 5% of genes) ==",
           df[cols].to_string(index=False, float_format="%.3f"), ""]
    text = "\n".join(per + lines)
    with open(os.path.join(ROB, "per_patient_summary.txt"), "w") as fh:
        fh.write(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
