"""Pre-specified fix for unstable held-out R² (e.g. NPY1R): floor each held-out sample's R² at 0
before computing shares. R2+ = max(R2, 0) per sample; averaged over held-out samples;
delta_morph+ = R2_full+ - R2_comp+; shares only where ceiling >= 0.2.
Reads the saved per-sample predictions (results/robustness/preds/<model>/*.npz) from 04 --save-predictions.
Writes results/robustness/floored/ (CSV + summary) and floored figures in results/robustness/figures/."""
import os, glob, importlib.util
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.environ.get("DH607_STORAGE", "/path/to/storage") + "/dh607"
ROB = os.path.join(ROOT, "results", "robustness")
OUT, FIG = os.path.join(ROB, "floored"), os.path.join(ROB, "figures")
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
MODELS = ["resnet50", "phikon2"]
rng = np.random.default_rng(0)


def load(model):
    """Per-sample, per-target R² vectors + ceiling, as arrays (n_samples, n_targets)."""
    files = sorted(glob.glob(os.path.join(ROB, "preds", model, "*.npz")))
    zs = [np.load(f, allow_pickle=True) for f in files]
    names = zs[0]["target_names"].astype(str)
    m = {k: np.vstack([z[k] for z in zs]) for k in ["r2_comp", "r2_img", "r2_full", "ceiling"]}
    return names, [str(z["patient"]) for z in zs], m


def table(names, m, idx):
    """Average floored R² over the samples in idx (may repeat, for bootstrap) and derive shares."""
    pos = {k: np.nanmean(np.clip(m[k][idx], 0, None), 0) for k in ["r2_comp", "r2_img", "r2_full"]}
    df = pd.DataFrame({"target": names, "ceiling": np.nanmean(m["ceiling"][idx], 0), **pos})
    df["delta_morph"] = df["r2_full"] - df["r2_comp"]
    ok = df["ceiling"].where(df["ceiling"] >= 0.2)
    df["comp_share"], df["morph_share"] = df["r2_comp"] / ok, df["delta_morph"] / ok
    df["is_programme"] = df["target"].str.startswith("PROG_")
    return df


def headline(df):
    g = df[~df["is_programme"]]
    s = g[g["r2_img"] > 0.05]
    return {"median_comp_over_img": float(np.median(s["r2_comp"] / s["r2_img"])) if len(s) else np.nan,
            "n_img_gt_005": len(s),
            "median_delta_morph": float(g["delta_morph"].median()),
            "n_morph_genes": int(((g["morph_share"] >= 0.10) & (g["r2_full"] > 0)).sum())}


lines, rows, rep_rows = [], [], []
for model in MODELS:
    names, pats, m = load(model)
    n = len(pats)
    full = table(names, m, np.arange(n))
    full.to_csv(os.path.join(OUT, f"decomposition_floored_{model}.csv"), index=False)
    h = headline(full)
    boot = pd.DataFrame([headline(table(names, m, rng.integers(0, n, n))) for _ in range(1000)])
    # genes with morph_share >= 0.10 and R2_full > 0 in >= 3 of 4 held-out patients
    per = [table(names, m, np.array([i])) for i in range(n)]
    hits = sum((((p["morph_share"] >= 0.10) & (p["r2_full"] > 0) & ~p["is_programme"]).astype(int)).values for p in per)
    stable = list(names[hits >= 3])
    rep_rows.extend({"model": model, "target": t, "n_patients_morph": int(h)}
                    for t, h in zip(names, hits) if h >= 1)
    lines.append(f"== {model} (R2 floored at 0 per held-out sample) ==")
    for k, v in h.items():
        lo, hi = np.nanpercentile(boot[k], [2.5, 97.5])
        lines.append(f"  {k}: {v:.4g}  95% CI [{lo:.4g}, {hi:.4g}]")
        rows.append({"model": model, "metric": k, "estimate": v, "ci_low": lo, "ci_high": hi})
    for i, p in enumerate(per):
        hp = headline(p)
        rows.append({"model": model, "metric": "per_patient", "patient": pats[i], **hp})
        lines.append(f"  {pats[i]}: ratio {hp['median_comp_over_img']:.3f} (n={hp['n_img_gt_005']}), morph genes {hp['n_morph_genes']}")
    lines.append(f"  genes with morph_share>=0.10 & R2_full>0 in >=3/4 patients: {len(stable)} {stable}")
    for gname in ["ADIPOQ", "COLEC12", "NPY1R"]:
        r = full[full.target == gname].iloc[0]
        lines.append(f"  {gname}: delta_morph+ {r.delta_morph:.3f} morph_share {r.morph_share:.3f} R2_full+ {r.r2_full:.3f}")

pd.DataFrame(rows).to_csv(os.path.join(OUT, "floored_headline.csv"), index=False)
# genes with morph_share >= 0.10 and R2_full > 0 in >= 1 held-out patient; replicated = >= 3 of 4
rep = pd.DataFrame(rep_rows).sort_values(["model", "n_patients_morph"], ascending=[True, False])
rep["replicated_ge3of4"] = rep["n_patients_morph"] >= 3
rep.to_csv(os.path.join(OUT, "replicated_morph_genes.csv"), index=False)
open(os.path.join(OUT, "floored_summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

# ---- figures: fig2 (floored) via 05_figures' own fig2 function, and the floored headline-CI panel
spec = importlib.util.spec_from_file_location("f5", os.path.join(ROOT, "code", "05_figures.py"))
f5 = importlib.util.module_from_spec(spec); spec.loader.exec_module(f5)
f5.SEGMENTS["noise"] = "#EEEEEE"          # light grey; edges drawn below so it stays visible
_orig_bars = f5.stacked_bars
def stacked_bars(ax, df, positions, width=0.8, labels=True):
    _orig_bars(ax, df, positions, width, labels)
    for patch in ax.patches[-len(df):]:        # last segment drawn = noise: thin edge
        patch.set_edgecolor("#999999"); patch.set_linewidth(0.3)
f5.stacked_bars = stacked_bars
def save(fig, name):
    ax = fig.axes[0]
    if ax.get_legend() is not None:            # legend outside the plot area, below the bars
        ax.get_legend().remove()
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=11,
                   bbox_to_anchor=(0.5, -0.06))
    for ext in ["png", "pdf"]:
        fig.savefig(os.path.join(FIG, f"{name}_floored.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
f5.save = save
for model in MODELS:
    f5.fig2_decomposition(pd.read_csv(os.path.join(OUT, f"decomposition_floored_{model}.csv")), model)

hd = pd.DataFrame(rows)
fig, axes = plt.subplots(1, 2, figsize=(9, 4))
for ax, metric, label in zip(axes, ["median_comp_over_img", "n_morph_genes"],
                             ["median R²comp / R²img\n(genes with R²img > 0.05)",
                              "genes with morph share ≥ 0.10\nand R²full > 0"]):
    for i, model in enumerate(MODELS):
        e = hd[(hd.model == model) & (hd.metric == metric)].iloc[0]
        ax.errorbar(i, e.estimate, yerr=[[e.estimate - e.ci_low], [e.ci_high - e.estimate]],
                    fmt="D", color="black", capsize=6)
        pp = hd[(hd.model == model) & (hd.metric == "per_patient")][metric]
        ax.scatter(np.full(len(pp), i) + np.linspace(-0.12, 0.12, len(pp)), pp,
                   color=["#4C72B0", "#DD8452"][i], alpha=0.7, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(MODELS); ax.set_xlim(-0.5, 1.5); ax.set_ylabel(label)
    if metric == "median_comp_over_img":
        ax.axhline(1, ls="--", color="grey", lw=1)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
fig.suptitle("Headline numbers, R² floored at 0 (diamond: pooled, 95% CI over patients; dots: patients)", fontsize=10)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(FIG, f"fig_headline_ci_floored.{ext}"), dpi=300, bbox_inches="tight")
