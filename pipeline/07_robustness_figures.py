"""
Step 07: robustness-round figures (fig5, fig_headline_ci; fig6 is a stub).
Reads results/robustness/*; writes 300 dpi PNG + PDF to results/robustness/figures/.
Style (colours, rcParams, save helper) follows code/05_figures.py; results/core is not touched.
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROB_DIR = os.path.join(PROJECT, "results", "robustness")
FIG_DIR = os.path.join(ROB_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

ENCODERS = ["resnet50", "phikon2"]
COLOURS = {"resnet50": "#4C72B0", "phikon2": "#DD8452"}   # same blue / orange as 05_figures segments
Y_FLOOR = -0.3          # clip extreme negative R² so the bulk is readable (as in 05_figures)
ENCODES_THRESHOLD = 0.3  # R² above this = the encoder "encodes" the cell type
# interpretable cell types only; the other five have near-constant targets (R² undefined/exploding)
SHOWN_TYPES = ["Epithelial_tumour", "T_cell", "Fibroblast", "Myeloid", "Endothelial", "Other"]

plt.rcParams.update({"font.size": 12, "axes.titlesize": 14, "axes.labelsize": 13,
                     "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {name}", flush=True)


def fig5_image_to_composition():
    """Image-only R² for predicting each cell-type fraction: bars = mean over held-out patients,
    dots = individual patients."""
    df = pd.read_csv(os.path.join(ROB_DIR, "image_to_composition.csv"))
    pcols = ["r2_P1", "r2_P2", "r2_P3", "r2_P4"]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    width = 0.38
    rng = np.random.default_rng(0)
    for k, enc in enumerate(ENCODERS):
        sub = df[df["encoder"] == enc].set_index("cell_type").loc[SHOWN_TYPES]
        x = np.arange(len(SHOWN_TYPES)) + (k - 0.5) * width
        ax.bar(x, sub["r2_mean"].clip(lower=Y_FLOOR), width * 0.95, color=COLOURS[enc],
               alpha=0.55, label=enc, linewidth=0)
        for xi, (_, row) in zip(x, sub.iterrows()):
            vals = row[pcols].astype(float).clip(lower=Y_FLOOR).values
            ax.scatter(xi + rng.uniform(-0.06, 0.06, len(vals)), vals, s=22,
                       color=COLOURS[enc], edgecolor="black", linewidth=0.5, zorder=3)
    ax.axhline(ENCODES_THRESHOLD, color="k", ls="--", lw=1)
    ax.text(len(SHOWN_TYPES) - 0.45, ENCODES_THRESHOLD + 0.015, "encodes cell type",
            ha="right", va="bottom", fontsize=10)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(np.arange(len(SHOWN_TYPES)))
    ax.set_xticklabels([t.replace("_", " ") for t in SHOWN_TYPES])
    ax.set_ylim(Y_FLOOR - 0.02, 1)
    ax.set_ylabel("held-out-patient R², image only")
    ax.set_title(f"R² below {Y_FLOOR} drawn at {Y_FLOOR}", fontsize=9, loc="right")
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.92))
    fig.text(0.01, -0.02, "Bars: mean over 4 held-out patients; dots: individual patients. "
             "Not shown: Mast, B_plasma, Perivascular_SMC, Myoepithelial, Adipocyte\n"
             "(near-constant targets; R² undefined or exploding).", fontsize=9, ha="left", va="top")
    save(fig, "fig5_image_to_composition")


def parse_summary(path):
    """{encoder: {metric: (point, lo, hi)}} from the patient-bootstrap blocks of the summary."""
    out, enc = {}, None
    for line in open(path):
        m = re.match(r"== (\w+): patient bootstrap", line)
        if m:
            enc = m.group(1)
            out[enc] = {}
            continue
        m = re.match(r"\s+(\w+)\s+point=(\S+)\s+95% CI \[(\S+), (\S+)\]", line)
        if m and enc:
            out[enc][m.group(1)] = tuple(float(m.group(i)) for i in (2, 3, 4))
    return out


def fig_headline_ci():
    """Two panels: point estimate with patient-bootstrap 95% CI, plus per-patient dots."""
    pp = pd.read_csv(os.path.join(ROB_DIR, "per_patient.csv"))
    pp = pp[pp["patient"] != "pooled"]
    ci = parse_summary(os.path.join(ROB_DIR, "per_patient_summary.txt"))
    panels = [("median_comp_over_img", "median R²(composition) / R²(image)\n(genes with R²(image) > 0.05)", 1.0),
              ("n_morph_share_ge_010_r2full_pos", "genes with morph share ≥ 0.10\nand R²(full) > 0", None)]
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.2))
    rng = np.random.default_rng(0)
    for ax, (metric, label, ref) in zip(axes, panels):
        for i, enc in enumerate(ENCODERS):
            point, lo, hi = ci[enc][metric]
            ax.errorbar(i, point, yerr=[[point - lo], [hi - point]], fmt="D", ms=9,
                        color=COLOURS[enc], ecolor="black", elinewidth=1.8, capsize=6, zorder=3)
            vals = pp[pp["model"] == enc][metric].values
            ax.scatter(i + 0.22 + rng.uniform(-0.04, 0.04, len(vals)), vals, s=22,
                       color=COLOURS[enc], alpha=0.6, edgecolor="black", linewidth=0.5, zorder=2)
        if ref is not None:
            ax.axhline(ref, color="k", ls="--", lw=1)
            ax.text(1.55, ref, "1.0", va="bottom", ha="right", fontsize=9)
        ax.set_xticks(range(len(ENCODERS)))
        ax.set_xticklabels(ENCODERS)
        ax.set_xlim(-0.5, 1.6)
        ax.set_ylabel(label, fontsize=11)
        if ref is None:
            ax.set_ylim(bottom=0)
    fig.suptitle("Diamond + bar: pooled estimate, 95% CI over patients; dots: single patients",
                 fontsize=9, y=0.02)
    fig.tight_layout()
    save(fig, "fig_headline_ci")


def other_percent():
    """% of spots typed 'Other' per sample: before (per-sample labels, 'old' rows) and after
    (joint labels, 'Other' column) from data/spots_joint/old_vs_new.csv."""
    ct = pd.read_csv(os.path.join(PROJECT, "data", "spots_joint", "old_vs_new.csv"))
    rows = []
    for sample, g in ct.groupby("sample"):
        total = g.drop(columns=["sample", "old"]).values.sum()
        rows.append((sample, 100 * g.loc[g["old"] == "Other"].drop(columns=["sample", "old"]).values.sum() / total,
                     100 * g["Other"].sum() / total))
    return pd.DataFrame(rows, columns=["sample", "before", "after"])


def fig6_joint_vs_persample():
    """Per-sample cell typing (before) vs joint cell typing (after): composition R², morphology
    gain, and share of spots labelled Other."""
    fig = plt.figure(figsize=(14, 5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.8], wspace=0.35)
    axA, axB, axC = (fig.add_subplot(g) for g in gs)
    lo = Y_FLOOR - 0.02
    for enc in ENCODERS:
        before = pd.read_csv(os.path.join(PROJECT, "results", "core", f"decomposition_{enc}.csv")).set_index("target")
        after = pd.read_csv(os.path.join(ROB_DIR, "joint", f"decomposition_{enc}.csv")).set_index("target")
        after = after.loc[before.index.intersection(after.index)]
        before = before.loc[after.index]
        for ax, col in ((axA, "r2_comp"), (axB, "delta_morph")):
            ax.scatter(before[col].clip(lower=Y_FLOOR), after[col].clip(lower=Y_FLOOR), s=14,
                       alpha=0.55, color=COLOURS[enc], label=enc, linewidth=0)
    for ax, xl, yl, tt in ((axA, "R² composition, per-sample labels", "R² composition, joint labels", "A  composition"),
                           (axB, "Δ morphology, per-sample labels", "Δ morphology, joint labels", "B  morphology")):
        ax.plot([lo, 1], [lo, 1], "k--", lw=1)
        ax.set_xlim(lo, 0.6 if ax is axB else 1)
        ax.set_ylim(lo, 0.6 if ax is axB else 1)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(tt, loc="left", fontsize=12)
    axA.legend(frameon=False, loc="upper left")
    axA.text(1, lo + 0.02, f"axes clipped at {Y_FLOOR}", ha="right", fontsize=9)
    axB.text(0.6, lo + 0.02, f"axes clipped at {Y_FLOOR}", ha="right", fontsize=9)
    oth = other_percent()
    x = np.arange(len(oth))
    axC.bar(x - 0.2, oth["before"], 0.4, color="#CCCCCC", label="per-sample")
    axC.bar(x + 0.2, oth["after"], 0.4, color="#555555", label="joint")
    axC.set_xticks(x)
    axC.set_xticklabels(oth["sample"], rotation=30, ha="right")
    axC.set_ylabel("% spots labelled Other")
    axC.set_title("C  'Other' share", loc="left", fontsize=12)
    axC.legend(frameon=False)
    fig.suptitle("Joint cell typing: composition not stronger, morphology unchanged", fontsize=14, y=1.0)
    save(fig, "fig6_joint_vs_persample")
    print(oth.round(1).to_string(), flush=True)


if __name__ == "__main__":
    fig5_image_to_composition()
    fig_headline_ci()
    fig6_joint_vs_persample()
