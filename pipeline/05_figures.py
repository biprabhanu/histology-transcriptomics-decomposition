"""
Step 05: poster figures from results/core/decomposition_<model>.csv.
Writes 300 dpi PNG + PDF to results/core/figures/.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config

FIG_DIR = os.path.join(config.RESULT_DIR, "figures")
MODELS = ["resnet50", "phikon2"]
# segment name -> colour, in stacking order (bottom to top)
SEGMENTS = {"composition": "#4C72B0", "morphology": "#DD8452",
            "reachable, unexplained": "#CCCCCC", "noise": "#555555"}

plt.rcParams.update({"font.size": 12, "axes.titlesize": 14, "axes.labelsize": 13,
                     "axes.spines.top": False, "axes.spines.right": False})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {name}", flush=True)


Y_FLOOR = -0.3   # scatter axes: clip extreme negative cross-patient R² so the bulk is readable


def segments(df):
    """The four stacked-bar segments. Composition and morphology are clipped at >= 0 and
    to what fits under the ceiling, so every bar sums to exactly 1 (noise = 1 - ceiling)."""
    ceiling = df["ceiling"].clip(0, 1)
    comp = df["r2_comp"].clip(lower=0).clip(upper=ceiling)
    morph = df["delta_morph"].clip(lower=0).clip(upper=ceiling - comp)
    return pd.DataFrame({
        "composition": comp,
        "morphology": morph,
        "reachable, unexplained": ceiling - comp - morph,
        "noise": 1 - ceiling,
    })


def stacked_bars(ax, df, positions, width=0.8, labels=True):
    seg = segments(df)
    bottom = np.zeros(len(df))
    for name, colour in SEGMENTS.items():
        ax.bar(positions, seg[name].values, width, bottom=bottom, color=colour,
               label=name if labels else None, linewidth=0)
        bottom += seg[name].values


def short(name):
    return name.replace("PROG_", "")


def fig1_noise(df, model):
    """Ceiling vs R2_img: apparent failure is often noise."""
    fig, ax = plt.subplots(figsize=(6, 6))
    genes, progs = df[~df["is_programme"]], df[df["is_programme"]]
    ax.scatter(genes["ceiling"], genes["r2_img"].clip(lower=Y_FLOOR), s=14, alpha=0.6, color="#4C72B0", label="gene")
    ax.scatter(progs["ceiling"], progs["r2_img"].clip(lower=Y_FLOOR), s=60, color="#C44E52", label="programme")
    for _, row in progs.iterrows():
        ax.annotate(short(row["target"]), (row["ceiling"], row["r2_img"]),
                    xytext=(5, 5), textcoords="offset points", color="#C44E52")
    ax.plot([0, 1], [0, 1], "k--", lw=1, label="y = x (ceiling)")
    ax.set_xlim(0, 1)
    ax.set_ylim(Y_FLOOR - 0.02, 1)   # R2 below Y_FLOOR drawn at the floor (outliers)
    ax.set_xlabel("noise ceiling (max achievable R²)")
    ax.set_ylabel(f"R² image only ({model})")
    ax.set_title(f"R² below {Y_FLOOR} shown at {Y_FLOOR}", fontsize=9, loc="right")
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.88))
    save(fig, f"fig1_noise_{model}")


def fig2_decomposition(df, model):
    """Stacked bars for all genes sorted by R2_full, plus top 20 by morph_share."""
    genes = df[~df["is_programme"]].sort_values("r2_full", ascending=False)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={"width_ratios": [2, 1]})
    stacked_bars(ax1, genes, np.arange(len(genes)), width=1.0)
    ax1.set_xlim(-0.5, len(genes) - 0.5)
    ax1.set_ylim(0, 1)
    ax1.set_xticks([])
    ax1.set_xlabel(f"{len(genes)} genes, sorted by R² full")
    ax1.set_ylabel("fraction of variance")
    ax1.set_title(f"Per-gene decomposition ({model})")
    ax1.legend(frameon=False, loc="upper right", fontsize=10)

    top = genes[genes["r2_full"] > 0].nlargest(20, "morph_share")   # skip negative-R² artefacts
    stacked_bars(ax2, top, np.arange(len(top)), labels=False)
    ax2.set_xticks(np.arange(len(top)))
    ax2.set_xticklabels(top["target"], rotation=90, fontsize=10)
    ax2.set_ylim(0, 1)
    ax2.set_title("Top 20 genes by morphology share (R² full > 0)")
    save(fig, f"fig2_decomposition_{model}")


def fig3_programmes(tables):
    """Stacked bars for the 4 programmes, one bar per model side by side."""
    progs = [t for t in next(iter(tables.values())).query("is_programme")["target"]]
    width = 0.8 / len(tables)
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, (model, df) in enumerate(tables.items()):
        rows = df.set_index("target").loc[progs].reset_index()
        pos = np.arange(len(progs)) + (i - (len(tables) - 1) / 2) * width
        stacked_bars(ax, rows, pos, width=width * 0.95, labels=(i == 0))
        for p in pos:
            ax.text(p, 1.01, model, ha="center", va="bottom", fontsize=9, rotation=90)
    ax.set_xticks(np.arange(len(progs)))
    ax.set_xticklabels([short(p) for p in progs])
    ax.set_ylim(0, 1.25)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_ylabel("fraction of variance")
    ax.set_title("Hypothesis programmes")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=4, fontsize=10)
    save(fig, "fig3_programmes")


def fig4_encoders(tables):
    """Does a pathology model see more morphology, or just more composition?"""
    a, b = tables["resnet50"], tables["phikon2"]
    merged = a.merge(b, on="target", suffixes=("_resnet50", "_phikon2"))
    merged = merged[~merged["is_programme_resnet50"]]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    x = merged["delta_morph_resnet50"].clip(lower=Y_FLOOR)
    y = merged["delta_morph_phikon2"].clip(lower=Y_FLOOR)
    lo, hi = min(x.min(), y.min()), max(x.max(), y.max())
    ax1.scatter(x, y, s=14, alpha=0.6, color="#4C72B0")
    ax1.plot([lo, hi], [lo, hi], "k--", lw=1)
    ax1.set_xlabel("ΔR² morphology, resnet50")
    ax1.set_ylabel("ΔR² morphology, phikon2")
    ax1.set_title("Morphology beyond composition")

    for model, colour in zip(MODELS, ["#4C72B0", "#DD8452"]):
        ax2.scatter(merged[f"r2_comp_{model}"].clip(lower=Y_FLOOR), merged[f"r2_img_{model}"].clip(lower=Y_FLOOR), s=14, alpha=0.6,
                    color=colour, label=model)
    hi = max(merged[[f"r2_comp_{m}" for m in MODELS] + [f"r2_img_{m}" for m in MODELS]].max())
    ax2.plot([0, hi], [0, hi], "k--", lw=1)
    ax2.set_xlabel("R² composition only")
    ax2.set_ylabel("R² image only")
    ax2.set_title(f"Image score vs composition score (values < {Y_FLOOR} at {Y_FLOOR})")
    ax2.legend(frameon=False)
    save(fig, "fig4_encoders")


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    tables = {}
    for model in MODELS:
        path = os.path.join(config.RESULT_DIR, f"decomposition_{model}.csv")
        if os.path.exists(path):
            tables[model] = pd.read_csv(path)
    if not tables:
        raise SystemExit("no decomposition CSVs found")

    for model, df in tables.items():
        fig1_noise(df, model)
        fig2_decomposition(df, model)
    fig3_programmes(tables)
    if len(tables) == 2:
        fig4_encoders(tables)


if __name__ == "__main__":
    main()
