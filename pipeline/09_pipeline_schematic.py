"""Figure 0: pipeline schematic (boxes + arrows). Writes results/figures/fig0_pipeline.{png,pdf}."""
import os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = os.environ.get("DH607_STORAGE", "/path/to/storage") + "/dh607"
OUT = os.path.join(ROOT, "results", "figures")
os.makedirs(OUT, exist_ok=True)

fig, ax = plt.subplots(figsize=(12.5, 5.4))
ax.set_xlim(0, 12.1); ax.set_ylim(0.2, 5.2); ax.axis("off")

def box(x, y, w, h, text, colour):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05,rounding_size=0.12",
                                fc=colour, ec="#333333", lw=1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9, wrap=True)

def arrow(x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.2))

BLUE, ORANGE, GREY, GREEN = "#DCE6F2", "#F8E0CC", "#EEEEEE", "#E2EFDA"
# inputs (left), processing (middle), models (right); no arrow crosses a box
box(0.1, 3.75, 2.6, 1.1, "Xenium transcripts\n280 genes, qv ≥ 20", GREY)
box(0.1, 0.4, 2.6, 1.1, "H&E patch\n224 px, 100 µm spot", GREY)
box(3.2, 3.75, 2.5, 1.1, "Binomial thinning\nA 50% · B1 25% · B2 25%", BLUE)
box(3.2, 0.4, 2.5, 1.1, "Image encoders\nResNet-50 · Phikon-v2", ORANGE)
box(6.2, 3.9, 2.7, 1.0, "Part A: Leiden cell types\n→ composition per spot", BLUE)
box(6.2, 2.6, 2.7, 1.0, "B1 vs B2:\nnoise ceiling (Spearman–Brown)", BLUE)
box(6.2, 1.3, 2.7, 1.0, "Target B = B1 + B2\nlog-normalised, centred", BLUE)
box(9.4, 1.3, 2.5, 3.6, "Ridge regression,\nleave-one-patient-out\n\ncomposition → B\nimage → B\ncomposition + image\non residuals\n\nΔmorph = R²full − R²comp\nshares of noise ceiling", GREEN)

arrow(2.7, 4.3, 3.2, 4.3)
arrow(2.7, 0.95, 3.2, 0.95)
arrow(5.7, 4.5, 6.2, 4.4); arrow(5.7, 4.1, 6.2, 3.1); arrow(5.7, 3.85, 6.2, 1.8)
arrow(8.9, 4.4, 9.4, 4.4); arrow(8.9, 3.1, 9.4, 3.1); arrow(8.9, 1.8, 9.4, 1.8)
ax.plot([5.7, 10.65], [0.95, 0.95], color="#333333", lw=1.2)     # encoder features route below
arrow(10.65, 0.95, 10.65, 1.3)
ax.set_title("Separating noise, composition and morphology per gene", fontsize=12, loc="left")
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(OUT, f"fig0_pipeline.{ext}"), dpi=300, bbox_inches="tight")
print("saved")
