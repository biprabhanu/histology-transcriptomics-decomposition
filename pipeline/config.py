"""
Shared settings for the DH607 core pipeline.

Question: for each gene, how much of an H&E-to-expression model's score is
  (1) unreachable measurement noise, (2) cell-type composition, (3) real morphology?
Data: HEST-1k "Xenium v1 Human Breast FFPE with Biomarkers & Housekeeping Genes
      Custom Panel" - 11 sections from 4 patients, 280-gene panel.
"""
import os
import numpy as np

# ---------------- paths ----------------
NAS = os.environ.get("DH607_STORAGE", "/path/to/storage")
PROJECT = os.path.join(NAS, "dh607")
HEST_DIR = os.path.join(PROJECT, "data", "hest")          # raw HEST files
SPOT_DIR = os.path.join(PROJECT, "data", "spots")         # output of step 02
FEATURE_DIR = os.path.join(PROJECT, "features")           # output of step 03
RESULT_DIR = os.path.join(PROJECT, "results", "core")     # output of steps 04-05

# ---------------- samples ----------------
# sample id -> patient (from HEST metadata). Splits are always by patient.
SAMPLE_TO_PATIENT = {
    "TENX191": "P1", "TENX192": "P1", "TENX193": "P1",
    "TENX195": "P2", "TENX196": "P2",
    "TENX197": "P3", "TENX198": "P3", "TENX199": "P3",
    "TENX200": "P4", "TENX201": "P4", "TENX202": "P4",
}

# ---------------- spots ----------------
SPOT_SIZE_UM = 100.0          # HEST pseudo-spots for Xenium: 100 um squares, 100 um apart
MIN_QV = 20                   # standard Xenium transcript quality cut-off
MIN_CELLS_PER_SPOT = 3        # keep spots that actually contain tissue
MIN_COUNTS_PER_SPOT = 50      # total target (half B) counts per spot
MIN_TRANSCRIPTS_PER_CELL = 10 # in half A, needed to call a cell type

# ---------------- data thinning ----------------
# Every transcript is randomly sent to one part (binomial thinning of Poisson counts):
#   A  (50%) -> defines cell types / composition
#   B1 (25%) and B2 (25%) -> together form the prediction target B = B1 + B2;
#                            their agreement gives each gene's noise ceiling
THIN_A, THIN_B1 = 0.50, 0.25
SEED = 0

# ---------------- cell typing ----------------
LEIDEN_RESOLUTION = 1.0
CELL_TYPE_MARKERS = {
    "Epithelial_tumour": ["EPCAM", "KRT19", "CEACAM6", "CLDN4", "FOXA1", "GATA3", "ERBB2", "TACSTD2", "S100A14", "CD24"],
    "Myoepithelial":     ["KRT5", "KRT14", "KRT17", "TP63", "OXTR", "KRT15"],
    "Fibroblast":        ["FBLN1", "SFRP4", "PDGFRB", "MMP2", "CCDC80", "LAMA2", "COLEC12", "CXCL12"],
    "Perivascular_SMC":  ["MYH11", "RGS5", "CNN1", "ACTA2", "TAGLN"],
    "Endothelial":       ["PECAM1", "VWF", "CLEC14A", "PLVAP", "AQP1", "FLT1"],
    "T_cell":            ["CD3E", "TRAC", "CD8A", "IL7R", "CCL5", "TIGIT", "CTLA4"],
    "B_plasma":          ["MS4A1", "BANK1", "TNFRSF13C", "TENT5C", "CD27"],
    "Myeloid":           ["CD68", "CD163", "LYZ", "ITGAX", "APOC1", "MMP12", "CHIT1"],
    "Mast":              ["CPA3", "KIT"],
    "Adipocyte":         ["ADIPOQ", "LPL", "LEPR"],
}

# ---------------- hypothesis programmes ----------------
PROGRAMMES = {
    "Proliferation": ["MKI67", "TOP2A", "CENPF", "BIRC5", "CCNB1", "UBE2C", "MYBL2", "TYMS", "RRM2", "AURKA",
                      "CDC20", "KIF2C", "NDC80", "NUF2", "CEP55", "MELK", "PTTG1", "ANLN", "EXO1", "UBE2T"],
    "Hypoxia":       ["VEGFA", "SLC2A1", "PGK1", "LDHA", "NDRG1", "ENO1", "NDUFA4L2", "PKM"],
    "EMT":           ["MMP2", "MMP14", "TAGLN", "PDGFRB", "TIMP1", "LOXL1", "SFRP4", "FBLN1", "CXCL12",
                      "ITGB1", "TGFB2", "LAMC2", "COL4A1", "LAMA2", "VCL"],
    "Interferon":    ["ISG15", "CXCL10", "CD74", "CCL5"],
}
MIN_GENES_PER_PROGRAMME = 3

# ---------------- genes ----------------
CONTROL_PREFIXES = ("NegControl", "Unassigned", "BLANK", "Deprecated", "Intergenic", "antisense")
MIN_SPOT_DETECTION = 0.10     # keep genes detected (B > 0) in >= 10% of spots

# ---------------- models ----------------
RIDGE_ALPHAS = np.logspace(-1, 9, 21)   # G3 fix: superset of logspace(-1,6,15); 1e6 was a boundary (89% of residual fits)
IMAGE_MODELS = ["resnet50", "uni2h"]   # uni2h is used only if gated access works
BATCH_SIZE = 128
