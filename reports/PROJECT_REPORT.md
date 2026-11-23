# DH607 — Project report: noise vs composition vs morphology in H&E → gene-expression prediction

## Changelog
- **2026-11-24 — Final fixes:** fig2_floored legend moved below the plot and noise segment made visible (light grey, thin edge). Added `floored/replicated_morph_genes.csv`. Reworded decision points 1 and 2 and the morphology-gene statement; no new analysis.
- **2026-11-24 — Robustness round (pre-registered add-on analyses; main method unchanged):**
  - Image → composition (06).
  - Benchmark sanity check.
  - Per-patient consistency with bootstrap CIs (06b; 04 gained an opt-in `--save-predictions` flag, and default output is byte-identical).
  - Joint cell typing with Harmony (02b, 04b).
  - Robustness figures (07).
- **2026-11-24 — Method change (pre-specified by the user):**
  - All share calculations floor held-out R² at 0 per sample before averaging. This fixes unstable held-out R², e.g. NPY1R Δmorph 2.24 in P4.
  - Script `08_floored_shares.py` writes `results/robustness/floored/` and regenerates fig2 and the headline-CI panel (`*_floored`).
  - `results/core` itself is untouched.
- **2026-11-24 — G4 passed:**
  - Independent reviewer: PASS on patient separation, alpha selection, residual fitting (note: residuals are in-sample, training only), ceiling and shuffled control. All RESULTS_SUMMARY numbers match the CSVs.
  - Applied 3 wording corrections: Phikon-v2 Δmorph is −0.0000, the Interferon verdict note, and the NPY1R artefact.
  - Figures fixed (presentation only; no change to the analysis):
    - stacked bars now sum to exactly 1 (composition and morphology capped under the ceiling);
    - the fig2 inset shows only genes with R²_full > 0;
    - scatter axes floor negative R² at −0.3;
    - fig1 legend moved.
  - Results copied to the Mac.
- **2026-11-24 — G3 passed with caveat:** Final 4-patient run of 04 for both encoders.
  - Changed the alpha grid from `logspace(-1,6,15)` to `logspace(-1,9,21)` (a superset): 1e6 was the upper boundary, chosen for 89% of residual fits.
  - Shares are now computed only where ceiling ≥ 0.2.
  - G3(c) reached 91.4% / 80.0%, not 95%. The user accepted this with a caveat.
- **2026-11-24 — data:** The three missing files arrived on the server through a retry loop, all SHA-256 verified. Watcher `fin` ran 02 (TENX202, TENX195) and 03 (TENX195), all exit 0.
- **2026-11-23 — scope change:** Cut to one section per patient (fast path) because the campus network kept resetting Hugging Face downloads.
- **2026-11-23 — G1 tweak:** Lowered the "Other" cell-type threshold from 0.5 to 0.3, because TENX191 had 32% Other.
- **2026-11-23 — environment fixes:**
  - Pinned scipy to 1.13.1 and numba to 0.61.2.
  - Step 02 uses the sklearn exact kNN (`transformer="sklearn"`) because of a umap-learn 0.6 break.
  - Installed the ResNet-50 V2 weights from the HF mirror `timm/resnet50.tv2_in1k`, because pytorch.org was blocked.

## 1. Question & hypothesis
The per-gene score of an H&E → expression model mixes three things:
1. measurement noise,
2. cell-type composition,
3. genuine morphology signal.

We separate them per gene, each as a fraction of a per-gene noise ceiling.

The hypothesis had three parts:
- much apparent failure is noise;
- much apparent success is composition;
- a few programmes (proliferation, hypoxia, EMT, interferon) are truly visible in morphology.

## 2. Dataset
- **Source:** HEST-1k v1.3.0 (`MahmoodLab/hest` on Hugging Face). Study: *"Xenium v1 Human Breast FFPE with Biomarkers & Housekeeping Genes Custom Panel"* (10x Genomics Xenium public datasets, https://www.10xgenomics.com/datasets).
- **Patients and sections:** 4 patients, 11 sections:

  | Patient | Sections |
  |---|---|
  | P1 | TENX191, 192, 193 |
  | P2 | TENX195, 196 |
  | P3 | TENX197, 198, 199 |
  | P4 | TENX200, 201, 202 |

- **Sections used:** one per patient: TENX191 (P1), TENX195 (P2), TENX199 (P3), TENX202 (P4). The campus network made full downloads infeasible, and one section per patient keeps leave-one-patient-out intact.
- **Panel and image:** 280-gene panel; all 280 genes pass detection (≥10% of spots). H&E pixel ≈ 0.274 µm.
- **Spots:** HEST pseudo-spots, 100 µm squares on a 100 µm grid (365 px). Only spots with an H&E patch are used: 224 px at 0.5 µm/px, centred on the spot.

| Sample | Patient | Spots (patch) | Kept (≥3 cells, ≥50 B counts) | Typed cells | Transcripts in spots |
|---|---|---|---|---|---|
| TENX191 | P1 | 4,949 | 4,605 | 193,448 | 99.8% |
| TENX195 | P2 | 6,254 | 5,971 | 319,213 | 99.9% |
| TENX199 | P3 | 2,259 | 2,247 | 146,234 | 99.5% |
| TENX202 | P4 | 4,291 | 4,253 | 197,820 | 99.8% |

## 3. Processing pipeline
1. **01_download_core.py**: downloads metadata, st (h5ad), patches (h5) and transcripts (parquet) per sample. Replacements:
   - `dl_parallel.py`, a wrapper with 8 workers;
   - `.refetch.py`, which retries per file with `hf_hub_download`.
2. **02_build_spots.py SAMPLE**: builds spots and cell types.
   - **QC:** transcripts with qv ≥ 20 and `is_gene` only; control probes (NegControl, BLANK, Unassigned, Deprecated, Intergenic, antisense) removed.
   - **Binomial thinning** of each transcript into A (50%), B1 (25%) and B2 (25%), seed 0.
   - **Transcript → spot:** H&E coordinates are assigned to the 100 µm square around each spot centre.
   - **Cell typing on part A:** cells with ≥10 A transcripts are clustered with Leiden (resolution 1.0, kNN = 15, sklearn exact kNN). Each cluster gets the best-scoring type from a 10-type marker table, or "Other" if its best score is < 0.3.
   - **Composition per spot:** type fractions plus n_cells.
3. **03_extract_features.py**: computes image features on GPU 0.
   - ResNet-50 (ImageNet V2), global average pool → 2048 features.
   - Phikon-v2 (owkin), CLS token → 1024 features.
   - UNI2-h, H-optimus, Virchow2 and CONCH are gated, with no access.
4. **04_decompose.py MODEL**: modelling (§4).
5. **05_figures.py**: figures 1–4, 300 dpi PNG + PDF.

**Cell-type summary** (cells, % of typed cells):

| Sample | Epithelial | Other | Notable |
|---|---|---|---|
| TENX191 | 38,649 (20%) | 42,832 (22%) | two weak "Other" clusters (26k, 13k cells) |
| TENX195 | 79,105 (25%) | 30,082 (9%) | T cells 65k, B/plasma 25k |
| TENX199 | 78,356 (54%) | 0 | only 38 myoepithelial cells, no B/plasma |
| TENX202 | 99,745 (50%) | 18,990 (10%) | "Mast" cluster with 13.7k cells (score 1.18), which is implausible |

Cluster tables: `data/spots/<SAMPLE>_clusters.csv`.

## 4. Modelling
- **Targets:**
  - for C ∈ {B, B1, B2}: `log1p(C / C.sum(1) · 100)`;
  - programme scores: genes z-scored within sample, then averaged over each programme's genes (needs ≥ 3 genes);
  - B is centred per sample.
- **Features:**
  - composition = fractions + log1p(n_cells);
  - image = encoder features;
  - each column z-scored within sample. Constant columns become 0.
- **Noise ceiling:** per sample and target, r = Pearson(B1, B2). Ceiling = 2r/(1+r) (Spearman–Brown), clipped to [0, 1], then averaged over samples.
- **Models:** SVD ridge with no intercept.
  - `comp`: composition only.
  - `img`: image only.
  - `full`: comp fit plus an image ridge fitted on the training comp residuals.
  - `img_shuffled`: image rows permuted within each training sample (ResNet-50 only).
- **Validation:**
  - outer leave-one-patient-out (4 folds);
  - alpha chosen per target by inner leave-one-patient-out on the 3 training patients;
  - alpha grid `logspace(-1, 9, 21)`;
  - metrics computed per held-out sample, then averaged.
- **Metrics:**
  - R²_comp, R²_img, R²_full;
  - naive r;
  - Δmorph = R²_full − R²_comp;
  - partial R²;
  - comp_share = R²_comp / ceiling and morph_share = Δmorph / ceiling, only where ceiling ≥ 0.2.
- **Categories:**
  - `noise_limited` if ceiling < 0.2;
  - else `morphology_visible` if morph_share ≥ 0.10;
  - else `composition_explained` if R²_img ≥ 0.10 · ceiling;
  - else `not_predicted`.

## 5. Quality gates
- **G1 (cell types): pass.**
  - Transcripts inside spots: 99.5–99.9%.
  - "Other": 22% / 9% / 0% / 10%, all ≤ 25% after the 0.5 → 0.3 tweak.
  - Epithelial is the largest named type in every sample. In TENX191, "Other" is the single largest label.
- **G2 (features): pass.** 100.0% of feature barcodes match spot barcodes, for all 4 samples and both encoders.
- **G3 (modelling): pass with caveat.**

  | Check | Target | ResNet-50 | Phikon-v2 |
  |---|---|---|---|
  | (a) median ceiling by mean-count quartile | increasing | 0.588 / 0.702 / 0.706 / 0.819 | same (same targets) |
  | (b) shuffled-image median abs(R²) | < 0.01 | 0.0007 | not run (control is ResNet-50 only) |
  | (c) R²_full ≥ R²_comp − 0.01 | > 95% of genes | 91.4% (64% before the grid fix) | 80.0% (51% before) |
  | (d) 4 patients, 4 held-out samples, each scored once | all | yes | yes |

  - **Why (c) misses:** the residual image model now picks the largest alpha for 75% of genes (ResNet-50) and 61% (Phikon-v2). The remaining failures are genes where the inner CV over 3 patients found image signal that did not transfer to the held-out patient. ResNet-50: 24 genes fail, median gap 0.020.
- **G4: pass.**
  - All 4 patients present.
  - Both CSVs, both summaries and the figures exist.
  - The reviewer agrees with every headline number. Its flags are folded in: in-sample residuals; negative-R² outliers such as NPY1R; figure fixes.
  - One remaining note: in the fig2 main panel, NPY1R (far right) still shows a morphology bar that is an artefact of negative R².

## 6. Results
Numbers come from `results/core/decomposition_*.csv` and `summary_*.txt`.

- **Noise:**
  - median ceiling 0.707;
  - 3 of 280 genes are noise-limited (1.1%);
  - 48 genes have ceiling < 0.5. (`fig1_noise_*`)
- **Categories:**

  | Encoder | not predicted | composition explained | morphology visible | noise limited |
  |---|---|---|---|---|
  | ResNet-50 | 226 | 51 | 0 | 3 |
  | Phikon-v2 | 181 | 90 | 6 | 3 |

  Phikon-v2's six morphology-visible genes are NPY1R, ADIPOQ, ADH1B, LPL, LEPR and COLEC12.
  - NPY1R is an artefact: R²_comp −10.37, R²_full −9.87.
  - Only ADIPOQ (R²_full 0.023) and COLEC12 (0.100) have R²_full > 0.
- **Overall level:** median R²_comp −0.0025, R²_img −0.001, R²_full −0.005. 104 genes have R²_comp > 0.05.
  - Many genes have negative cross-patient R², for example ESR1 R²_comp −0.608. This reflects patient-specific expression–composition relations (e.g. receptor status) and tiny within-sample variance in some held-out samples.
- **Composition share of apparent success:** among genes with R²_img > 0.05:
  - median R²_comp / R²_img = 1.778 (ResNet-50, n = 59) and 1.061 (Phikon-v2, n = 100);
  - median comp_share = 0.344 / 0.246. (`fig2_decomposition_*`)
- **Morphology:**
  - median Δmorph −0.0001 (ResNet-50) and 0.0000 (Phikon-v2);
  - Δmorph > 0.02 in 8 / 19 genes;
  - top morph_share genes are adipose: ADIPOQ, ADH1B, LPL, LEPR (ResNet-50 morph_share ≤ 0.096).
- **Programmes** (`fig3_programmes`):

  | Programme | Ceiling | R²_comp | R²_img (R / P) | Δmorph (R / P) | Verdict |
  |---|---|---|---|---|---|
  | Proliferation | 0.957 | 0.265 | 0.132 / 0.231 | −0.008 / 0.006 | composition |
  | EMT | 0.971 | 0.729 | 0.468 / 0.610 | 0.013 / 0.027 | composition |
  | Interferon (4 genes) | 0.926 | 0.320 | 0.075 / 0.166 | −0.000 / −0.004 | composition |
  | Hypoxia | 0.833 | −0.072 | −0.037 / −0.086 | −0.000 / −0.027 | not predicted (not noise) |

- **Encoders** (`fig4_encoders`):
  - Phikon-v2 beats ResNet-50 on R²_img in 52% of genes;
  - 100 vs 59 genes have R²_img > 0.05;
  - Δmorph correlation between encoders is r = 0.19.
  - Interpretation: more composition, little extra morphology.
- **Hypothesis check:**
  1. "Failure is noise": **not supported**. Ceilings are high on this deep panel.
  2. "Success is composition": **supported**.
  3. "Programmes visible in morphology": **not supported**. All programme morph_share < 0.03.

## 6b. Robustness analyses
All of these were decided before the numbers were seen. Outputs are in `results/robustness/`.
- **Image → composition** (`06_image_to_composition.py`, `image_to_composition.csv`, `fig5`).
  - Ridge from image features to the per-spot cell-type fractions, same LOPO and inner alpha as 04.
  - Mean R², ResNet-50 / Phikon-v2:

    | Cell type | R² ResNet-50 | R² Phikon-v2 |
    |---|---|---|
    | Epithelial | 0.42 | 0.53 |
    | T cell | 0.37 | 0.53 |
    | Fibroblast | 0.25 | 0.40 |
    | Endothelial | 0.07 | 0.38 |
    | Myeloid | 0.14 | 0.26 |
    | Other | 0.02 | 0.06 |

  - Adipocyte R² is negative (−0.81 / −1.16).
  - Mast, B/plasma, perivascular SMC and myoepithelial fractions are near-constant in some held-out samples, so their R² is undefined or explodes and is not interpreted.
  - Pre-set reading (R² ≥ 0.3): the image encodes epithelial and T cells with both encoders, plus fibroblast and endothelial with Phikon-v2.
- **Benchmark sanity** (`benchmark_check.md`).
  - Naive r, median [IQR]:

    | Encoder | All 280 genes | Top-50 genes |
    |---|---|---|
    | ResNet-50 | 0.197 [0.047–0.377] | 0.279 [0.087–0.447] |
    | Phikon-v2 | 0.297 [0.113–0.488] | 0.353 [0.113–0.568] |

  - Published HEST-Benchmark IDC: ResNet-50 ≈ 0.473, Phikon (v1) ≈ 0.533 (github.com/mahmoodlab/HEST README, read from snippets).
  - The image models score below published random-split benchmarks (top-50 r 0.28 / 0.35 vs ~0.47 / 0.53), as expected under leave-one-patient-out; encoder ranking agrees.
- **Per-patient consistency** (`06b_per_patient.py`, `per_patient.csv`, `per_patient_summary.txt`).
  - The rerun of 04 with `--save-predictions` reproduced `results/core` exactly (max abs diff 0).
  - Per held-out patient, median R²_comp/R²_img:
    - ResNet-50: P1 1.79, P2 1.52, P3 1.57, P4 2.12;
    - Phikon-v2: 1.07, 0.98, 1.17, 1.17.
  - The headline (ratio ≥ 0.9 and morphology genes ≤ 5%) holds in **4 of 4 patients** for both encoders.
  - ADIPOQ, COLEC12 and NPY1R Δmorph change sign between patients.
- **Floored shares** (`08_floored_shares.py`, `floored/floored_summary.txt`). Bootstrap over patients, 1000 resamples, seed 0:

  | Encoder | comp/img ratio [95% CI] | n | Median Δmorph | Morphology genes [95% CI] | Which | In ≥3/4 patients |
  |---|---|---|---|---|---|---|
  | ResNet-50 | 1.771 [1.557–2.003] | 105 | 9.2e-6 | 2 [0–5] | ADIPOQ, LPL | 0 |
  | Phikon-v2 | 1.118 [0.998–1.179] | 142 | 4.0e-5 | 5 [2–9] | ADIPOQ, LPL, MYH11, ADH1B, APOC1 | 1 (MYH11) |

  Interpretation: these are adipose and smooth-muscle genes, most likely **missed composition (adipocyte/stromal segmentation limits)**.
- **Joint cell typing** (`02b_joint_celltypes.py`, `04b_decompose_joint.py`, `data/spots_joint/`, `results/robustness/joint/`, `fig6`).
  - Method:
    - same part-A counts (barcodes and B counts verified identical);
    - pooled cells, ≤100k per sample, Harmony (harmonypy 2.1.0) by sample, Leiden 1.0, the same marker rule;
    - kNN label transfer to all cells.
  - "Other", TENX191/195/199/202: before 22.1 / 9.4 / 0.0 / 9.6% → after 5.8 / 11.4 / 9.5 / 6.8%.
  - TENX202 "Mast" cells became Other 51%, Epithelial 37%, Myoepithelial 5%. The per-sample Mast call was wrong.
  - TENX191 "Other" became Myoepithelial 56%, Epithelial 44%. The cluster is 95% TENX191 with a marginal score of 0.30, which indicates a residual batch effect.
  - R²_comp: median −0.0025 → −0.0030; genes with R²_comp > 0.05: 104 → 91.
  - Median Δmorph about 0 in every run. Morphology genes 0 → 0 (ResNet-50) and 2 → 1 (Phikon-v2).
- **Decision point 1: headline is ROBUST.**
  - It holds in 4/4 patients for both encoders, and the image encodes the major cell types.
  - The benchmark condition concerns absolute image performance, not the claim: the image models score below published random-split benchmarks (top-50 r 0.28 / 0.35 vs ~0.47 / 0.53), as expected under leave-one-patient-out; encoder ranking agrees.
- **Morphology genes:** only 2 (ResNet-50) / 5 (Phikon-v2) genes show image signal beyond composition, all adipose/lipid or smooth-muscle markers (ADIPOQ, LPL, ADH1B, APOC1; MYH11) — cell types our composition features capture poorly (also unstable in fig5). Only MYH11 replicates in ≥3/4 patients (`floored/replicated_morph_genes.csv`). No cell-state or programme gene shows reproducible morphology signal.
- **Decision point 2 (pre-set rule): joint cell typing did NOT strengthen composition.**
  - 91 vs 104 genes with R²_comp > 0.05; residual batch effect (cluster 0 = 95% TENX191).
  - Per the rule, per-sample labels stay the main result. Joint clustering was not iterated further.
  - Δmorph is unchanged: morphology genes 0 → 0 and 2 → 1, within the unfloored per-patient bootstrap CIs [0, 4] and [1, 9].
  - Under joint labels the comp/img ratio falls: ResNet-50 1.778 → 1.408, Phikon-v2 1.061 → 0.961.
  - Robustness statement: **the morphology conclusion is unchanged under an alternative joint labelling.**

## 7. Decisions log
| Decision | Why |
|---|---|
| Xenium breast study in HEST-1k (fixed in the master prompt) | Single-cell transcripts allow thinning and cell typing on the same tissue as the target |
| ResNet-50 + Phikon-v2; UNI2-h, Virchow2, H-optimus and CONCH skipped | Gated, no access; Phikon-v2 has open weights |
| Fast path: one section per patient | Campus network reset Hugging Face downloads (0–7 MB/s, repeated `[FAILED]`) |
| Library pins and sklearn kNN | scanpy 1.11.5 incompatible with scipy ≥ 1.14 and umap-learn 0.6; compatibility only, same algorithm |
| ResNet-50 weights from the HF timm mirror | pytorch.org connection resets; identical IMAGENET1K_V2 weights, all keys matched |
| "Other" threshold 0.5 → 0.3 | The only G1 tweak the master prompt allows; TENX191 Other 32% → 22% |
| Alpha grid extended to 1e9 | 1e6 was a boundary for 89% of residual fits; superset of the original grid |
| Shares only where ceiling ≥ 0.2 | Avoids dividing by near-zero ceilings; those genes are noise_limited anyway |
| G3(c) accepted at 91.4% / 80.0% (user's choice, option A of A/B/C) | Not a bug: non-transferable image signal; the alternatives were a nested residual fit or stopping |
| Not used: multiple GPUs | IT policy (one GPU) and no benefit; the remaining compute is CPU ridge |
| Robustness round pre-registered: thresholds fixed before results (R² ≥ 0.3; headline in ≥3/4 patients; within CI) | Avoids tuning after seeing the numbers |
| Floor held-out R² at 0 for shares (user-specified) | Single-patient R² outliers (NPY1R Δmorph 2.24 in P4) dominated the shares |
| Joint typing kept as robustness only | It did not strengthen composition; a batch effect remains |

## 8. Limitations & future work
**Limitations**
- **Inconsistent cell-type labels** (TENX202 "Mast" cluster, TENX191 large "Other"). They make composition look weaker, so **Δmorph is an upper bound**.
- **Interferon** has only 4 genes.
- **Scale:** one study, one section per patient (n = 4), 280-gene panel.
- **Thinning:** Poisson/binomial, not negative-binomial.
- **Encoders:** gated foundation models not tested.
- **G3(c)** below 95%.

**Future work**
- UNI2-h and Virchow2.
- All 11 sections.
- TCGA bulk validation.
- Visium HD.
- Negative-binomial thinning.
- Harmonised cell typing with better batch correction, or segmentation-aware adipocyte/stromal typing.

## 9. Reproducibility
- **Server:** `<KCDH server>`, project `$DH607_STORAGE/dh607`.
- **Python:** `.venv/bin/python` (3.11). Key pins: scipy 1.13.1, numba 0.61.2, scanpy 1.11.5, umap-learn 0.6.0.
- **Caches:**
  - `HF_HOME=$DH607_STORAGE/cache/huggingface`
  - `TORCH_HOME=$DH607_STORAGE/cache/torch`
- **Commands:**
  ```
  .venv/bin/python code/02_build_spots.py TENX191        # per sample, OMP_NUM_THREADS=16
  CUDA_VISIBLE_DEVICES=0 .venv/bin/python code/03_extract_features.py TENX191 TENX195 TENX199 TENX202
  .venv/bin/python code/04_decompose.py resnet50 --samples TENX191 TENX195 TENX199 TENX202   # OMP_NUM_THREADS=24, ~5 min
  .venv/bin/python code/04_decompose.py phikon2  --samples TENX191 TENX195 TENX199 TENX202   # ~1 min
  .venv/bin/python code/05_figures.py
  ```
- **Seed:** `config.SEED = 0`, used for thinning, Leiden and the shuffle.
- **Logs:** `logs/` (`spots_*.log`, `decompose_*_final2.log`, `figures_final.log`).
- **Copies on the Mac:** results in `~/Projects/DH607-core/results/core/`, code in `~/Projects/DH607-core/pipeline/`.
