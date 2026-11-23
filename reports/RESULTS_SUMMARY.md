# DH607 — Results summary (for slides / poster)

*Sources:*
- *Category counts, ceilings, R² and programme values: `results/core/decomposition_{resnet50,phikon2}.csv` and `summary_*.txt` (final run, 2026-11-24).*
- *Share-based numbers (comp/img ratio, Δmorph, morphology genes): R² floored at 0 per held-out sample, a pre-specified fix for unstable held-out R² (`results/robustness/floored/`).*
- *Robustness: `results/robustness/`.*

*Design: 280 genes, 4 patients, leave-one-patient-out. 95% CIs come from 1000 bootstrap resamples of the patients.*

## One-sentence takeaway
On a deep Xenium panel the noise ceiling is high (median 0.707), so failure is **not** mostly noise. What H&E models predict across patients is mostly **cell-type composition**. Morphology on top of composition adds almost nothing, and in all 4 patients. The few exceptions are adipose and smooth-muscle genes, which most likely reflect cell types that the composition features miss.

## Headline findings

1. **Noise explains very little of the failure.**
   - Only **3 of 280 genes (1.1%)** are noise-limited (ceiling < 0.2).
   - The median ceiling is **0.707**, rising with expression: quartile medians 0.588, 0.702, 0.706, 0.819.
   - Yet **226 of 280 (ResNet-50)** and **181 of 280 (Phikon-v2)** genes are *not predicted*.
   - Figure: `fig1_noise_resnet50.png` (and `_phikon2`).
2. **Apparent image success is mostly composition.**
   - Among genes the image model "predicts" (R²_img > 0.05), the median R²_comp / R²_img is **1.771 [95% CI 1.557–2.003]** for ResNet-50 (n = 105) and **1.118 [0.998–1.179]** for Phikon-v2 (n = 142).
   - So composition alone roughly matches (Phikon-v2, CI lower bound just below 1) or clearly beats (ResNet-50) the image.
   - Figure: `fig2_decomposition_resnet50_floored.png` (and `_phikon2`).
3. **Morphology beyond composition is rare, small, and not consistent across patients.**
   - Median Δmorph is about 0: 9.2e-6 for ResNet-50 and 4.0e-5 for Phikon-v2.
   - Genes with morph_share ≥ 0.10 and R²_full > 0:
     - ResNet-50: **2 [0–5]**, ADIPOQ and LPL;
     - Phikon-v2: **5 [2–9]**, ADIPOQ, LPL, MYH11, ADH1B and APOC1.
   - Only 2 (ResNet-50) / 5 (Phikon-v2) genes show image signal beyond composition, all adipose/lipid or smooth-muscle markers (ADIPOQ, LPL, ADH1B, APOC1; MYH11) — cell types our composition features capture poorly (also unstable in fig5). Only MYH11 replicates in ≥3/4 patients. No cell-state or programme gene shows reproducible morphology signal (`results/robustness/floored/replicated_morph_genes.csv`).
   - Figure: `fig_headline_ci_floored.png`.
4. **Phikon-v2 mostly sees more composition, not more morphology.**
   - It "predicts" more genes than ResNet-50: R²_img > 0.05 in 100 vs 59 (core, unfloored; 142 vs 105 when floored), composition_explained 90 vs 51.
   - Its morphology share stays small, and Δmorph agrees only weakly between encoders (r = **0.19**).
   - Figure: `fig4_encoders.png`.

## Robustness
- **The headline is robust.** It holds in 4/4 held-out patients for both encoders, and the image encodes the major cell types. The benchmark condition concerns absolute image performance, not the claim.
- **The image encodes composition.** The R² for predicting cell-type fraction from the image is above 0.3 for epithelial (0.42 / 0.53) and T cells (0.37 / 0.53) with both encoders (ResNet-50 / Phikon-v2). Fibroblast (0.40) and endothelial (0.38) pass with Phikon-v2 only. → `fig5_image_to_composition.png`
- **Benchmark sanity.** The image models score below published random-split benchmarks (top-50 r 0.28 / 0.35 vs ~0.47 / 0.53), as expected under leave-one-patient-out; encoder ranking agrees. (Source: github.com/mahmoodlab/HEST.)
- **Consistency.** The headline holds in **4 of 4** held-out patients for both encoders. The per-patient comp/img ratio is 1.52–2.12 (ResNet-50) and 0.98–1.17 (Phikon-v2).
- **Joint cell typing** (Harmony across samples) did NOT strengthen composition: 91 vs 104 genes with R²_comp > 0.05, with a residual batch effect (cluster 0 = 95% TENX191). Per the pre-set rule, per-sample labels stay the main result. Morphology genes go 0 → 0 (ResNet-50) and 2 → 1 (Phikon-v2): **the morphology conclusion is unchanged under an alternative joint labelling.** → `fig6_joint_vs_persample.png`

## Programme verdicts (ResNet-50 / Phikon-v2)
*Core (unfloored) values. With R² floored at 0, Hypoxia's R²_comp becomes 0 and Interferon's ResNet-50 R²_img becomes 0.099; the verdicts are unchanged.*


| Programme | Genes on panel | Verdict | Key number |
|---|---|---|---|
| Proliferation | 20 / 20 | Composition | R²_comp 0.265 (28% of ceiling 0.957); Δmorph −0.008 / 0.006 |
| EMT | 15 / 15 | Composition | R²_comp 0.729 (75% of ceiling 0.971); Δmorph 0.013 / 0.027 |
| Interferon | 4 / 4 | Composition. ResNet-50's image alone is below 0.1×ceiling, so its category is not_predicted. | R²_comp 0.320 (35% of ceiling 0.926); Δmorph −0.000 / −0.004 |
| Hypoxia | 8 / 8 | Not predicted (not noise) | ceiling 0.833, but R²_comp −0.072 and R²_img −0.037 / −0.086 |

No programme is morphology-visible: every morph_share is below 0.03. → `fig3_programmes.png`

## ResNet-50 vs Phikon-v2
Phikon-v2 raises image-only R² (100 vs 59 genes above 0.05, core; 142 vs 105 floored), but almost entirely by reproducing composition: R²_comp/R²_img is 1.118 vs 1.771 (floored).
Its extra morphology is limited to a few adipose and stromal genes, and Δmorph agrees only weakly between encoders (r = 0.19).

## Methods (5 bullets)
- Data: HEST-1k Xenium breast, 280 genes. 4 patients, one representative section each. 100 µm pseudo-spots with 224 px H&E patches at 0.5 µm/px.
- Transcripts with qv ≥ 20 are binomially thinned 50/25/25 into A (cell typing → composition), B1 and B2. The target is B = B1 + B2.
- Noise ceiling = Spearman–Brown of corr(B1, B2), per sample.
- Three ridge models: composition; image (ResNet-50 or Phikon-v2); image fitted on the composition residuals. Leave-one-patient-out, with alpha chosen by inner leave-one-patient-out.
- Δmorph = R²_full − R²_comp. Shares are fractions of the ceiling (genes with ceiling ≥ 0.2), with held-out R² floored at 0. Categories and the programme table use core values.

## Limitations (4 bullets)
- Cell-type labels are inconsistent across samples (the TENX202 "Mast" cluster; 22% "Other" in TENX191), so **Δmorph is an upper bound**. Joint typing did not fix this; there is a residual batch effect.
- One study, one section per patient (n = 4), 280-gene panel. The interferon programme has only 4 genes.
- R²_full ≥ R²_comp − 0.01 holds for 91.4% (ResNet-50) and 80.0% (Phikon-v2) of genes. Absolute image performance is below random-split benchmarks.
- Thinning is Poisson (binomial), not negative-binomial. Gated foundation models (UNI2-h, Virchow2, H-optimus, CONCH) were not tested.

## Figures and suggested titles
| File | Suggested title |
|---|---|
| `results/core/figures/fig1_noise_resnet50.png` | "Failure is not noise: high ceilings, low image R²" |
| `results/robustness/figures/fig2_decomposition_resnet50_floored.png` / `_phikon2_floored` | "Where each gene's predictability comes from" |
| `results/core/figures/fig3_programmes.png` | "Biological programmes: composition, not morphology" |
| `results/core/figures/fig4_encoders.png` | "Pathology foundation model: more composition, little more morphology" |
| `results/robustness/figures/fig5_image_to_composition.png` | "H&E image features encode the major cell types" |
| `results/robustness/figures/fig_headline_ci_floored.png` | "The headline holds in every patient" |
| `results/robustness/figures/fig6_joint_vs_persample.png` | "Joint cell typing: composition not stronger, morphology unchanged" |
