# Benchmark sanity check (naive_r, non-programme genes)
Ours (median [IQR]), all 280 genes / top-50 by mean_count: ResNet-50 0.197 [0.047-0.377] / 0.279 [0.087-0.447]; Phikon-v2 0.297 [0.113-0.488] / 0.353 [0.113-0.568].
Published HEST-Benchmark IDC (Pearson, top-50 HVG, ridge on PCA-256 features; README snapshot 08.30.24): ResNet50 ~0.473, Phikon (v1) ~0.533, CONCH ~0.536, UNI ~0.572.
Source: https://github.com/mahmoodlab/HEST (README benchmark table; Jaume et al. 2024, arXiv 2406.16192). Values read from README snippets, not re-verified against the NeurIPS paper tables; I could not confirm the IDC task is Xenium-based.
Comparability: loose only. Ours is leave-one-patient-out within one study on 100 um pseudo-spots with our own gene set (not top-50 HVG), and their splits, patch size and regressor differ; ours are median-over-genes, theirs are mean-over-genes.
Verdict: our top-50 values (0.28 ResNet-50, 0.35 Phikon-v2) sit below the published ~0.47-0.53, but the ordering (foundation model > ResNet-50) matches.
