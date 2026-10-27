# Decomposing the Predictability of Spatial Gene Expression from Histopathology Images

This repository holds the Phase 1 baseline: H&E patches → frozen ImageNet ResNet-50 features → Ridge regression → expression of 50 genes on the HEST-Benchmark IDC task.

Run all commands from the repository root.

## Setup

```bash
conda env create -f environment.yml
conda activate computational_multiomics
```

## Get the data

```bash
hf download MahmoodLab/hest-bench --repo-type dataset --include "IDC/*" --exclude "IDC/old/*" --local-dir data/hest-bench
```



## Run

Open `baseline.ipynb` with the `computational_multiomics` kernel and run all cells.

## Data credit

- Dataset: [HEST-Benchmark, IDC task](https://huggingface.co/datasets/MahmoodLab/hest-bench/tree/main/IDC) on Hugging Face
- Source: HEST — Jaume et al., *HEST-1k: A Dataset for Spatial Transcriptomics and Histology Image Analysis*, NeurIPS 2024 ([paper](https://arxiv.org/abs/2406.16192), [code](https://github.com/mahmoodlab/HEST), [documentation](https://hest.readthedocs.io/))
- License: [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). The data are not included in this repository.
