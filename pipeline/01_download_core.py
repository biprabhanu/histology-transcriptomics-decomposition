"""
Download the DH607 core dataset from HEST-1k (MahmoodLab/hest).

Dataset: "Xenium v1 Human Breast FFPE with Biomarkers & Housekeeping Genes Custom Panel"
         11 sections from 4 patients (one study, one panel, one pixel size).
Order:   one section per patient first, so the pipeline can start early.
Files:   metadata/, st/, patches/, transcripts/  (no WSIs, no segmentation needed)
"""
import os
import time
from huggingface_hub import snapshot_download

DATA_DIR = os.environ.get("DH607_STORAGE", "/path/to/storage") + "/dh607/data/hest"
FOLDERS = ["metadata", "st", "patches", "transcripts"]

FIRST = ["TENX199", "TENX191", "TENX202", "TENX195"]          # one per patient
REST = ["TENX197", "TENX200", "TENX196", "TENX193", "TENX201", "TENX192", "TENX198"]


def download_sample(sample_id, max_tries=8):
    patterns = [f"{folder}/{sample_id}[_.]*" for folder in FOLDERS]
    for attempt in range(1, max_tries + 1):
        try:
            snapshot_download("MahmoodLab/hest", repo_type="dataset",
                              local_dir=DATA_DIR, allow_patterns=patterns, max_workers=4)
            print(f"[done] {sample_id}", flush=True)
            return
        except Exception as error:          # network drops happen on this server
            print(f"[retry {attempt}] {sample_id}: {error}", flush=True)
            time.sleep(20 * attempt)
    print(f"[FAILED] {sample_id}", flush=True)


if __name__ == "__main__":
    for sample_id in FIRST + REST:
        download_sample(sample_id)
    print("ALL DOWNLOADS FINISHED", flush=True)
