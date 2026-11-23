"""
Step 03 - image features for every H&E patch (one patch per spot).

Two frozen encoders, no fine-tuning:
  resnet50 : ImageNet ResNet-50, global average pool -> 2048 numbers per patch
  phikon2  : Phikon-v2 (Owkin), ViT-L trained with DINOv2 on pathology tiles,
             CLS token -> 1024 numbers per patch. Open weights, no gated access.

Uses ONE GPU (CUDA_VISIBLE_DEVICES). Saves features/<model>/<sample>.npz with the
patch barcodes so they can be matched to spots later.

Usage: python 03_extract_features.py TENX199 TENX191 ...
"""
import os
import sys
import time

import h5py
import numpy as np
import torch
import torchvision

import config

IMAGENET_MEAN = torch.tensor([0.485, 0.456, 0.406]).view(1, 3, 1, 1)
IMAGENET_STD = torch.tensor([0.229, 0.224, 0.225]).view(1, 3, 1, 1)


def load_encoder(name, device):
    if name == "resnet50":
        net = torchvision.models.resnet50(weights=torchvision.models.ResNet50_Weights.IMAGENET1K_V2)
        net.fc = torch.nn.Identity()                      # keep the pooled 2048-d vector
        encode = lambda x: net(x)
    elif name == "phikon2":
        from transformers import AutoModel
        net = AutoModel.from_pretrained("owkin/phikon-v2")
        encode = lambda x: net(pixel_values=x).last_hidden_state[:, 0]   # CLS token
    else:
        raise ValueError(name)
    net.eval().to(device)
    return encode


@torch.inference_mode()
def extract(sample_id, model_name, encode, device):
    out_dir = os.path.join(config.FEATURE_DIR, model_name)
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{sample_id}.npz")
    if os.path.exists(out_path):
        print(f"skip {model_name}/{sample_id} (exists)")
        return

    with h5py.File(os.path.join(config.HEST_DIR, "patches", f"{sample_id}.h5"), "r") as h5:
        images = h5["img"]                                 # N x 224 x 224 x 3, uint8, 0.5 um/px
        barcodes = np.array([b[0].decode() for b in h5["barcode"][:]])
        features = []
        for start in range(0, len(barcodes), config.BATCH_SIZE):
            batch = torch.from_numpy(images[start:start + config.BATCH_SIZE]).permute(0, 3, 1, 2)
            batch = (batch.float() / 255 - IMAGENET_MEAN) / IMAGENET_STD
            with torch.autocast("cuda", dtype=torch.float16):
                features.append(encode(batch.to(device)).float().cpu().numpy())
    features = np.concatenate(features)
    np.savez_compressed(out_path, barcodes=barcodes, features=features)
    print(f"{model_name}/{sample_id}: {features.shape}", flush=True)


def main(sample_ids):
    device = "cuda"
    for model_name in ["resnet50", "phikon2"]:
        start = time.time()
        encode = load_encoder(model_name, device)
        for sample_id in sample_ids:
            extract(sample_id, model_name, encode, device)
        print(f"{model_name} done in {time.time() - start:.0f}s", flush=True)
        torch.cuda.empty_cache()


if __name__ == "__main__":
    main(sys.argv[1:])
