"""
Run step 04 on the joint-clustering spots (data/spots_joint) -> results/robustness/joint.

Usage: python 04b_decompose_joint.py resnet50     (any 04 arguments are passed through)
"""
import os
import sys
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import config

config.SPOT_DIR = os.path.join(config.PROJECT, "data", "spots_joint")
config.RESULT_DIR = os.path.join(config.PROJECT, "results", "robustness", "joint")

spec = importlib.util.spec_from_file_location("decompose", os.path.join(HERE, "04_decompose.py"))
s04 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s04)
assert s04.config is config and s04.config.SPOT_DIR.endswith("spots_joint")

if __name__ == "__main__":
    s04.main()
