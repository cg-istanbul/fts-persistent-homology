# scripts/report/outputs.sage

import os, sys

def current_example_name():
    base = os.path.basename(sys.argv[0])
    return base.replace(".sage.py", "").replace(".sage", "")

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path

def out_dir(kind):
    # kind: "figures", "barcodes", "mapper", "meta"
    root = os.path.join("examples", kind)
    ensure_dir(root)
    return root

def scatter_png_path(name):
    return os.path.join(out_dir("figures"), f"{name}.png")

def mapper_png_path(name):
    return os.path.join(out_dir("mapper"), f"{name}_mapper.png")

