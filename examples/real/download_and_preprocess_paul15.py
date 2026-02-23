# examples/real/download_and_preprocess_paul15.py

import os
import json
import numpy as np
import pandas as pd
import scanpy as sc

# ----------------------------
# Configuration (edit if you want)
# ----------------------------
N_TARGET = 1000
SEED = 0
N_HVG = 500
N_PCS = 30

# ----------------------------
# Paths (repo-root relative)
# ----------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

OUT_DIR = os.path.join(ROOT, "data", "paul15")
OUT_X = os.path.join(OUT_DIR, f"paul15_n{N_TARGET}.csv")
OUT_META = os.path.join(OUT_DIR, f"paul15_n{N_TARGET}_meta.csv")
OUT_MANIFEST = os.path.join(OUT_DIR, f"paul15_n{N_TARGET}_manifest.json")

os.makedirs(OUT_DIR, exist_ok=True)
rng = np.random.default_rng(SEED)

# ----------------------------
# 1) Load dataset (Paul15; ~2730 cells)
# ----------------------------
adata = sc.datasets.paul15()

# Keep raw counts in a layer (safer for HVG)
# (Paul15 comes as counts in .X in most scanpy versions; we store anyway.)
adata.layers["counts"] = adata.X.copy()

# ----------------------------
# 2) Normalize + log (standard)
# ----------------------------
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ----------------------------
# 3) HVG selection using counts layer (Seurat v3 style)
# ----------------------------
# Note: seurat_v3 expects counts; passing layer="counts" is the safest choice.
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_HVG,
    flavor="seurat_v3",
    layer="counts"
)

# ----------------------------
# 4) Subset to HVGs
# ----------------------------
if "highly_variable" not in adata.var:
    raise RuntimeError("HVG selection failed: adata.var['highly_variable'] missing.")
adata = adata[:, adata.var["highly_variable"]].copy()

# ----------------------------
# 5) Scale (on log-normalized HVGs)
# ----------------------------
sc.pp.scale(adata, max_value=10)

# ----------------------------
# 6) PCA
# ----------------------------
# random_state makes runs more reproducible across machines/versions.
sc.tl.pca(adata, n_comps=N_PCS, svd_solver="arpack", random_state=SEED)

# ----------------------------
# 7) Downsample to N_TARGET cells (stratified if possible)
# ----------------------------
n_obs = int(adata.n_obs)
if N_TARGET >= n_obs:
    idx = np.arange(n_obs)
    cluster_key = None
else:
    # Prefer a known Paul15 cluster column if present, else fall back to common ones.
    cluster_key = None
    for k in ["paul15_clusters", "clusters", "leiden", "louvain"]:
        if k in adata.obs:
            cluster_key = k
            break

    if cluster_key is None:
        idx = rng.choice(n_obs, size=N_TARGET, replace=False)
    else:
        groups = adata.obs[cluster_key].astype(str).to_numpy()
        uniq, counts = np.unique(groups, return_counts=True)

        # Proportional allocation
        raw = counts * (N_TARGET / n_obs)
        base = np.floor(raw).astype(int)
        base = np.maximum(base, 1)  # keep at least 1 per group (when possible)

        # Adjust to exact total using largest remainders
        diff = int(N_TARGET - base.sum())
        rema = raw - np.floor(raw)
        order = np.argsort(-rema)  # descending

        j = 0
        while diff != 0 and j < 10_000:
            g = order[j % len(order)]
            if diff > 0:
                base[g] += 1
                diff -= 1
            else:
                if base[g] > 1:
                    base[g] -= 1
                    diff += 1
            j += 1

        # Sample within each group
        idx_list = []
        for g, take in zip(uniq, base):
            members = np.flatnonzero(groups == g)
            take = min(int(take), len(members))
            if take > 0:
                idx_list.append(rng.choice(members, size=take, replace=False))

        idx = np.concatenate(idx_list) if idx_list else np.array([], dtype=int)

        # Final fixup to exact N_TARGET
        idx = np.unique(idx)
        if len(idx) > N_TARGET:
            idx = rng.choice(idx, size=N_TARGET, replace=False)
        elif len(idx) < N_TARGET:
            remaining = np.setdiff1d(np.arange(n_obs), idx, assume_unique=False)
            extra = rng.choice(remaining, size=(N_TARGET - len(idx)), replace=False)
            idx = np.concatenate([idx, extra])

adata_sub = adata[idx].copy()

# ----------------------------
# 8) Export
# ----------------------------
X = np.asarray(adata_sub.obsm["X_pca"][:, :N_PCS], dtype=float)  # (N_TARGET, N_PCS)

# Meta: include cluster label if available + original obs_name
meta = pd.DataFrame({
    "obs_name": adata_sub.obs_names.astype(str)
})

if cluster_key is not None:
    meta["label"] = adata_sub.obs[cluster_key].astype(str).to_numpy()
else:
    meta["label"] = "NA"

pd.DataFrame(X).to_csv(OUT_X, index=False, header=False)
meta.to_csv(OUT_META, index=False)

manifest = {
    "dataset": "scanpy.datasets.paul15()",
    "N_TARGET": N_TARGET,
    "SEED": SEED,
    "N_HVG": N_HVG,
    "N_PCS": N_PCS,
    "cluster_key": cluster_key,
    "output_X": os.path.relpath(OUT_X, ROOT),
    "output_meta": os.path.relpath(OUT_META, ROOT),
    "scanpy_version": getattr(sc, "__version__", "unknown"),
}
with open(OUT_MANIFEST, "w", encoding="utf-8") as f:
    json.dump(manifest, f, indent=2)

print(f"Wrote:\n  {OUT_X}\n  {OUT_META}\n  {OUT_MANIFEST}\nShape: {X.shape}")

