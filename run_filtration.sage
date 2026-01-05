# run_filtration.sage

import numpy as np

load("fts_basis.sage")
load("homology.sage")

np.random.seed(0)

def blobs_on_circle(num_blobs, points_per_blob, radius, noise):
    X = []
    for i in range(num_blobs):
        theta = 2*np.pi*i/num_blobs
        center = np.array([radius*np.cos(theta), radius*np.sin(theta)])
        blob = np.random.normal(center, noise, (points_per_blob,2))
        X.append(blob)
    return np.vstack(X)

X = blobs_on_circle(
    num_blobs=6,
    points_per_blob=6,
    radius=1.0,
    noise=0.3
)

def alpha_range(X, num):
    dists = []
    n = len(X)
    for i in range(n):
        for j in range(i+1, n):
            dists.append(np.linalg.norm(X[i]-X[j]))
    a, b = min(dists), max(dists)
    return [a + (b-a)*k/(num-1) for k in range(num)]

alphas = alpha_range(X, 10)

print("\n=== Persistent Homology Summary (FTS / Crosscut) ===")
print("alpha\t| #V\t| #S\t| β0\t| β1")
print("-"*45)

for alpha in alphas:
    B0, poset, crosscut = fts_pipeline(X, alpha)
    betti = homology_ranks_from_simplices(crosscut, 1)

    print(
        f"{alpha:.4f}\t| "
        f"{len(B0):>2}\t| "
        f"{len(crosscut):>3}\t| "
        f"{betti[0]:>2}\t| "
        f"{betti[1]:>2}"
    )
