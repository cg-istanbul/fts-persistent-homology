# Persistent Homology via Finite Topological Spaces (FTS)

This repository contains a **research prototype** implementing persistent homology
via **finite topological spaces (FTS)**, as introduced in the paper

> **Persistent Homology via Finite Topological Spaces**  
> Selçuk Kayacan  
> arXiv:2512.23348

The goal of this project is to provide an alternative viewpoint on persistent
homology in which filtrations are built from **topologies on a fixed finite set**,
rather than from nested simplicial complexes. The emphasis is on understanding
**how dense regions of data interact across density scales**, not on geometric
shape recovery.

---

## Conceptual overview

Let $(X,d)$ be a finite metric space (typically a point cloud).

Instead of building a filtration of **nested simplicial complexes**, we construct
a filtration of **topologies on the same underlying set**:

$$
\mathcal T_1 \supseteq \mathcal T_2 \supseteq \dots \supseteq \mathcal T_m,
$$

where each topology is **coarser** than the previous one. The identity map on $X$
is continuous
$$
(X,\mathcal T_i) \to (X,\mathcal T_j)\quad (i<j),
$$
which makes the construction **functorial**.

For each topology $(X,\mathcal T_i)$:

1. Form the associated **$T_0$-space** by identifying points with identical
   minimal open neighborhoods.
2. Encode the $T_0$-space as a **finite poset**, ordered by set inclusion of minimal open neighborhoods.
3. Build the **order complex** of this poset and compute homology on its
   **2-skeleton** (sufficient for $\beta_0$ and $\beta_1$).

Although the simplicial complexes are **not nested**, the simplicial maps between
steps are induced by the identity on $X$, and the resulting persistence module is
well-defined.

---

## Density-guided topology construction (implemented)

The current implementation supports a single topology mode:

```
density_guided
```

The goal is to track how **regions of varying density interact** as the density
threshold changes.

### Local density scale

Fix an integer parameter `k_density = k`.

For each point $x\in X$, let $r_k(x)$ be the distance from $x$ to its
$k$-th nearest neighbor. This defines a **local scale** at each point.

A simple density proxy is defined by

$$
\rho(x)=\frac{1}{r_k(x)+\varepsilon}.
$$

### Anchors and neighborhoods

For a threshold $\tau$, points with $\rho(x)\ge \tau$ are called **anchors**.

Each anchor $a$ generates a neighborhood

$$
N(a)=\{y\in X \mid d(a,y)\le \lambda\cdot r_k(a)\},
$$

where `lambda_scale = λ > 1`.

Anchor neighborhoods generate the topology at that step. If some points are not
covered by any anchor neighborhood, the pipeline introduces **controlled fallback
generators** so that the filtration remains well-behaved and the induced maps
remain meaningful.

---

## Filtration parameters

The following parameters are explicitly used in the current code and appear in the
tables:

- `m` (int, ≥2): number of filtration steps.
- `k_density` (int, ≥1): k-th nearest neighbor used to define $r_k(\cdot)$.
- `lambda_scale` (float, >1): neighborhood expansion factor $\lambda$.
- `schedule` (string): threshold schedule (currently `"quantile"`).
- `tau_q_start` (float in [0,1]): start quantile for $\tau$.
- `tau_q_end` (float in [0,1]): end quantile for $\tau$.
- `tau_q_exponent` (float ≥1): nonlinear spacing of quantiles.

---

## Provided outputs (inspect before running)

This repository intentionally includes precomputed outputs so that readers can
inspect results **without running anything**:

- `examples/tables/*_tables.txt`  
  Filtration summary tables (parameters + per-step diagnostics + $\beta_0,\beta_1$).

- `examples/barcodes/*.svg`  
  Rendered barcode plots.

- `examples/barcodes/*.jsonl`  
  Machine-readable barcodes (one JSON object per bar).

- `examples/figures/*.png`  
  Dataset figures for synthetic examples.

- `examples/mapper/*.png`  
  Mapper visualizations (only for context; not part of the PH pipeline).

---

## Repository layout (current)

```
.
├── examples/
│   ├── barcodes/                 # *.svg and *.jsonl barcodes (β0/β1)
│   ├── figures/                  # synthetic dataset figures (*.png)
│   ├── mapper/                   # mapper figures (*.png), optional context
│   ├── tables/                   # *_tables.txt summary tables
│   ├── synthetic/                # synthetic example scripts (*.sage)
│   └── real/                     # real example scripts + dataset download
│       ├── download_and_preprocess_coil20.py
│       ├── download_and_preprocess_paul15.py
│       ├── coil20_n1000.sage
│       └── paul15_n1000.sage
└── scripts/
    ├── algebra/                  # H0/H1 over the 2-skeleton
    ├── pipeline/                 # filtration runner + validation
    ├── report/                   # tables + barcodes + svg + mapper helpers
    ├── topology/                 # anchors, schedule, T0-poset, order complex
    └── util/
```

---

## Real datasets: download and preprocess

Real datasets are downloaded and preprocessed by scripts in `examples/real/`.
For example, to download COIL20:

```bash
python examples/real/download_and_preprocess_coil20.py
```

This will create a `data/` folder in the repository root and place the processed
dataset files there. After that, you can run the corresponding Sage example, e.g.

```bash
sage examples/real/coil20_n1000.sage
```

Similarly for Paul15:

```bash
python examples/real/download_and_preprocess_paul15.py
sage examples/real/paul15_n1000.sage
```

---

## Interpreting the summary tables (`examples/tables/*_tables.txt`)

Each table file contains:

1. The command that produced it.
2. A **parameter summary line**, e.g.
   ```
   m = 25 | k_density = 7 | lambda_scale = 1.650000 | schedule = quantile | tau_q_start = 0.080000 | tau_q_end = 0.950000 | tau_q_exponent = 2.000000
   ```
3. A per-step table with header like:
   ```
     i           τ      β₀      β₁   cyc/vert    cyc/edge   anchors  covered  uncovered      |V|       |E|         |T|
   ```

### Column-by-column legend

- **`i`**  
  Filtration step index (1..m). Step 1 is the initial “finest” stage.

- **`τ`**  
  Density threshold at this step (printed as `-` at step 1).

- **`β₀`**  
  Number of connected components of the order-complex 2-skeleton $K_i$.  
  Interpretation: the number of **distinct interacting regions** at this density scale.

- **`β₁`**  
  Number of 1-dimensional homology classes of $K_i$.  
  Interpretation in this pipeline: **cyclic interaction patterns between regions**,
  not geometric holes in the ambient space.

- **`cyc/vert`**  
  Cycle density normalized by the number of vertices: $\beta_1/|V|$.

- **`cyc/edge`**  
  Cycle density normalized by the number of edges: $\beta_1/|E|$.

- **`anchors`**  
  Number of anchors (points with $\rho(x)\ge \tau$).

- **`covered`**  
  Number of points covered by at least one anchor neighborhood.

- **`uncovered`**  
  Number of points not covered by anchor neighborhoods (handled by fallback generators).

- **`|V|`**  
  Number of vertices of $K_i$. These are the $T_0$-classes (poset elements).
  This is expected to be **non-increasing** along the filtration.

- **`|E|`**  
  Number of edges in the 1-skeleton of $K_i$.

- **`|T|`**  
  Number of triangles (2-simplices) in the 2-skeleton of $K_i$.

### How to read β0 and β1 “persistence” here

- **β₀ behavior**: regions merge as $\tau$ increases (topologies coarsen).  
  Long-lived β₀ features indicate region separations that are stable across many
  density thresholds.

- **β₁ behavior**: β₁ features represent stable **multi-region interaction cycles**
  (in the order-complex model).  
  This is not “shape recovery”: persistent β₁ does **not** claim a geometric loop in
  the ambient space. It indicates that the interaction pattern among dense regions
  consistently contains a cycle across density scales.

---

## Interpreting the barcode files (`examples/barcodes/`)

For each example and each dimension (β0 / β1) we provide:

- `*_beta0.svg`, `*_beta1.svg`  
  Barcode plots.

- `*_beta0.jsonl`, `*_beta1.jsonl`  
  One JSON object per bar. 

### What the JSONL contains (high level)

Each bar includes:

- dimension (`dim`),
- birth/death step and birth/death $\tau$,
- `meta.birth` with dimension-specific metadata.

For β0, `meta.birth` includes information such as component support size and a
canonical representative point.

For β1, `meta.birth` may include a sparse description of a chosen cycle representative
at its birth step (cycle edges in the order complex 2-skeleton), which can be useful
for interpretation and debugging.

---

## Mapper (optional context only)

Each example script contains a **commented** block that can be uncommented to run
Mapper / Mapper-comparison code.

Mapper is included **only for context**. This repository implements a persistent
homology pipeline; Mapper is not.

---

## Running the examples

Run examples from the repository root. Synthetic examples, e.g.

```bash
sage examples/synthetic/figure8_blobs.sage
```

Real examples require dataset download first (see above).

---

## Status

⚠️ Research prototype.

- Emphasis is on mathematical structure and interpretability.
- Not optimized for large-scale datasets.
- Interfaces and outputs may evolve.

---

## License

Released under the MIT License.
See the LICENSE file for details.
