
# Persistent Homology via Finite Topological Spaces (FTS)

This repository contains a **research prototype** implementing persistent homology
via **finite topological spaces (FTS)**, as introduced in the paper

> **Persistent Homology via Finite Topological Spaces**  
> Selçuk Kayacan  
> arXiv:2512.23348  
> https://arxiv.org/abs/2512.23348

The goal of this project is to provide an alternative viewpoint on persistent
homology in which filtrations are built from **topologies on a fixed finite set**,
rather than from nested simplicial complexes. The emphasis is on understanding
**how dense regions of data interact across scales**, not on geometric shape
reconstruction.

---

## High-level description of the method

Let (X,d) be a finite metric space.

1. **Local density scale**  
   Fix an integer `k_density = k`.  
   For each point x ∈ X, let r_k(x) be the distance from x to its
   k-th nearest neighbor. This determines a local scale at x.

2. **Dense points (anchors)**  
   A simple density proxy is defined by  
   ρ(x) = 1 / (r_k(x) + ε).  
   For a threshold value τ, points with ρ(x) ≥ τ are regarded
   as *dense* and are called **anchors**.

3. **Anchor neighborhoods**  
   Each anchor a generates a neighborhood  
   N(a) = { y ∈ X | d(a,y) ≤ λ r_k(a) },  
   where `lambda` is a fixed expansion factor. These neighborhoods are adaptive:
   their sizes depend on local density.

4. **Finite topology**  
   The anchor neighborhoods generate a topology on X. Points that are not covered
   by any anchor neighborhood are handled by controlled fallback rules to avoid
   degenerate behavior.

5. **Associated combinatorial complex**  
   From the resulting finite topological space, a simplicial complex is built
   (using maximal regions as vertices). Homology is computed on the 2-skeleton of
   this complex.

6. **Filtration**  
   Repeating the above construction for a sequence of thresholds τ yields
   a filtration of topologies  
   T1 ⊇ T2 ⊇ ⋯ ⊇ Tm,  
   and hence a persistence summary.


---

## Interpreting the output tables

Each example produces a table indexed by filtration level, displaying:

- Betti numbers β0 and β1
- The number of maximal regions
- Basic combinatorial data of the associated complex
- Diagnostic information on coverage and monotonicity


---

## Repository structure

```text
.
├── fts_pipeline.sage
├── LICENSE
├── README.md
└── examples/
    ├── disjoint_clusters.sage
    ├── figure8_blobs.sage
    ├── six_blobs_on_circle.sage
    ├── six_bumps_on_ring_band.sage
    ├── tree_like_clusters.sage
    ├── figures/
    │   └── *.png
    └── tables/
        └── *.txt
```

---

## Running the code

All examples are run from the **root directory**.

For example:
```bash
sage examples/figure8_blobs.sage
```

Each example:
- generates a dataset plot (saved in `examples/figures/`),
- runs the density-anchor filtration,
- prints a summary table.

---

## Status

⚠️ **Research prototype**

- The code prioritizes clarity and mathematical structure over performance.
- It is not optimized for large-scale data.
- Parameter defaults reflect a conservative, density-driven interpretation.

---

## Citation

If you use this code, please cite:

```bibtex
@article{FTS2025,
  title   = {Persistent Homology via Finite Topological Spaces},
  author  = {Kayacan, Selçuk},
  journal = {arXiv preprint arXiv:2512.23348},
  year    = {2025}
}
```

---

## License

Released under the **MIT License**.
See the LICENSE file for details.
