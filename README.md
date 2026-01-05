# Finite Topological Spaces (FTS) Filtration and Persistence

This repository implements a **finite-topological-space–based filtration** for point cloud data and computes homological invariants via **crosscut complexes**.

The code accompanies ongoing research on persistence theory for finite topological spaces (FTS), providing a computational pipeline parallel to (but distinct from) Vietoris–Rips persistence.

---

## Overview

Given a finite metric space (e.g. a point cloud in Euclidean space), the pipeline:

1. Constructs a family of **FTS topologies** parameterized by a scale parameter `α`
2. Reduces each topology to its **T₀ quotient**
3. Forms the associated **specialization poset**
4. Builds the **crosscut complex** (2-skeleton)
5. Computes **homology ranks** (β₀, β₁)
6. Outputs a **barcode-like table** indexed by `α`

Unlike classical persistent homology, the simplicial complexes here are **not nested by inclusion**, so homology is computed *independently at each scale*.

---

## Mathematical Background (Informal)

- Each scale `α` defines a family of open sets via **Rips-like maximal cliques**
- These open sets generate a finite topology on the point set
- T₀-identification collapses indistinguishable points
- The specialization order defines a poset
- The **crosscut complex** captures the homotopy type of the finite space
- Homology of the crosscut complex represents the homology of the FTS

This construction supports:
- Stability results (proved in the accompanying paper)
- Local and global topological analysis
- Extensions beyond Euclidean point clouds

---

## Repository Structure

```text
.
├── fts_basis.sage        # FTS construction, minimal basis, T0 quotient, crosscut
├── homology.sage         # Homology computation via Sage
├── run_filtration.sage   # End-to-end example (filtration + table)
├── LICENSE
└── README.md

---

Requirements

- SageMath (tested with Sage ≥ 9.x)
- Python packages (usually bundled with Sage):
  - `numpy`

This code is intended to be run using:
```bash
sage run_filtration.sage

---

## Example: Running a Filtration

- The script run_filtration.sage demonstrates:
- Synthetic datasets (e.g. blobs on a circle)
- Automatic selection of α values
- Computation of Betti numbers β₀ and β₁

Typical output:

=== Persistent Homology Summary (FTS / Crosscut) ===
alpha   | #V | #S | β0 | β1
---------------------------------------------
0.3970  | 27 | 63 | 13 |  0
1.7961  | 30 | 590|  1 |  1

---

## Notes on Interpretation

- The number of vertices corresponds to T₀-equivalence classes
- Large α values may increase combinatorial complexity
- The appearance of β₁ reflects nontrivial loops in the induced topology
- Unlike Rips persistence, **support and localization of cycles** are meaningful and recoverable in FTS

---

## Future Directions

- Enriched barcodes tracking support of homology classes
- Localized FTS persistence
- Visualization tools
- Extensions to general finite metric spaces
- Comparison with Rips and Čech persistence

---

## Citation

If you use this code in your research, please cite:

**Selçuk Kayacan**,  
*Persistent Homology via Finite Topological Spaces*,  
arXiv:2512.23348 (2025).  
https://arxiv.org/abs/2512.23348

---

## License

This project is released under the **MIT License**.
See the LICENSE file for details.

---

## Disclaimer

This is **experimental research code**.
Interfaces and mathematical definitions may change as theory develops.

Contributions, discussions, and feedback are welcome.

---

## Authors

Selçuk Kayacan






