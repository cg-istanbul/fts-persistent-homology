# Persistent Homology via Finite Topological Spaces (FTS)

This repository contains a **research prototype** implementing persistent homology via **finite topological spaces (FTS)**, as introduced in the paper

> **Persistent Homology via Finite Topological Spaces**  
> arXiv:2512.23348  
> https://arxiv.org/abs/2512.23348

The code provides an alternative to classical Rips filtrations by constructing
topologies, specialization posets, and crosscut complexes from finite metric data.

---

## Overview

Given a finite point cloud \( X \subset \mathbb{R}^d \) and a scale parameter \( \alpha \),

1. We construct **Rips-like open sets** as maximal \( \alpha \)-cliques.

2. These open sets generate a **finite topology** on \( X \).

3. The topology is reduced to a **\(T_0\)-space** via Kolmogorov identification.

4. The resulting **specialization poset** is used to construct a **crosscut complex**.

5. Homology is computed on the **2-skeleton** of the crosscut complex.

6. Repeating this over increasing \( \alpha \) yields a **persistent homology summary**.

This approach naturally supports **local analysis** and **cycle localization**, which are harder to access in standard simplicial filtrations.

---

## Repository Structure

```text
.
├── fts_basis.sage        # FTS construction, minimal basis, T0 quotient, crosscut complex, homology computation
├── run_filtration.sage   # End-to-end example (filtration + table)
├── LICENSE
└── README.md

---

## Requirements

- SageMath (tested with Sage ≥ 9.x)
- Python packages (usually bundled with Sage):
  - `numpy`

This code is intended to be run using:
```bash
sage run_filtration.sage

---

## Quick Start

Run a complete filtration experiment on synthetic data:

- Generate a synthetic point cloud (e.g. blobs on a circle)

- Construct FTS filtrations for multiple α-values

- Compute β0 and β1

- Output a barcode-like summary table

Example output:

=== Persistent Homology Summary (FTS / Crosscut) ===
alpha   | #V | #S | β0 | β1
---------------------------------------------
0.3970  | 26 |110 | 1  | 0
1.7961  | 30 |590 | 1  | 1
2.1459  | 30 |4525| 1  | 0

---

## Notes on Interpretation

- The filtration **is not simplicial by inclusion**, so matrix-reduction-based persistent homology does not directly apply.

- Homology is computed **one scale at a time**.

- The appearance of β1 = 1 corresponds to **macroscopic cycles** emerging from overlapping finite open sets.

- The framework naturally supports tracking the **support of homology classes** inside the finite space (future extension).

---

## Status

⚠️ **Research prototype**

- Code prioritizes clarity over performance

- Not optimized for large datasets

- API may change

Contributions, experiments, and discussions are welcome.

---

## Citation

If you use this code, please cite:

@article{FTS2025,
  title={Persistent Homology via Finite Topological Spaces},
  author={Kayacan, Selçuk},
  journal={arXiv preprint arXiv:2512.23348},
  year={2025}
}

---

## License

This project is released under the **MIT License**.
See the LICENSE file for details.








