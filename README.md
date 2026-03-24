# Elastic Cavitation Interactions

FEniCS codes for the paper:

> **"The effect of interactions on elastic cavitation"**
> Ali Saeedi, S. Chockalingam, Mrityunjay Kothari
> *arXiv:submit/7386717 [cond-mat.soft], March 2026*
> Department of Mechanical Engineering, University of New Hampshire & Purdue University

---

## Overview

This repository contains all simulation code used to study how proximity to a **rigid interface** or a **neighboring identical cavity** shifts the critical cavitation pressure in an incompressible neo-Hookean solid under remote hydrostatic tension.

Two interaction problems are studied:

| Problem | Bottom boundary | Script |
|---------|----------------|--------|
| Cavity near a rigid interface | Fixed BC (u₁ = u₂ = 0) | `BVP_Bulk_2Cavity_fixedBC_2.py` |
| Two interacting identical cavities | Symmetry BC (u₂ = 0 only) | `BVP_Bulk_2Cavity_symmBC_2.py` |

Key results from the paper:
- **Cavity–interface**: cavitation pressure increases monotonically from the bulk value p̄c,bulk = 2.5μ as the cavity approaches the interface, reaching p̄c,int ≈ 3.5μ at d/R ≈ 1.
- **Cavity–cavity**: cavitation pressure varies non-monotonically with inter-cavity spacing, peaking at p̄c ≈ 2.8μ around d/R ≈ 5.7.

---

## Repository Structure

```
elastic-cavitation-interactions/
│
├── geo_file_creator_2.py          # generates Gmsh .geo geometry file for a given d/R
├── 2d_gmsh_convert.py             # converts Gmsh .msh to XDMF format for FEniCS
├── BulkandSurfaceElasticity_edited.py  # constitutive model class (neo-Hookean + surface energy)
├── BVP_Bulk_2Cavity_fixedBC_2.py  # BVP solver: cavity-interface interaction (rigid wall)
├── BVP_Bulk_2Cavity_symmBC_2.py   # BVP solver: cavity-cavity interaction (symmetry plane)
├── run_sim.sh                     # single-run pipeline: mesh → convert → solve
├── batch_run.sh                   # batch runner over all d/R values from the paper
│
├── meshes/                        # generated mesh files (.geo, .msh, .xdmf) — created at runtime
└── results/                       # simulation output files (.mat, .xdmf) — created at runtime
```

---

## Dependencies

| Package | Purpose |
|---------|---------|
| [FEniCS](https://fenicsproject.org/) | finite element solver (legacy version, tested with FEniCS 2019.1) |
| [Gmsh](https://gmsh.info/) | mesh generation from .geo files |
| [meshio](https://github.com/nschloe/meshio) | mesh format conversion (.msh → .xdmf) |
| [NumPy](https://numpy.org/) | numerical array operations |
| [SciPy](https://scipy.org/) | nonlinear solve for mesh grading + .mat file output |

## Usage

### Running a single simulation

```bash
./run_sim.sh <y_c> <mesh_quality> <Delta_Vbarf> <steps> <BC>
```

| Argument | Description |
|----------|-------------|
| `y_c` | Dimensionless cavity distance d/R (e.g. `1.10`) |
| `mesh_quality` | Element size level: `Coarse`, `Medium`, `Fine`, or `Finest` |
| `Delta_Vbarf` | Target volume expansion (V−V₀)/V₀ (e.g. `3374`) |
| `steps` | Number of incremental load steps (e.g. `15000`) |
| `BC` | `fixed` = cavity-interface problem; `symm` = cavity-cavity problem |

**Example — cavity near a rigid interface at d/R = 1.10:**
```bash
./run_sim.sh 1.10 Finest 3374 15000 fixed
```

**Example — two interacting cavities at d/R = 5.35:**
```bash
./run_sim.sh 5.35 Finest 3374 15000 symm
```

`run_sim.sh` executes the following pipeline automatically:
1. Generate the `.geo` geometry file (`geo_file_creator_2.py`)
2. Move it to `meshes/` and run Gmsh to produce a `.msh` file
3. Convert `.msh` to `.xdmf` format (`2d_gmsh_convert.py`)
4. Launch the FEniCS solver (`BVP_Bulk_2Cavity_fixedBC_2.py` or `BVP_Bulk_2Cavity_symmBC_2.py`)

### Running all d/R values from the paper

```bash
./batch_run.sh
```

This loops over all 27 d/R values used in the paper and runs each through the full pipeline. By default it uses `fixed` BC (cavity-interface problem). To run the cavity-cavity problem instead, change `fixed` to `symm` in `batch_run.sh`.

---

## Output

Each simulation produces a `.mat` file in `results/` with the following columns:

| Column | Variable |
|--------|----------|
| 0 | Remote displacement applied at top boundary |
| 1 | Current cavity volume V |
| 2 | Bulk strain energy |
| 3 | Deformed cavity surface area |
| 4 | Remote pressure p̄ (average normal traction on outer boundary) |
| 5 | Inner pressure (average normal traction on cavity surface, should be ≈ 0) |

If `save_xdmf=1` is passed to the BVP scripts, displacement and stress fields are also saved as `.xdmf` files viewable in [ParaView](https://www.paraview.org/).

The cavitation pressure p̄c for each d/R value is extracted from the `.mat` files by fitting the pressure–volume curve using Laguerre polynomials and extrapolating to V → ∞ (see Section 4.2 of the paper).

---

## Simulation Details

- **Material model**: incompressible neo-Hookean solid, strain energy W = (μ/2)(Ī₁ − 3)
- **Incompressibility**: enforced via a perturbed Lagrangian mixed u-p formulation
- **Element type**: Taylor-Hood (P2/P1) — quadratic displacement, linear pressure
- **Solver**: Newton-Raphson with MUMPS direct linear solver
- **Loading**: displacement-controlled at the outer boundary (equivalent to remote hydrostatic tension)
- **Domain size**: 100R × (d + 100R) (verified to approximate infinite medium behavior)
- **Mesh**: spatially adaptive, refined near cavity surface (minimum element size R/100 at Finest quality)
- **Simulation software**: [FEniCS](https://fenicsproject.org/) with mesh generation via [Gmsh](https://gmsh.info/)

---

## Acknowledgements

Many thanks to **Chockalingam Senthilnathan** for providing the FEniCS code for spherical growth, which served as the foundation for this research. The weak-form setup in FEniCS parallels the implementation in:

> C. Senthilnathan, *Understanding the mechanics of growth: A large deformation theory for coupled swelling-growth and morphogenesis of soft biological systems*, Ph.D. thesis, Massachusetts Institute of Technology (2024).

---

## Citation

If you use this code in your research, please cite:

```bibtex
@article{saeedi2026cavitation,
  title   = {The effect of interactions on elastic cavitation},
  author  = {Saeedi, Ali and Chockalingam, S. and Kothari, Mrityunjay},
  journal = {arXiv preprint arXiv:submit/7386717},
  year    = {2026}
}
```

---

## Contact

- **Mrityunjay Kothari** — mrityunjay.kothari@unh.edu  
  Department of Mechanical Engineering, University of New Hampshire
