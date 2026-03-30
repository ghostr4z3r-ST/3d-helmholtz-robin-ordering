# 3d-helmholtz-robin-ordering

Companion reproducibility repository for the manuscript **Reproducible Low-Mode Ordering in a 3D Helmholtz--Robin Eigenproblem on Cube-Like and Lattice-Linked Geometries** by **Steven Trümpert**.

This repository is intentionally focused on one task: enabling external readers to reproduce the numerical analyses underlying the manuscript. It is **not** intended to serve as a public research diary or a complete reconstruction archive of every intermediate project phase.

## What this repository contains

This repository provides:

- curated source code in `src/paper1/` for the reusable numerical core,
- runnable analysis scripts in `scripts/`,
- reference outputs in `results/` for the analyses reported in the manuscript,
- a minimal archive of legacy source files in `archive/` only where current wrapper scripts still depend on them,
- concise reader-facing documentation in `docs/`.

## Manuscript scope covered here

The repository covers the numerical analyses corresponding to the manuscript sections on:

- low-mode ordering in the cubic kernel,
- emergent center readout and Helmholtz/Robin role separation,
- geometric hardening across cube-like geometries,
- lattice-linked readout and family signatures,
- null models and selective hardening,
- face-diagnosis, field-information, and ordering-principle support analyses where they are used as manuscript support rather than as standalone narratives.

## Quickstart

Use Python 3.10+.

```bash
python -m pip install -e .
```

Run the main manuscript-facing reproduction paths:

```bash
./scripts/reproduce_core.sh
./scripts/reproduce_extended.sh
./scripts/reproduce_public_minimum.sh
```

Additional branch-specific runs:

```bash
./scripts/main_strand/nullmodel_blindtests/reproduce_nullmodel_digest.sh
./scripts/side_branches/face_diagnosis/reproduce_face_diagnosis_v1.sh
./scripts/side_branches/field_information/reproduce_field_information_v1.sh
```

## Start here

For external readers, the most useful files are:

- `docs/manuscript_link.md`
- `docs/repo_map.md`
- `docs/reproducibility_status.md`
- `docs/manuscript_figure_map.md`
- `docs/forschungshistorie.md`

## Repository and manuscript linkage

This repository belongs to the manuscript:

> **Reproducible Low-Mode Ordering in a 3D Helmholtz--Robin Eigenproblem on Cube-Like and Lattice-Linked Geometries**

The manuscript should be read as the primary narrative document. This repository supplies the code, analysis scripts, and reference outputs needed to reproduce the numerical results reported there.

## Public release record

- **GitHub repository:** `https://github.com/ghostr4z3r-ST/3d-helmholtz-robin-ordering`
- **Zenodo DOI (latest released archive, v1.1.0):** `[10.5281/zenodo.19236767](https://doi.org/10.5281/zenodo.19335898)`

## Licensing

- **Code:** MIT (`LICENSE`)
- **Documentation, CSV files, and figures:** CC BY 4.0 (`LICENSE-data`)
