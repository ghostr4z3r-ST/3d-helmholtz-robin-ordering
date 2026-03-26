# 3d-helmholtz-robin-ordering

Companion reproducibility repository for the paper **Artifact-Hardened Ordering in a 3D Helmholtz-Robin Eigenproblem on Cube-Like and Lattice-Linked Geometries** by **Steven Trümpert**.

Public reproducibility repository for **Paper 1** on the numerical 3D Helmholtz–Robin core.

This repository is meant to do one job well: provide a clear, auditable backbone for the later paper without forcing the manuscript itself to carry every technical detail.

## Scope

This release covers the full Paper-1 package:

- the **main strand** from the cube solver to the supercell readouts,
- the **null-model / blind-test block**,
- the **face-diagnosis** side branch,
- the **field-information** side branch,
- and the **explicit ordering-principle** side branch.

The repository is intentionally split into a small number of layers:

- `src/paper1/` contains curated reusable code for the core numerical pieces.
- `scripts/` contains runnable entry points and wrappers.
- `results/` contains versioned reference outputs that anchor the paper claims.
- `archive/` contains recovered historical source files that are still useful for provenance.
- `docs/` explains what is direct, historical, reconstructed, and intended for the paper.

## Reproducibility levels

This repo distinguishes three levels of material:

1. **Curated core**  
   Clean modules and scripts intended for direct use.

2. **Recovered historical scripts**  
   Original files from the research history, wrapped so that old local paths are not required.

3. **Reconstructed support files**  
   Only used where a historical source survived merely as a stub. These files are explicitly marked.

## Quickstart

Create a Python environment with Python 3.10+ and install the package:

```bash
python -m pip install -e .
```

Run the main curated core:

```bash
./scripts/reproduce_core.sh
./scripts/reproduce_extended.sh
```

Run the null-model digest and branch-specific historical wrappers:

```bash
./scripts/main_strand/nullmodel_blindtests/reproduce_nullmodel_digest.sh
./scripts/side_branches/face_diagnosis/reproduce_face_diagnosis_v1.sh
./scripts/side_branches/field_information/reproduce_field_information_v1.sh
```

## Paper-facing documentation

Start with these files:

- `docs/reproducibility_status.md`
- `docs/repo_map.md`
- `docs/paper_figure_map.md`
- `docs/provenance_and_reconstructions.md`

## Public provenance policy

The public repository includes code, reference outputs, and concise provenance notes.
The large chat-export archive used during reconstruction is **not bundled** in this public release. It is kept outside the public repo so that the GitHub presentation stays focused on the scientific material.

## Licensing

This repository uses a **split license model**:

- **Code** is licensed under **MIT** (`LICENSE`).
- **Documentation, CSV result tables, and images** are licensed under **CC BY 4.0** (`LICENSE-data`).

That split is deliberate: Creative Commons is reasonable for non-code research material, but it is usually not the best choice for executable source code.

## Citation

A starter `CITATION.cff` is included.
Before the public Zenodo release, fill in the final author metadata and repository DOI.

## Suggested citation pairing

Use this repository together with the companion arXiv preprint of the same title. Once the Zenodo release is created, add the DOI badge and cite the tagged repository release in addition to the paper.
