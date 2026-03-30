# Repository map

## Core directories

- `src/paper1/` — curated reusable numerical code
- `scripts/` — runnable manuscript-facing entry points
- `results/` — reference outputs used to anchor the manuscript analyses
- `docs/` — concise reader-facing guidance
- `archive/` — minimal legacy material still required by some wrapper scripts

## Main manuscript-facing scripts

- `scripts/01_cube_solver_basics.py`
- `scripts/02_geometry_center_readout.py`
- `scripts/03_tetrahedron_control.py`
- `scripts/04_hexahedron_shear_study.py`
- `scripts/05_center_emergence.py`
- `scripts/06_primitive_cubic_supercell.py`
- `scripts/07_fcc_supercell_readout.py`

## Extended support scripts

- `scripts/main_strand/nullmodel_blindtests/`
- `scripts/side_branches/face_diagnosis/`
- `scripts/side_branches/field_information/`
- `scripts/side_branches/ordering_principle/`

## Reader workflow

1. Read `docs/manuscript_link.md`
2. Check `docs/reproducibility_status.md`
3. Run `scripts/reproduce_public_minimum.sh`
4. Inspect relevant reference outputs under `results/`
