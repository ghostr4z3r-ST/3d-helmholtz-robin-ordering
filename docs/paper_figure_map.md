# Draft paper figure / table map

This is a manuscript-assembly aid for Paper 1.
It is intentionally a **draft map**, not a claim that the final figure numbering is already fixed.

## Main strand

### Figure A — cube spectrum / β scan / convergence
- script: `scripts/01_cube_solver_basics.py`
- outputs: `results/reference/` and user-generated run outputs from `scripts/reproduce_core.sh`
- paper role: minimal numerical kernel and first spectral ordering

### Figure B — center emergence and geometry comparison
- scripts: `scripts/02_geometry_center_readout.py`, `scripts/05_center_emergence.py`
- paper role: center signature emerges from the solution rather than being added externally

### Figure C — tetrahedron control and hexahedral contrast
- scripts: `scripts/03_tetrahedron_control.py`, `scripts/04_hexahedron_shear_study.py`
- paper role: artifact hardening across contrast geometries

### Figure D — primitive cubic vs. fcc supercell readouts
- scripts: `scripts/06_primitive_cubic_supercell.py`, `scripts/07_fcc_supercell_readout.py`
- reference outputs: `results/reference/primitive_cubic_supercell_ordered_3x3x3.csv`
- paper role: finite-q / family-level ordering readout

### Table E — null-model and blind-test summary
- directory: `results/main_strand/nullmodel_blindtests/reference/`
- helper: `scripts/main_strand/nullmodel_blindtests/01_reference_digest.py`
- paper role: selective robustness against shuffle and ablation baselines

## Side branches / supplement candidates

### Face diagnosis
- scripts: `scripts/side_branches/face_diagnosis/*`
- results: `results/side_branches/face_diagnosis/reference/`
- role: shows that natural face observables remain weaker than volume / lattice observables

### Field information
- scripts: `scripts/side_branches/field_information/*`
- results: `results/side_branches/field_information/reference/`
- role: supplements the volumetric reading with entropy / anisotropy style diagnostics

### Explicit ordering principle
- scripts: `scripts/side_branches/ordering_principle/*`
- results: `results/side_branches/ordering_principle/reference/`
- role: compact geometric synthesis, especially the phase-map view
