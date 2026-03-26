# Repository map

## Core numerical kernel

- `src/paper1/cube_fem.py`
- `src/paper1/hexa_fem.py`
- `src/paper1/tetra_fem.py`
- `src/paper1/octant.py`
- `src/paper1/center_emergence.py`

## Main-strand scripts

- `scripts/01_cube_solver_basics.py`
- `scripts/02_geometry_center_readout.py`
- `scripts/03_tetrahedron_control.py`
- `scripts/04_hexahedron_shear_study.py`
- `scripts/05_center_emergence.py`
- `scripts/06_primitive_cubic_supercell.py`
- `scripts/07_fcc_supercell_readout.py`

## Main-strand hardening block

- `scripts/main_strand/nullmodel_blindtests/01_reference_digest.py`
- `scripts/main_strand/nullmodel_blindtests/02_label_shuffle_permutation_tests.py`
- `scripts/main_strand/nullmodel_blindtests/03_feature_ablation_tests_reduced.py`
- `scripts/main_strand/nullmodel_blindtests/04_field_shuffle_cube_null_tests.py`
- `scripts/main_strand/nullmodel_blindtests/05_field_shuffle_supercell_reconstructed.py`

Reference tables:

- `results/main_strand/nullmodel_blindtests/reference/`

## Side branch: face diagnosis

- `scripts/side_branches/face_diagnosis/01_face_information_first_test.py`
- `scripts/side_branches/face_diagnosis/02_face_information_extended.py`
- `scripts/side_branches/face_diagnosis/03_face_patch_scan_cube.py`
- `scripts/side_branches/face_diagnosis/04_face_patch_scan_supercell.py`
- `results/side_branches/face_diagnosis/reference/`

## Side branch: field information

- `scripts/side_branches/field_information/01_field_information_first_test.py`
- `scripts/side_branches/field_information/02_field_information_supercell.py`
- `results/side_branches/field_information/reference/`

## Side branch: explicit ordering principle

- `scripts/side_branches/ordering_principle/01_geometry_continuation_study.py`
- `scripts/side_branches/ordering_principle/02_geometry_phase_map.py`
- `results/side_branches/ordering_principle/reference/`

## Historical provenance retained publicly

- `archive/original_scripts/`
- `archive/main_strand/nullmodel_blindtests/original_scripts/`
- `archive/side_branches/*/original_scripts/`

Large chat-export logs are intentionally not part of the public release.
