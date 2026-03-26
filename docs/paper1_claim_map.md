# Paper 1 claim map

This document freezes the intended **claim structure** for Paper 1 and maps each claim to the
supporting repository material.

It is meant to keep the manuscript narrow:
- **main text** should carry only the load-bearing claims,
- **supplement** may carry compact secondary support,
- **repo only** contains full reproducibility, extended diagnostics, and provenance.

## Framing rule

Paper 1 should be written as a **numerical / geometric evidence paper**.
The target is **not** to claim a complete physical theory, but to document a reproducible,
artifact-hardened ordering structure of the 3D Helmholtz-Robin problem on cube-like and lattice-linked geometries.

---

## Claim C0 — Minimal scope / non-overclaiming statement

**Claim wording**

> The paper studies a geometric 3D Helmholtz-Robin eigenproblem and reports a reproducible ordering structure inside that model class. The paper does not claim a complete physical theory; it claims a robust numerical-geometric evidence framework.

**Paper placement**
- main text: Introduction + Discussion

**Support in repo**
- `README.md`
- `docs/provenance_and_reconstructions.md`
- `docs/reproducibility_status.md`

**Purpose**
- protects the paper from overclaiming
- aligns the paper with what the repo actually demonstrates

---

## Claim C1 — The 3D Helmholtz-Robin kernel carries an ordered mode landscape

**Claim wording**

> The 3D Helmholtz-Robin setup on the cube already exhibits a nontrivial ordered mode structure rather than an undifferentiated numerical spectrum.

**Paper status**
- main text

**Best paper support**
- Figure 1: basic spectrum / beta scan / convergence

**Primary scripts**
- `scripts/01_cube_solver_basics.py`
- core modules under `src/paper1/`

**Primary repo evidence**
- runtime outputs from `scripts/reproduce_core.sh`
- optional reference context from `results/reference/`

**Reason this belongs in main text**
- this is the entry claim; everything else depends on it

**Do not overstate as**
- “final physical spectrum”
- “proof of universal ordering”

---

## Claim C2 — The center signature emerges from the solution and is not manually inserted

**Claim wording**

> The characteristic center reading is an emergent property of the computed modes and is not introduced by an explicit center term.

**Paper status**
- main text

**Best paper support**
- Figure 2: cube / quader / hexa center-emergence comparison

**Primary scripts**
- `scripts/02_geometry_center_readout.py`
- `scripts/05_center_emergence.py`

**Primary repo evidence**
- reproducible outputs from the scripts above
- `src/paper1/octant.py`
- `src/paper1/center_emergence.py`

**Secondary support**
- `docs/paper_figure_map.md`

**Reason this belongs in main text**
- this is one of the sharpest distinctions between “emergent reading” and “built-in ansatz”

---

## Claim C3 — Helmholtz carries the volumetric structure; Robin modulates spectral selection rather than creating the structure

**Claim wording**

> The volumetric mode structure is borne by the Helmholtz problem, while the Robin term acts primarily as a selector or reweighter of spectral position rather than as the origin of the structure itself.

**Paper status**
- main text

**Best paper support**
- Figure 1 or Figure 3 depending on final layout
- compact beta comparison in Results

**Primary scripts**
- `scripts/01_cube_solver_basics.py`
- `scripts/05_center_emergence.py`

**Primary repo evidence**
- public runs from `scripts/reproduce_public_minimum.sh`
- core FEM modules under `src/paper1/`

**Secondary support**
- side-branch ordering synthesis in
  `results/side_branches/ordering_principle/reference/geometry_phase_map_summary.csv`
- `results/side_branches/ordering_principle/reference/geometry_phase_map_compact.csv`

**Reason this belongs in main text**
- this is the conceptual center of Paper 1

**Critical wording discipline**
- say “supports the interpretation that …” or “is consistent with …”
- do not imply analytic proof beyond the numerical evidence shown

---

## Claim C4 — Cube-like / hexahedral geometries support mixed volumetric ordering more robustly than anisotropic contrast geometries

**Claim wording**

> The ordering structure is not confined to the perfect cube; it persists across cube-like hexahedral deformations and weakens or changes character under anisotropic contrast geometries.

**Paper status**
- main text

**Best paper support**
- Figure 3: tetrahedron / quader / hexahedral contrast and shear hardening

**Primary scripts**
- `scripts/03_tetrahedron_control.py`
- `scripts/04_hexahedron_shear_study.py`
- `scripts/05_center_emergence.py`

**Primary repo evidence**
- script outputs from the three scripts above

**Secondary support**
- `results/side_branches/ordering_principle/reference/geometry_continuation_study_first_test.csv`
- `results/side_branches/ordering_principle/reference/geometry_phase_map_full.csv`
- `results/side_branches/ordering_principle/reference/geometry_phase_map_overview.png`

**Reason this belongs in main text**
- this is the main artifact-hardening geometry block

---

## Claim C5 — The ordering extends from single-cell geometry to lattice-linked / finite-q family readout

**Claim wording**

> The observed ordering is not restricted to isolated-cell reading; it extends into lattice-linked supercell diagnostics, where family-level distinctions become readable in finite-q style observables.

**Paper status**
- main text, but kept compact

**Best paper support**
- Figure 4: primitive cubic vs. fcc readout

**Primary scripts**
- `scripts/06_primitive_cubic_supercell.py`
- `scripts/07_fcc_supercell_readout.py`

**Primary repo evidence**
- `results/reference/primitive_cubic_supercell_ordered_3x3x3.csv`
- `results/reference/cubic_family_signature_matrix.csv`
- `results/reference/cubic_family_structure_factor_proxy_results.csv`
- `results/reference/cubic_family_resolution_robustness.csv`
- `results/reference/cubic_family_physical_readout_results.csv`

**Secondary support**
- `results/side_branches/ordering_principle/reference/geometry_phase_map_summary.csv`

**Reason this belongs in main text**
- it shows that the structure scales beyond the isolated cubic cell and starts to acquire a family / signature language

**Limit to maintain in paper**
- do not let this turn into a full crystallographic paper inside Paper 1
- keep it as a support block for the main ordering claim

---

## Claim C6 — Natural face observables are weaker than volume / lattice observables; the ordering is volumetric-primary

**Claim wording**

> Within the tested diagnostics, natural face-based observables remain systematically weaker than volume or lattice observables, indicating that the ordering is volumetric-primary rather than face-native.

**Paper status**
- supplement or compressed main-text paragraph

**Best paper support**
- one compact summary table or one sentence with supplement pointer

**Primary scripts**
- `scripts/08_face_information_sidebranch.py`
- `scripts/side_branches/face_diagnosis/02_face_information_extended.py`
- `scripts/side_branches/face_diagnosis/03_face_patch_scan_cube.py`
- `scripts/side_branches/face_diagnosis/04_face_patch_scan_supercell.py`

**Primary repo evidence**
- `results/side_branches/face_diagnosis/reference/face_information_first_test.csv`
- `results/side_branches/face_diagnosis/reference/face_information_extended_test.csv`
- `results/side_branches/face_diagnosis/reference/face_information_extended_cube_robustness.csv`
- `results/side_branches/face_diagnosis/reference/face_patch_scan_first_test.csv`
- `results/side_branches/face_diagnosis/reference/face_patch_scan_supercell_test.csv`

**Reason this is not a main load-bearing figure**
- it strengthens the reading of the kernel, but the paper can stand without giving the full branch center stage

---

## Claim C7 — Full 3D ordering is not reducible to naive “maximal delocalization”; field-information diagnostics show a more specific structure

**Claim wording**

> Field-information diagnostics indicate that the relevant 3D ordering is not captured by a naive “most spread out” interpretation; it is better described through a joint reading of effective support, entropy-style measures, and anisotropy / family structure.

**Paper status**
- supplement or compressed main-text paragraph

**Primary scripts**
- `scripts/side_branches/field_information/01_field_information_first_test.py`
- `scripts/side_branches/field_information/02_field_information_supercell.py`

**Primary repo evidence**
- `results/side_branches/field_information/reference/field_information_distribution_first_test.csv`
- `results/side_branches/field_information/reference/field_information_distribution_supercell_test.csv`

**Reason this is secondary in the manuscript**
- useful for interpretation discipline
- not needed as a central paper figure if the manuscript must stay lean

---

## Claim C8 — The explicit ordering principle can be summarized geometrically as symmetry-supported mixed ordering that tips under anisotropy and is modulated by Robin coupling

**Claim wording**

> A compact geometric ordering principle is supported by the continuation and phase-map diagnostics: symmetry-equivalent 3D directions sustain mixed volumetric ordering, anisotropy pushes the system toward axis-bound families, and Robin coupling modulates where this transition appears.

**Paper status**
- main text discussion paragraph or supplement figure

**Primary scripts**
- `scripts/side_branches/ordering_principle/01_geometry_continuation_study.py`
- `scripts/side_branches/ordering_principle/02_geometry_phase_map.py`

**Primary repo evidence**
- `results/side_branches/ordering_principle/reference/geometry_continuation_study_first_test.csv`
- `results/side_branches/ordering_principle/reference/geometry_phase_map_summary.csv`
- `results/side_branches/ordering_principle/reference/geometry_phase_map_compact.csv`
- `results/side_branches/ordering_principle/reference/geometry_phase_map_overview.png`

**Reason this may deserve a short main-text appearance**
- it gives the paper a concise conceptual synthesis
- but it should not swallow the main numerical narrative

---

## Claim C9 — The ordering survives hardening selectively under null models; it is therefore not an obvious label, feature, or shuffle artifact

**Claim wording**

> The strongest classification and ordering signals remain selectively stable under label shuffling, field shuffling, and feature ablation, while weaker channels collapse as they should. This selective survival argues against the main findings being obvious numerical or feature-construction artifacts.

**Paper status**
- main text, near the end of Results

**Best paper support**
- Table 1 or Figure 5: null-model / blind-test digest

**Primary scripts**
- `scripts/main_strand/nullmodel_blindtests/01_reference_digest.py`
- `scripts/main_strand/nullmodel_blindtests/02_label_shuffle_permutation_tests.py`
- `scripts/main_strand/nullmodel_blindtests/03_feature_ablation_tests_reduced.py`
- `scripts/main_strand/nullmodel_blindtests/04_field_shuffle_cube_null_tests.py`
- `scripts/main_strand/nullmodel_blindtests/05_field_shuffle_supercell_reconstructed.py`

**Primary repo evidence**
- `results/main_strand/nullmodel_blindtests/reference/label_shuffle_permutation_tests.csv`
- `results/main_strand/nullmodel_blindtests/reference/feature_ablation_tests_reduced.csv`
- `results/main_strand/nullmodel_blindtests/reference/feature_ablation_summary_reduced.csv`
- `results/main_strand/nullmodel_blindtests/reference/field_shuffle_null_tests_50.csv`
- `results/main_strand/nullmodel_blindtests/reference/field_shuffle_supercell_null_tests.csv`
- `results/main_strand/nullmodel_blindtests/reference/supercell_qdominated_field_shuffle.csv`
- `results/main_strand/nullmodel_blindtests/reference_digest.md`

**Why this is a load-bearing final claim**
- it is the main anti-artifact closure for Paper 1
- it justifies why the earlier ordering claims should be taken seriously

**Important wording discipline**
- say “not an obvious artifact” or “unlikely to be reducible to the tested artifacts”
- do not claim that every possible artifact class has been eliminated

---

## Claim C10 — Paper 1 closes as a geometric evidence framework, not as a complete downstream application paper

**Claim wording**

> Paper 1 establishes a geometric and reproducible evidence base. It closes before the detector / double-slit application and therefore serves as the numerical-geometric backbone for later work rather than as the full application paper itself.

**Paper status**
- Discussion / Conclusion

**Support in repo**
- overall repository structure
- `docs/repo_map.md`
- `docs/reproducibility_status.md`
- `docs/public_release_checklist.md`

**Purpose**
- preserves scope
- prevents the manuscript from collapsing into a “first everything paper”

---

# Suggested main-text spine

If the manuscript must stay tight, the main text can be built primarily on:
- C0
- C1
- C2
- C3
- C4
- C5
- C9
- C10

And use these as supplement / compressed support:
- C6
- C7
- C8

---

# Suggested figure / table loadout

- **Figure 1**: basic cube spectral ordering / beta scan
- **Figure 2**: center emergence across core geometries
- **Figure 3**: tetrahedron / quader / hexahedral contrast
- **Figure 4**: primitive cubic vs. fcc family readout
- **Table 1**: null-model and blind-test digest
- **Supplement Figure S1**: face-diagnosis summary
- **Supplement Figure S2**: field-information summary
- **Supplement Figure S3**: ordering-principle phase map

---

# Scope guardrails

Things that should stay **out of the central Paper 1 claim set** unless absolutely necessary:
- broad physical interpretation beyond the tested model class
- downstream detector / double-slit material
- full research-history narration
- every auxiliary diagnostic in full width
- any wording that implies formal proof where the evidence is numerical

