# Paper 1 arXiv section map

This document translates the frozen claim set in `docs/paper1_claim_map.md` into a compact
**arXiv-oriented manuscript plan**.

The goal is a manuscript that is:
- readable without the repo,
- narrow enough for a first arXiv paper,
- strong enough to stand on the repository as its reproducibility backbone,
- disciplined in scope: **numerical-geometric evidence paper**, not full physical theory.

---

## Working paper title candidates

### Conservative title
**Ordering Structure in a 3D Helmholtz-Robin Eigenproblem on Cube-Like and Lattice-Linked Geometries**

### Slightly stronger but still safe
**Artifact-Hardened Ordering in a 3D Helmholtz-Robin Eigenproblem: From Cube-Like Domains to Lattice-Linked Readout**

### Very compact title
**Geometric Ordering in a 3D Helmholtz-Robin Eigenproblem**

Final paper title: **Artifact-Hardened Ordering in a 3D Helmholtz-Robin Eigenproblem on Cube-Like and Lattice-Linked Geometries**.

---

## Target manuscript profile

### Suggested length
- **main text**: ~12–18 pages before references
- **figures/tables**: 4 main figures + 1 main table
- **appendix**: compact, not a second paper

### Main-text load-bearing claims
- C0
- C1
- C2
- C3
- C4
- C5
- C9
- C10

### Secondary claims
- C6
- C7
- C8

These should appear as short support, appendix support, or repo support.

---

## Abstract shape

The abstract should do exactly four things:

1. state the model class
2. state the core finding of ordered structure
3. state the anti-artifact hardening
4. state the scope limit

### Abstract skeleton

> We study a three-dimensional Helmholtz eigenproblem with Robin boundary coupling on cube-like and lattice-linked geometries. Across a sequence of core geometries, supercell readouts, and diagnostic branches, we find a reproducible ordering structure that is not well described as an undifferentiated numerical spectrum. The results support a separation in which the Helmholtz problem carries the volumetric mode structure while Robin coupling modulates spectral selection. The main signals persist selectively under geometry variation, null-model shuffling, and feature-ablation tests, arguing against a straightforward discretization or feature-construction artifact. The paper is presented as a numerical-geometric evidence framework and does not claim a complete physical theory.

This should be tightened later, but its logic should remain unchanged.

---

## Main-text section map

## 1. Introduction

### Purpose
Open the problem narrowly and freeze the scope.

### Must achieve
- define the paper as a **numerical-geometric evidence paper**
- state the model class: 3D Helmholtz + Robin boundary coupling
- explain why ordering vs. artifact is the core question
- explicitly say this paper stops **before** later application work

### Claims carried
- C0
- C10

### Keep out
- long project history
- detector / double-slit discussion
- broad metaphysical or physical interpretation

### Suggested length
~1–1.5 pages

---

## 2. Model and numerical setup

### Purpose
Give only the minimum method needed for an arXiv reader to understand the results.

### Must include
- PDE / eigenproblem statement
- Robin boundary condition notation
- brief geometry families used in the paper
- brief description of discretization / numerical pipeline
- short definition of the main observables used in the main text

### Keep compact
Anything detailed enough for full reproduction should point to the repo.

### Claims carried
- supports C1–C5 and C9 indirectly

### Repo companions
- `README.md`
- `scripts/reproduce_public_minimum.sh`
- `docs/repo_map.md`

### Suggested length
~1.5–2 pages

---

## 3. Ordered structure already appears in the cubic kernel

### Purpose
Establish the first load-bearing result: the cube is not just a noisy numerical spectrum.

### Core message
The base 3D Helmholtz-Robin kernel already shows a nontrivial ordered mode landscape.

### Claims carried
- C1

### Main figure
- **Figure 1**: cube spectrum / beta scan / ordering onset

### Suggested evidence
- lowest-mode ordering
- compact beta comparison
- convergence or robustness note if needed in caption/text

### Repo support
- `scripts/01_cube_solver_basics.py`

### Suggested length
~1.5 pages

---

## 4. Emergent center reading and Helmholtz-Robin role separation

### Purpose
Show both the center-emergence result and the conceptual split between volumetric structure and spectral selection.

### Core messages
- center reading emerges from the computed solutions
- the Robin term modulates spectral position rather than creating the volumetric structure from scratch

### Claims carried
- C2
- C3

### Main figure
- **Figure 2**: center-emergence comparison (cube / quader / hexa) with a compact beta-role reading

### Text emphasis
Use cautious wording such as:
- “supports the interpretation that ...”
- “is consistent with a separation between ...”

### Repo support
- `scripts/02_geometry_center_readout.py`
- `scripts/05_center_emergence.py`

### Suggested length
~2 pages

---

## 5. Geometric hardening: from cube to contrast geometries and deformations

### Purpose
Show the structure survives the move away from the perfect cube and changes in a geometrically meaningful way.

### Core message
Cube-like / hexahedral geometries support mixed volumetric ordering more robustly than anisotropic contrast geometries.

### Claims carried
- C4

### Main figure
- **Figure 3**: tetrahedron / quader / hexahedral contrast and shear hardening

### Keep visible
- not all deformations are equivalent
- anisotropy changes or weakens the ordering family
- the result is selective, not “everything stays the same”

### Repo support
- `scripts/03_tetrahedron_control.py`
- `scripts/04_hexahedron_shear_study.py`
- `scripts/05_center_emergence.py`

### Suggested length
~2 pages

---

## 6. Extension to lattice-linked readout

### Purpose
Show the structure is not confined to an isolated cell.

### Core message
The ordering extends into lattice-linked supercell diagnostics and becomes readable in finite-q style family observables.

### Claims carried
- C5

### Main figure
- **Figure 4**: primitive cubic vs. fcc family readout

### Wording guardrail
This section should *not* expand into a full crystallographic or scattering paper.
It only needs to show that the ordering acquires a reproducible family-level readout beyond the single cell.

### Repo support
- `scripts/06_primitive_cubic_supercell.py`
- `scripts/07_fcc_supercell_readout.py`
- signature/reference CSVs in `results/reference/`

### Suggested length
~1.5–2 pages

---

## 7. Diagnostic hardening under null models and blind tests

### Purpose
Close the paper with the strongest anti-artifact argument.

### Core message
The strongest signals survive selectively under label-shuffle, field-shuffle, and feature-ablation tests, while weaker channels collapse.

### Claims carried
- C9

### Main table / figure
- **Table 1**: null-model / blind-test digest
- optional compact **Figure 5** only if really needed

### What must be visible
- selective survival, not blanket survival
- some channels fail, which strengthens the result
- the paper rules out obvious tested artifact classes, not every conceivable artifact class

### Repo support
- `scripts/main_strand/nullmodel_blindtests/...`
- all reference CSVs under `results/main_strand/nullmodel_blindtests/reference/`

### Suggested length
~1.5–2 pages

---

## 8. Discussion

### Purpose
Condense the meaning of the results without overreaching.

### Core messages
- the paper supports a geometric ordering structure in the tested Helmholtz-Robin model class
- the evidence favors a separation in which Helmholtz carries volumetric structure and Robin modulates selection
- the ordering is volumetric-primary and survives anti-artifact hardening selectively
- the paper is a **backbone paper**, not the full downstream application paper

### Claims carried
- C0
- C3
- C10
- optional short appearance of C8

### One sentence that can anchor the discussion
> The main contribution of Paper 1 is not a complete physical theory but a reproducible numerical-geometric evidence framework for an artifact-hardened ordering structure in the 3D Helmholtz-Robin problem.

### Suggested length
~1–1.5 pages

---

## 9. Conclusion

### Purpose
End cleanly and narrowly.

### Should do
- restate the main evidence result
- restate the anti-artifact result
- state that the paper forms the backbone for later work

### Should not do
- launch into the detector / double-slit story in detail
- claim finality beyond the demonstrated evidence

### Suggested length
~0.5 page

---

## Appendix map

Appendices should stay compact and serve the main text, not replace it.

### Appendix A — Numerical details and implementation notes
Include only what is needed for arXiv completeness:
- mesh/discretization summary
- parameter conventions
- small reproducibility notes

### Appendix B — Compact secondary diagnostics
Use for compressed support of:
- C6 face diagnostics
- C7 field-information diagnostics
- C8 explicit ordering principle / phase-map summary

### Appendix C — Additional tables
Only if needed for null-model or supercell digests.

---

## Main figure / table order

### Figure 1
**Cube kernel ordering**

Supports:
- C1
- partly C3

### Figure 2
**Center emergence and role separation**

Supports:
- C2
- C3

### Figure 3
**Geometry contrast and hexahedral hardening**

Supports:
- C4

### Figure 4
**Lattice-linked family readout: primitive cubic vs. fcc**

Supports:
- C5

### Table 1
**Null-model and blind-test digest**

Supports:
- C9

### Supplement Figure S1
**Face-diagnosis summary**

Supports:
- C6

### Supplement Figure S2
**Field-information summary**

Supports:
- C7

### Supplement Figure S3
**Ordering-principle phase map / continuation digest**

Supports:
- C8

---

## What stays repo-only

The following should be cited as available in the repository but should not be expanded in the main text unless a specific referee-style need later appears:

- full research-history material
- archival provenance layers
- full side-branch breadth
- all intermediate CSVs and auxiliary tables
- reconstructed-vs-original details beyond one short note
- extended family-signature detail beyond what is needed for Figure 4

---

## arXiv formatting guardrails

### Tone
Use restrained wording:
- “we report”
- “we find”
- “supports the interpretation that”
- “consistent with”
- “unlikely to be reducible to the tested artifact classes”

Avoid:
- “proves”
- “establishes final physics”
- “explains reality”
- “demonstrates universal law”

### Citations
The manuscript should be self-contained and conventionally referenced.
The repo is a reproducibility backbone, not a substitute for normal scholarly citation.

### Repo references in text
Keep them sparse and functional:
- one repo statement in Methods
- one reproducibility statement near the end
- possibly one data/code availability statement before references

### Visual economy
Prefer:
- fewer stronger figures
- fewer but clearer captions
- one strong null-model table

---

## Recommended writing order

To keep the arXiv draft disciplined, write in this order:

1. Section 3 — cube ordering
2. Section 4 — center emergence + role split
3. Section 5 — geometric hardening
4. Section 7 — null-model closure
5. Section 6 — lattice-linked readout
6. Introduction
7. Discussion
8. Conclusion
9. Methods
10. Appendix
11. Abstract
12. Title

This order follows the load-bearing numerical spine first and delays framing until the evidence wording is already stable.

---

## Immediate next documents after this section map

The next three paper-facing documents should be:

1. `docs/paper1_figure_caption_drafts.md`
2. `docs/paper1_methods_map.md`
3. `docs/paper1_main_text_outline.md`

That sequence will take the project from repository coordination into actual manuscript drafting.
