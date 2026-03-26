# Reproducibility status

## Directly runnable from this public repo

- cube solver, convergence scan, and basic β scan,
- geometry-center readout for cube / quader / deformed hexahedron,
- tetrahedron control,
- first hexahedral shear study,
- primitive-cubic supercell readout,
- fcc supercell readout,
- face-diagnosis branch,
- field-information branch,
- geometry continuation and phase-map branch,
- null-model reference digest.

## Historical wrappers included

Several later steps are run through wrappers around recovered historical scripts.
This keeps the public paths clean while preserving the original numerical logic.

## Reconstructed or partial-provenance items

Two files remain provenance-qualified rather than fully original-source-complete:

- `field_shuffle_supercell_null_tests.py` survived historically only as a stub, so the repo ships an explicitly marked reconstruction.
- `supercell_qdominated_field_shuffle.py` is preserved only as a stub; the scientifically important artifact here is the versioned reference CSV.

## Practical reading for the paper

For Paper 1, this is now strong enough to serve as the reproducibility backbone.
The remaining caveat is provenance transparency, not missing numerical substance.
