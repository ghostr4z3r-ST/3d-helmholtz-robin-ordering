# Reproducibility status

This repository accompanies the manuscript **Reproducible Low-Mode Ordering in a 3D Helmholtz--Robin Eigenproblem on Cube-Like and Lattice-Linked Geometries**.

## Status summary

The repository is intended to reproduce the numerical analyses reported in the manuscript.

### Directly runnable curated core

- low-mode cube solver and center readout
- tetrahedron control
- geometric hardening on deformed hexahedra
- primitive-cubic and fcc lattice-linked readouts

### Supported by wrapper scripts

Some manuscript-supporting analyses still run through thin wrappers around legacy source files preserved in `archive/`.
These wrappers are kept only where they remain useful for reproduction.

### Reconstructed component

One supercell field-shuffle component survives publicly only as a reconstruction because the recovered original script was a stub. This is explicitly marked in the relevant script and documentation.
