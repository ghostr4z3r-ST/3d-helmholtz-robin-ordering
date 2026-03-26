# Provenance and reconstructions

This repository does **not** treat all files as equal.

## Curated code

Files in `src/paper1/` and the main wrapper scripts are the clean public layer.
They are the first place a reader should look.

## Recovered historical code

Files under `archive/.../original_scripts/` are retained because they capture how later parts of the project were actually computed.
Where needed, wrapper scripts run them in a repo-local staging environment so that obsolete absolute paths do not leak into the public release.

## Reconstructed files

Two items are explicitly provenance-qualified:

1. `field_shuffle_supercell_null_tests_reconstructed.py`  
   Included because the historical source survived only as a stub, while the result CSV survived.

2. `supercell_qdominated_field_shuffle.py`  
   The script survives only as a stub. The versioned CSV is the meaningful scientific artifact.

## Omitted from the public repo

Bulk chat-export logs used during reconstruction are kept outside the public GitHub repository.
They were useful for rebuilding the research history, but they do not improve the public scientific readability of the repo itself.
