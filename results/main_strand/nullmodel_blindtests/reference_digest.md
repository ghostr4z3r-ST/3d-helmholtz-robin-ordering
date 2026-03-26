# Nullmodell-/Blindtest-Referenzdigest

Diese Datei fasst **mitgelieferte Referenz-CSV** zusammen. Sie reproduziert die historischen Tests nicht neu, sondern bietet einen kuratierten Leseeinstieg in Phase J des Hauptstrangs.

## Feature-Ablation (reduzierter erster Lauf)
- cube_raw_1F: voll=0.643, wichtigste Gruppe=abs_coeffs, beste Einzelgruppe=abs_coeffs (0.619)
- cube_raw_6F: voll=0.714, wichtigste Gruppe=face_means, beste Einzelgruppe=face_means (0.762)
- cube_scan_1F: voll=0.484, wichtigste Gruppe=patch_std, beste Einzelgruppe=patch_mean (0.492)
- cube_volume: voll=0.952, wichtigste Gruppe=abs_pair3d, beste Einzelgruppe=abs_axis (1.000)
- supercell_raw_1F: voll=0.417, wichtigste Gruppe=abs_coeffs, beste Einzelgruppe=abs_coeffs (0.639)
- supercell_scan_1F: voll=0.361, wichtigste Gruppe=patch_std, beste Einzelgruppe=patch_std (0.486)
- supercell_volume: voll=0.917, wichtigste Gruppe=uniformity_tail, beste Einzelgruppe=uniformity_tail (1.000)

## Feld-Shuffle auf Würfel-Aufgaben
- cube_raw_1F: echt=0.744, Shuffle=0.346 ± 0.051, p=0.020
- cube_raw_6F: echt=0.786, Shuffle=0.381 ± 0.127, p=0.020
- cube_volume: echt=1.000, Shuffle=0.344 ± 0.115, p=0.020

## Feld-Shuffle auf Superzelle / q-Diagnostik
- supercell_raw_1F: echt=0.593, Shuffle=0.491 ± 0.076, p=0.048
- supercell_raw_6F: echt=0.500, Shuffle=0.481 ± 0.110, p=0.524
- supercell_volume: echt=0.889, Shuffle=0.958 ± 0.044, p=1.000

## q-Metriken mit positiver Shuffle-Robustheit
- primitive_q:best_q_amp: echt=0.592, Shuffle=0.423 ± 0.020, z=8.422, p=0.048
- fcc:fcc_face_uniformity: echt=0.378, Shuffle=0.326 ± 0.031, z=1.708, p=0.095
- fcc:fcc_face_qamp: echt=0.601, Shuffle=0.490 ± 0.017, z=6.734, p=0.048
