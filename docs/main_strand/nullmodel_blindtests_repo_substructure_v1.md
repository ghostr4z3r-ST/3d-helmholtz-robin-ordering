# Nullmodell-/Blindtest-Block – Repo-Unterstruktur v1

## Ziel
Diese Unterstruktur ordnet den Nullmodell-/Blindtest-Block als **letzte methodische Härtungsphase des Paper-1-Hauptstrangs** ein.

Sie soll sichtbar machen:
- dass der Block **zum Hauptstrang** gehört,
- welche Referenzdaten bereits vorliegen,
- welche Teile erst historisch belegt sind,
- und wie man diesen Abschnitt im öffentlichen Repo schlank, aber ehrlich abbildet.

## Empfohlene Struktur

```text
paper1-github-repo/
├── scripts/
│   └── main_strand/
│       └── nullmodel_blindtests/
│           ├── 01_reference_digest.py
│           └── reproduce_nullmodel_digest.sh
├── results/
│   └── main_strand/
│       └── nullmodel_blindtests/
│           ├── README.md
│           ├── reference/
│           │   ├── feature_ablation_tests_reduced.csv
│           │   ├── feature_ablation_summary_reduced.csv
│           │   ├── field_shuffle_null_tests_50.csv
│           │   └── field_shuffle_supercell_null_tests.csv
│           ├── reference_digest.md
│           └── historical_pending/
│               └── README.md
├── docs/
│   └── main_strand/
│       ├── nullmodel_blindtests_forschungshistorie_v1.md
│       ├── nullmodel_blindtests_repo_substructure_v1.md
│       └── nullmodel_blindtests_status.md
└── archive/
    └── main_strand/
        └── nullmodel_blindtests/
            ├── chatexports/
            │   ├── chatexport_b_1_...
            │   ├── ...
            │   └── chatexport_b_6_...
            ├── original_scripts/
            │   └── field_shuffle_null_tests.py
            └── history_only_missing_artifacts.md
```

## Rollen der Ebenen

### 1. `scripts/main_strand/nullmodel_blindtests/`
Hier liegt **kein vorgetäuschter Re-Run** der historischen Tests, sondern ein schlanker Einstieg in die vorhandenen Referenzdaten.

Empfohlen ist zunächst nur:
- `01_reference_digest.py`

Diese Ebene dient dazu, die Phase J für Außenstehende lesbar zu machen, ohne mehr Reproduzierbarkeit zu behaupten, als aktuell vorhanden ist.

### 2. `results/main_strand/nullmodel_blindtests/reference/`
Hier liegen die **echten mitgelieferten Referenz-CSV**. Diese Dateien tragen bereits einen substanziellen Teil der historischen Härtungsphase.

### 3. `results/main_strand/nullmodel_blindtests/reference_digest.md`
Diese Datei ist eine **kuratierte Lesebrücke** über die Referenz-CSV. Sie ist nützlich für GitHub, weil Leser sofort sehen, was die Kernbefunde dieses Blocks sind.

### 4. `results/main_strand/nullmodel_blindtests/historical_pending/`
Hier sollten bewusst **keine erfundenen Ersatzdateien** liegen. Stattdessen dokumentiert ein README, welche historischen Artefakte noch fehlen.

### 5. `docs/main_strand/`
Hier liegt die **eigentliche Ordnungsdokumentation** dieser Phase:
- Forschungshistorie
- Substruktur-Plan
- Status-/Provenienznotiz

### 6. `archive/main_strand/nullmodel_blindtests/`
Hier bleibt die Rohprovenienz:
- geordnete Chat-Historie
- wiedergefundene Originalskripte
- Manifest fehlender historischer Artefakte

## Was aktuell als restauriert gelten sollte

### Restaurierter Kern
- Referenzdaten: `feature_ablation_*`, `field_shuffle_null_tests_50.csv`, `field_shuffle_supercell_null_tests.csv`
- Provenienz: `chatexport_b_1` bis `chatexport_b_6`
- historisches Original: `field_shuffle_null_tests.py`
- kuratierter GitHub-Einstieg: `01_reference_digest.py`

### Historisch belegt, aber noch nicht dateiseitig restauriert
- Label-Shuffle-Skript + Ergebnisdatei
- Ablationsskript
- Supercell-Feld-Shuffle-Skript
- q-dominierter Nachtest als Dateien
- historische Zusatzabhängigkeit `face_information_sidebranch_extended.py`

## Empfohlene Benennung im öffentlichen Repo
Ich würde diesen Block nach außen nicht „nullmodelle“ allein nennen, sondern zweistufig:

- technisch: `nullmodel_blindtests`
- inhaltlich: **Phase J – methodische Härtung**

So bleibt sichtbar, dass dieser Teil **nicht optionaler Nebenstoff**, sondern eine tragende methodische Schicht des Hauptstrangs ist.

## Reihenfolge für jemanden, der diesen Block nachvollziehen will

1. `docs/main_strand/nullmodel_blindtests_forschungshistorie_v1.md`
2. `results/main_strand/nullmodel_blindtests/reference_digest.md`
3. die vier Referenz-CSV im Unterordner `reference/`
4. `docs/main_strand/nullmodel_blindtests_status.md`
5. erst danach die Archivquellen

## Praktische Konsequenz
Diese Struktur erreicht genau die gewünschte Trennung:
- Das Paper bleibt schlank.
- Der Hauptstrang zeigt eine methodische Härtungsebene.
- Die historischen Dateilücken bleiben offen sichtbar.
- GitHub wirkt trotzdem geordnet und ernsthaft.
