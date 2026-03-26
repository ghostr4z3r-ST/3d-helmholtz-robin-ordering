# Feldinformationsdiagnose-Nebenast – Repo-Unterstruktur v1

## Ziel
Diese Unterstruktur ordnet den Feldinformationsdiagnose-Nebenast so ein, dass er
- als **Nebenast von Paper 1** erkennbar bleibt,
- einen **restaurierten reproduzierbaren Kern** hat,
- und gleichzeitig die **historische Vollständigkeit** über Chat-Backups und Artefakt-Manifeste bewahrt.

## Empfohlene Struktur

```text
paper1-github-repo/
├── src/paper1/
│   └── field_information.py
├── scripts/
│   └── side_branches/
│       └── field_information/
│           ├── 01_field_information_first_test.py
│           └── reproduce_field_information_v1.sh
├── results/
│   └── side_branches/
│       └── field_information/
│           ├── README.md
│           ├── reference/
│           │   └── field_information_distribution_first_test.csv
│           └── historical_pending/
│               └── README.md
├── docs/
│   └── side_branches/
│       ├── field_information_forschungshistorie_v1.md
│       ├── field_information_repo_substructure_v1.md
│       └── field_information_status.md
└── archive/
    └── side_branches/
        └── field_information/
            ├── original_scripts/
            │   └── field_information_distribution.py
            ├── chatexports/
            │   ├── chatexport_field_1_...
            │   ├── chatexport_field_2_...
            │   ├── chatexport_field_3_...
            │   └── chatexport_field_projektkarte.txt
            └── history_only_missing_artifacts.md
```

## Rollen der Ebenen

### 1. `src/paper1/field_information.py`
Das ist die **kuratierte, paketfähige Kernfassung** des restaurierten ersten Feldinformations-Tests.

### 2. `scripts/side_branches/field_information/`
Hier liegen nur die **schlanken Einstiegsskripte** für den Nebenast.

Empfohlen ist genau ein erster reproduzierbarer Einstieg:
- `01_field_information_first_test.py`

Weitere Skripte für Superzellen-Distribution sollten erst dann hier ergänzt werden, wenn die historischen Originale oder eine belastbare kuratierte Rekonstruktion vorliegen.

### 3. `results/side_branches/field_information/reference/`
Hier liegen die **historischen Referenz-CSV**, auf die sich die Forschungshistorie bezieht.

Aktuell gehört hier hinein:
- `field_information_distribution_first_test.csv`

### 4. `results/side_branches/field_information/historical_pending/`
Hier liegen **noch keine Fake-Dateien**. Stattdessen dokumentiert ein README, welche Ergebnisse historisch belegt sind, aber als Dateien noch fehlen.

Das schützt das Repo davor, mehr Reproduzierbarkeit zu behaupten, als tatsächlich vorhanden ist.

### 5. `docs/side_branches/`
Hier liegen die **lesbaren Ordnungsdokumente**:
- Forschungshistorie
- Substruktur-Plan
- Status-/Provenienznotiz

### 6. `archive/side_branches/field_information/`
Hier bleibt die **Rohprovenienz**:
- exportierte Chat-Historie,
- wiedergefundene Originalskripte,
- Manifest fehlender historischer Artefakte.

## Was aktuell als restauriert gelten sollte

### Restaurierter Kern
- kuratierter Code: `src/paper1/field_information.py`
- historisches Original: `archive/side_branches/field_information/original_scripts/field_information_distribution.py`
- Referenzdaten: `results/side_branches/field_information/reference/field_information_distribution_first_test.csv`

### Historisch belegt, aber noch nicht dateiseitig restauriert
- Superzellen-Distribution primitive vs. fcc
- zugehörige Ergebnis-CSV der 3×3×3-Auswertung

## Empfohlene Benennung im öffentlichen Repo
Ich würde den Nebenast im öffentlichen GitHub-Repo zweistufig führen:

- technisch: `field_information`
- in der Forschungshistorie: `Feldinformationsdiagnose`

Dadurch bleibt der Bezug zum wiedergefundenen Originalskript erhalten, während die Repo-Struktur nach außen knapp und verständlich bleibt.

## Reihenfolge für jemanden, der den Nebenast nachvollziehen will

1. `docs/side_branches/field_information_forschungshistorie_v1.md`
2. `scripts/side_branches/field_information/01_field_information_first_test.py`
3. `results/side_branches/field_information/reference/field_information_distribution_first_test.csv`
4. `docs/side_branches/field_information_status.md`
5. erst danach die Archivquellen

## Praktische Konsequenz
Diese Struktur erlaubt genau das, was du für Paper 1 brauchst:
- Das Paper bleibt schlank.
- Das Repo zeigt den restaurierten Kern klar.
- Die Forschungshistorie bleibt vollständig.
- Fehlende historische Dateien werden offen markiert statt versteckt.
