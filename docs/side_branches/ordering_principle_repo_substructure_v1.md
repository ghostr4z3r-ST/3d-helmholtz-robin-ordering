# Ordnungsprinzip-Nebenast – Repo-Unterstruktur v1

## Ziel
Diese Unterstruktur ordnet den Nebenast zum **expliziten Ordnungsprinzip** so ein, dass er
- als **Nebenast von Paper 1** erkennbar bleibt,
- zunächst **dokumentations- und provenienzstark** erscheint,
- und später bei Wiederfund oder Rekonstruktion der fehlenden Artefakte in einen reproduzierbaren Diagnosestrang überführt werden kann.

## Empfohlene Struktur

```text
paper1-github-repo/
├── docs/
│   └── side_branches/
│       ├── ordering_principle_forschungshistorie_v1.md
│       ├── ordering_principle_repo_substructure_v1.md
│       └── ordering_principle_status.md
├── results/
│   └── side_branches/
│       └── ordering_principle/
│           ├── README.md
│           └── historical_pending/
│               └── README.md
└── archive/
    └── side_branches/
        └── ordering_principle/
            ├── notes/
            │   ├── Prinzipien-Notiz.txt
            │   └── Ordnungsphasennotiz.txt
            ├── chatexports/
            │   ├── chatexport_21_...
            │   ├── chatexport_22_...
            │   ├── chatexport_23_...
            │   └── chatexport_24_...
            └── history_only_missing_artifacts.md
```

## Rollen der Ebenen

### 1. `docs/side_branches/ordering_principle_*.md`
Hier liegt die **lesbare Rekonstruktion** des Nebenasts:
- Forschungshistorie
- Substruktur-Plan
- Status-/Provenienznotiz

Das ist die wichtigste Ebene, solange die numerischen Originalartefakte noch fehlen.

### 2. `results/side_branches/ordering_principle/`
Hier sollten aktuell **keine Fake-CSV und keine künstlich erzeugten PNG-Dateien** liegen.

Stattdessen genügt ein README, das dokumentiert, welche historisch genannten Ergebnisse später hier hineingehören würden:
- Pilot-CSV der ersten Kontinuation
- volle / kompakte CSV der Phasenkarte
- Übersichtsbild der Ordnungsphasenkarte

### 3. `archive/side_branches/ordering_principle/notes/`
Hier liegen die **begrifflichen Verdichtungen**:
- `Prinzipien-Notiz.txt`
- `Ordnungsphasennotiz.txt`

Diese Dateien sind keine Rohdaten, aber sie halten den argumentativen Endstand des Nebenasts fest.

### 4. `archive/side_branches/ordering_principle/chatexports/`
Hier bleibt die **Rohprovenienz**:
- neuer Arbeitsblock,
- erster Pilot,
- Verschärfung zur Phase-Map,
- kompakte Ordnungsphasenkarte.

### 5. Spätere Ausbaustufe
Erst wenn die historischen Skripte/CSV/PNG wiedergefunden oder sauber kuratiert rekonstruiert sind, sollte zusätzlich eine produktive Ebene entstehen:

```text
src/paper1/ordering_principle.py
scripts/side_branches/ordering_principle/01_geometry_continuation_pilot.py
scripts/side_branches/ordering_principle/02_geometry_phase_map.py
results/side_branches/ordering_principle/reference/...
```

Bis dahin sollte das Repo hier bewusst **dokumentieren statt simulieren**.

## Was aktuell als restauriert gelten sollte

### Restaurierter Kern
Noch **kein** lauffähiger numerischer Kern.

Restauriert und sauber archiviert sind derzeit nur:
- die geordneten Chat-Exporte,
- die Prinzipien-Notiz,
- die Ordnungsphasen-Notiz,
- die daraus gebaute Forschungshistorie und Statusstruktur.

### Historisch belegt, aber noch nicht dateiseitig restauriert
- `geometry_continuation_study.py`
- `geometry_continuation_study_first_test.csv`
- `geometry_phase_map.py`
- `geometry_phase_map_full.csv`
- `geometry_phase_map_summary.csv`
- `geometry_phase_map_compact.csv`
- `geometry_phase_map_overview.png`

## Empfohlene Benennung im öffentlichen Repo
Ich würde den Nebenast zweistufig benennen:
- technisch nach außen: `ordering_principle`
- historisch/inhaltlich: `explizites Ordnungsprinzip`

Das hält die Repo-Namen knapp, ohne den inhaltlichen Charakter zu verlieren.

## Reihenfolge für jemanden, der den Nebenast nachvollziehen will

1. `docs/side_branches/ordering_principle_forschungshistorie_v1.md`
2. `docs/side_branches/ordering_principle_status.md`
3. `archive/side_branches/ordering_principle/notes/Prinzipien-Notiz.txt`
4. `archive/side_branches/ordering_principle/notes/Ordnungsphasennotiz.txt`
5. erst danach die einzelnen Chat-Exporte

## Praktische Konsequenz
Diese Struktur erlaubt genau das, was Paper 1 hier braucht:
- Das Paper bleibt schlank.
- Das Repo verliert den Syntheseast nicht.
- Fehlende numerische Zwischenartefakte werden offen markiert.
- Die Wiederherstellung kann später gezielt nachgerüstet werden, ohne die Struktur erneut umzubauen.
