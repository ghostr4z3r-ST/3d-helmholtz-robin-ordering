# Flächendiagnose-Nebenast – Repo-Unterstruktur v1

## Ziel
Diese Unterstruktur ordnet den Flächendiagnose-Nebenast so ein, dass er
- als **Nebenast von Paper 1** erkennbar bleibt,
- einen **restaurierten reproduzierbaren Kern** hat,
- und gleichzeitig die **historische Vollständigkeit** über Chat-Backups und Artefakt-Manifeste bewahrt.

## Empfohlene Struktur

```text
paper1-github-repo/
├── src/paper1/
│   └── face_information.py
├── scripts/
│   └── side_branches/
│       └── face_diagnosis/
│           ├── 01_face_information_first_test.py
│           └── reproduce_face_diagnosis_v1.sh
├── results/
│   └── side_branches/
│       └── face_diagnosis/
│           ├── README.md
│           ├── reference/
│           │   └── face_information_first_test.csv
│           └── historical_pending/
│               └── README.md
├── docs/
│   └── side_branches/
│       ├── face_diagnosis_forschungshistorie_v1.md
│       ├── face_diagnosis_repo_substructure_v1.md
│       └── face_diagnosis_status.md
└── archive/
    └── side_branches/
        └── face_diagnosis/
            ├── original_scripts/
            │   └── face_information_sidebranch.py
            ├── chatexports/
            │   ├── chatexport_f_1_...
            │   ├── ...
            │   └── chatexport_f_9_...
            └── history_only_missing_artifacts.md
```

## Rollen der Ebenen

### 1. `src/paper1/face_information.py`
Das ist die **kuratierte, paketfähige Kernfassung** des restaurierten ersten Face-Information-Tests.

### 2. `scripts/side_branches/face_diagnosis/`
Hier liegen nur die **schlanken Einstiegsskripte** für den Nebenast.

Empfohlen ist genau ein erster reproduzierbarer Einstieg:
- `01_face_information_first_test.py`

Weitere Skripte für erweiterte Features oder Patch-Scans sollten erst dann hier ergänzt werden, wenn die historischen Originale oder eine belastbare kuratierte Rekonstruktion vorliegen.

### 3. `results/side_branches/face_diagnosis/reference/`
Hier liegen die **historischen Referenz-CSV**, auf die sich die Forschungshistorie bezieht.

Aktuell gehört hier hinein:
- `face_information_first_test.csv`

### 4. `results/side_branches/face_diagnosis/historical_pending/`
Hier liegen **noch keine Fake-Dateien**. Stattdessen dokumentiert ein README, welche Ergebnisse historisch belegt sind, aber als Dateien noch fehlen.

Das schützt das Repo davor, mehr Reproduzierbarkeit zu behaupten, als tatsächlich vorhanden ist.

### 5. `docs/side_branches/`
Hier liegen die **lesbaren Ordnungsdokumente**:
- Forschungshistorie
- Substruktur-Plan
- Status-/Provenienznotiz

Diese Ebene ist wichtig, damit Außenstehende den Nebenast einordnen können, ohne die gesamte Chat-Historie lesen zu müssen.

### 6. `archive/side_branches/face_diagnosis/`
Hier bleibt die **Rohprovenienz**:
- exportierte Chat-Historie,
- wiedergefundene Originalskripte,
- Manifest fehlender historischer Artefakte.

Das Archiv ist kein Teil der schlanken Reproduktionsoberfläche, aber entscheidend für die Forschungshistorie.

## Was aktuell als restauriert gelten sollte

### Restaurierter Kern
- kuratierter Code: `src/paper1/face_information.py`
- historisches Original: `archive/side_branches/face_diagnosis/original_scripts/face_information_sidebranch.py`
- Referenzdaten: `results/side_branches/face_diagnosis/reference/face_information_first_test.csv`

### Historisch belegt, aber noch nicht dateiseitig restauriert
- erweiterte Face-Features und Robustheit
- Patch-Scan auf Würfel
- Patch-Scan auf primitive Superzelle / fcc

## Empfohlene Benennung im öffentlichen Repo
Ich würde den Nebenast im öffentlichen GitHub-Repo **nicht nur** „face information“ nennen, sondern zweistufig:

- technisch: `face_diagnosis`
- intern/historisch: `face_information`

Dadurch bleibt der Dateiname der wiedergefundenen Originalskripte erhalten, während die Repo-Struktur nach außen verständlicher wird.

## Reihenfolge für jemanden, der den Nebenast nachvollziehen will

1. `docs/side_branches/face_diagnosis_forschungshistorie_v1.md`
2. `scripts/side_branches/face_diagnosis/01_face_information_first_test.py`
3. `results/side_branches/face_diagnosis/reference/face_information_first_test.csv`
4. `docs/side_branches/face_diagnosis_status.md`
5. erst danach die Archivquellen

## Praktische Konsequenz
Diese Struktur erlaubt genau das, was du brauchst:
- Das Paper bleibt schlank.
- Das Repo zeigt den restaurierten Kern klar.
- Die Forschungshistorie bleibt vollständig.
- Fehlende historische Dateien werden offen markiert statt versteckt.
