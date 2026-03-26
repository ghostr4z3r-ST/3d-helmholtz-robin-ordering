# Nullmodell-/Blindtest-Block – Forschungshistorie v1

## Zweck
Diese Notiz rekonstruiert den **letzten methodischen Härtungsblock des Paper-1-Hauptstrangs**. Sie dient als Provenienz- und Strukturhilfe für das Repository und hält fest, wie aus positiven Befunden ein expliziter Selbsttest-Block wurde.

## Quellenbasis
Verwendet wurden die geordneten Chat-Exporte `chatexport_b_1` bis `chatexport_b_6` sowie der aktuell wiedergefundene Daten-/Skriptstand (`feature_ablation_*`, `field_shuffle_*`, `field_shuffle_null_tests.py`).

## Historischer Kernbefund in einem Satz
Der Block verschiebt den Hauptstrang von „wir finden Ordnung“ zu der stärkeren Aussage:

> Die belastbaren Befunde überstehen gezielte Nullmodelle und Selbsttests **nicht pauschal**, sondern **selektiv**: robuste Volumen- und q-Ordnungsdiagnostik bleibt stehen, während schwächere oder zu grob gemischte Readouts offen einbrechen.

---

## Phase B0 – Der Schritt zur harten Selbstwiderlegung
Nachdem Modenhierarchie, Zentrumsemergenz, hexaedrische Robustheit und Superzellenordnung bereits positiv entwickelt waren, wird die methodische Frage schärfer:

> Würde dieselbe Diagnostik auch dann noch „Ordnung“ finden, wenn in Wahrheit keine spezifische Ordnungsstruktur vorliegt?

Genau daraus entsteht der Nullmodell-/Blindtest-Block. Historisch wichtig ist, dass dieser Schritt **nicht** als Zusatzspielerei, sondern als bewusste Korrektur früherer theoretischer Schwächen begriffen wird: Was fehlt, sind jetzt harte Selbstwiderlegungstests.

## Phase B1 – Label-Shuffle als einfachster erster Falsifikationstest
Der Block beginnt mit dem methodisch einfachsten und zugleich stärksten Einstieg: **bestehende Features bleiben gleich, nur die Klassenlabels werden permutiert**.

### Getestete Aufgaben
- Würfel: rohe Einzelflächen-, Sechsflächen-, Scan- und Volumenlesung
- Superzelle primitive vs. fcc: rohe Einzelflächen-, Sechsflächen-, Scan- und Volumenlesung

### Historische Lesung
Der erste Label-Shuffle-Test bringt sofort die zentrale Differenzierung:
- **klar robust**: `cube_volume`, `cube_raw_6F`, `cube_raw_1F`, `cube_scan_1F`, `supercell_volume`, `supercell_raw_1F`, `supercell_scan_1F`
- **nicht überzeugend**: `cube_scan_6F`, `supercell_raw_6F`, `supercell_scan_6F`

Damit wird ein früher Verdacht ausdrücklich zurückgewiesen: Die starken Resultate entstehen **nicht bloß aus der Pipelineform**. Zugleich wird sichtbar, dass **nicht jeder Teil der Diagnostik gleich gut trägt**. Schon hier verdichtet sich die spätere Lesung:
- Volumenordnung ist primär,
- Einzelflächen tragen partiell,
- mehr Randflächen zusammen sind **nicht automatisch** besser.

## Phase B2 – Feature-Ablation statt diffuser Gesamtbestätigung
Nach dem Label-Shuffle folgt der zweite Härtungsschritt: **Welche Featuregruppen tragen die guten Treffer wirklich?**

### Methodik
Für repräsentative Würfel- und Superzellenaufgaben werden verglichen:
- volle Featuremenge,
- nur eine Featuregruppe,
- alles außer einer Featuregruppe.

### Historische Konsequenz
Der Befund wird inhaltlich schärfer:
- Würfel-Volumen lebt **nicht** an diffuser Information, sondern stark an strukturellen Gruppen wie `abs_axis` und `abs_pair3d`.
- Supercell-Volumen lebt in der ersten Fassung besonders an `uniformity_tail`, also einer Ordnungs-/Uniformitätsstruktur.
- Auf Flächen tragen vor allem **absolute Koeffizientenstrukturen**, nicht feine Vorzeichen- oder Korrelationsdetails.
- Patch-Scan-Ideen bleiben eher schwach und unreif.

Dieser Schritt ist wissenschaftlich entscheidend, weil er zeigt: Die starken Befunde sitzen auf **identifizierbaren strukturellen Merkmalen** und nicht auf beliebiger Merkmalssammlung.

## Phase B3 – Feld-Shuffle auf dem Würfel: Zerstörung räumlicher Ordnung
Danach wird der härtere Test gefahren: **nicht mehr Labels oder Features vertauschen, sondern die räumliche Feldstruktur selbst zerstören**, bei erhaltener Norm und grober Amplitudenverteilung.

### Erster Lauf
Getestet werden die stärksten Würfel-Aufgaben:
- `cube_raw_1F`
- `cube_raw_6F`
- `cube_volume`

### Historische Lesung
Alle drei Klassifikationsleistungen brechen unter Feld-Shuffle deutlich ein. Damit wird ein zweiter wichtiger Verdacht adressiert:

> Die Diagnostik liest nicht bloß Wertestatistik, sondern echte räumliche Ordnung.

Der Würfelblock wird dadurch methodisch deutlich glaubwürdiger.

## Phase B4 – Übertragung auf die Gitterebene: Supercell-Feld-Shuffle
Der nächste Schritt überträgt denselben Härtecheck auf die **3×3×3-Superzelle** für primitive finite-`q` vs. fcc.

### Getestet werden
- `supercell_raw_1F`
- `supercell_raw_6F`
- `supercell_volume`
- zusätzlich q-/Ordnungsmetriken wie `best_q_amp` und `fcc_face_qamp`

### Historische Konsequenz
Gerade dieser Test ist wichtig, weil er **nicht einfach nur bestätigt**, sondern trennt:
- `supercell_raw_1F` bleibt knapp über Shuffle und trägt damit noch echte räumliche Information.
- `supercell_raw_6F` ist **nicht überzeugend**.
- die bisherige gemischte `supercell_volume`-Klassifikation ist in dieser Form **nicht shuffle-robust**.
- die eigentlichen **q-Metriken** bleiben dagegen klar robust.

Damit verschiebt sich die Hauptlesung auf Gitterebene entscheidend:

> Der harte Kern der Supercell-Befunde liegt eher in **q-/Sublattice-/Peak-Signaturen** als in jeder beliebigen globalen Volumenklassifikation.

## Phase B5 – q-dominierter Nachschnitt der Supercell-Volumenfeatures
Aus dem Supercell-Shuffle folgt der nächste logische Schritt: **Volumenfeatures q-/ordnungsdominiert neu zuschneiden** und erneut testen.

### Historische Lesung
Auch dieser Nachtest bleibt ehrlich zweigeteilt:
- q-Metriken wie `center_qamp_mean` und `face_qamp_mean` bleiben klar feld-shuffle-robust,
- eine einzelne kompakte Gesamtklassifikation `supercell_volume_qdominated` ist aber immer noch **nicht überzeugend robust**.

Das ist kein Rückschritt, sondern eine weitere Präzisierung:
- Würfel: starke Volumenordnung
- Superzelle: starke q-Ordnung
- primitive vs. fcc: Unterschied sitzt derzeit stärker in der **Art der Ordnungslesung** als in einem einzigen globalen Kompaktklassifikator

## Phase B6 – Endformulierung des Härtungsblocks
Am Ende steht keine pauschale Siegererzählung, sondern eine geschärfte Bilanz:

> Die starken Befunde sind **nicht beliebig**, **nicht rein labelgetrieben**, **nicht bloß amplitudenstatistisch** und **nicht auf diffuse Features verteilt**.

Genauso wichtig ist aber der zweite Halbsatz:

> Nicht jede zusammengesetzte Rand- oder Volumenmischung trägt dieselbe Ordnung gleichermaßen robust.

Gerade diese kombinierte Aussage macht die Phase methodisch wertvoll.

---

## Was für das Repository daraus folgt
Dieser Block gehört im Repo **in den Hauptstrang**, nicht in einen Nebenast. Die richtige Lesart ist:

1. **Härtungsschicht des Hauptstrangs**
   - Label-Shuffle
   - Feature-Ablation
   - Würfel-Feld-Shuffle
   - Supercell-Feld-Shuffle
   - q-dominierter Nachtest

2. **Direkt vorhandene Referenzartefakte**
   - `feature_ablation_tests_reduced.csv`
   - `feature_ablation_summary_reduced.csv`
   - `field_shuffle_null_tests_50.csv`
   - `field_shuffle_supercell_null_tests.csv`

3. **Historisch belegt, noch nicht dateiseitig restauriert**
   - Label-Shuffle-Skript und Ergebnisdatei
   - Ablationsskript
   - Supercell-Feld-Shuffle-Skript
   - q-dominierter Nachtest als Dateien

## Offene Wiederherstellungen
Aktuell fehlen als Originaldateien noch:
- `label_shuffle_permutation_tests.py`
- `label_shuffle_permutation_tests.csv`
- `feature_ablation_tests_reduced.py`
- `field_shuffle_supercell_null_tests.py`
- `supercell_qdominated_field_shuffle.py`
- `supercell_qdominated_field_shuffle.csv`
- sowie die Abhängigkeit `face_information_sidebranch_extended.py`, die für das historische `field_shuffle_null_tests.py` relevant ist

## Empfehlung
Für die öffentliche GitHub-Fassung sollte diese Phase so erscheinen:
- als **klar dokumentierte Hauptstrang-Härtung**,
- mit echten Referenzdaten,
- mit offen markierten Dateilücken,
- und mit der klaren methodischen Pointe: **Der Block schärft, statt bloß zu bestätigen.**
