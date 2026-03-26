# Ordnungsprinzip-Nebenast – Forschungshistorie v1

## Zweck
Diese Notiz rekonstruiert den **Nebenast zum expliziten Ordnungsprinzip** als geordneten Forschungsverlauf innerhalb von Paper 1. Sie dient nicht als fertiges Manuskript, sondern als interne Provenienz-, Struktur- und Restaurationshilfe für das Repository.

## Quellenbasis
Verwendet wurden:
- `chatexport_21_neuer Arbeitsblock Ordnungsprinzip explizit machen.txt`
- `chatexport_22_erster Testlauf und Pilotbefund.txt`
- `chatexport_23_nächster Schritt.txt`
- `chatexport_24_ordnungsphasenkarten und befunde.txt`
- `Prinzipien-Notiz.txt`
- `Ordnungsphasennotiz.txt`

Wichtig: Für diesen Nebenast liegen aktuell **geordnete Protokolle und Notizen**, aber noch **keine wiedergefundenen Originalskripte oder CSV/PNG-Artefakte** vor.

## Historischer Kernbefund in einem Satz
Die zuvor getrennt stehenden Befunde aus Hauptstrang, Flächendiagnose und Feldinformationsdiagnose werden hier erstmals zu einem expliziten Arbeitsprinzip zusammengezogen:

> **3D-Helmholtz trägt die Moden- und Ordnungslandschaft, Robin selektiert die bevorzugte Ordnungsfamilie, und die geometrische Gleichwertigkeit der drei Raumrichtungen entscheidet wesentlich darüber, ob gemischte volumetrische oder achsengebundene Ordnung bevorzugt wird.**

---

## Phase OP0 – Verdichtung bereits vorhandener Befunde
Dieser Nebenast entsteht nicht aus dem Nichts. Er setzt auf bereits stabilisierte Ergebnisse auf:
- das volumetrische Eigenwertproblem als Strukturträger,
- Robin als Selektions- statt Erzeugungsmechanismus,
- Volumen-vs.-Fläche als robuste Trennlinie,
- primitive kubisch vs. fcc als erste beobachtungsnahe Gegenlesung,
- sowie die Einsicht, dass nicht ein Einzelwert, sondern eine ganze **Ordnungsdiagnostik** aus Peakfamilie, Ordnungsvektor, Sublattice-Polarisation und Volumen-vs.-Fläche-Lesbarkeit physikalisch lesbar wird.

Die spätere Prinzipien-Notiz fasst genau diese Verdichtung ausdrücklich zusammen.

## Phase OP1 – Neuer Arbeitsblock: das wiederkehrende Prinzip explizit machen
Danach wird zum ersten Mal offen formuliert, dass die vielen robusten Befunde auf ein gemeinsames Prinzip verweisen könnten.

### Arbeitsvermutung
Die Ordnung wird durch drei Größen gemeinsam bestimmt:
1. **3D-Helmholtz** erzeugt die Modenlandschaft.
2. **Robin \(eta\)** selektiert die Ordnungsfamilie.
3. **Geometrische Gleichwertigkeit der drei Raumrichtungen** entscheidet, ob die Ordnung eher
   - volumetrisch finite-\(q\),
   - flächengetragen,
   - oder anisotrop gebrochen erscheint.

### Neue Leitfrage
> Wie verändert sich die Ordnungsfamilie in 3D-Helmholtz-Robin-Systemen, wenn die Gleichwertigkeit der drei Raumrichtungen systematisch gebrochen wird?

Damit verschiebt sich die Forschung von der Fall-zu-Fall-Lesung zu einer **Kontinuations- und Übergangsfrage**.

## Phase OP2 – Erster Kontinuations-Pilot
Der erste konkrete Testlauf untersucht vier Geometrien:
- cubic \((1,1,1)\)
- tetragonal_mild \((1,1,1.2)\)
- tetragonal_strong \((1,1,1.5)\)
- orthorhombic \((1,1.2,1.5)\)

für
- \(eta = 0, 1, 5\)

und liest jeweils die lokal stärkste `xyz`-tragende Ordnung auf einer **3×3×3-Superzellen-Logik** aus.

### Historisch genannte Artefakte
- `geometry_continuation_study.py`
- `geometry_continuation_study_first_test.csv`

### Ausgewertete Größen
- dominanter \(q\)-Vektor
- \(q\)-Kontrast
- mittlere lokale `xyz`-Stärke
- Feld-Anisotropie
- Randanteil

### Erster Pilotbefund
Dieser erste Lauf trägt bereits den Kern des späteren Ordnungsprinzips:
- kubische Gleichwertigkeit hält gemischtere Familien offen,
- milde tetragonale Verzerrung kippt früh in achsengebundene Lesarten wie `Z1`,
- größere Robin-Kopplung lässt wieder gemischtere Familien erscheinen,
- der Randanteil sinkt mit wachsendem \(eta\).

Wichtig ist die historische Einordnung: Der Lauf wird selbst ausdrücklich als **Pilot einer diskreten Robin-Surrogat-Continuation** markiert, also als belastbarer Zwischenstand, nicht als endgültiger Phasensatz.

## Phase OP3 – Erweiterte Geometrie-Kontinuation und zweites Auswahlkriterium
Im nächsten Schritt wird der Pilot gezielt verschärft:
- mehr Geometriepunkte zwischen kubisch und orthorhombisch,
- dasselbe Problem mit **zwei** Auswahlkriterien,
- daraus eine erste **Übersichtskarte**.

### Historisch genannte Artefakte
- `geometry_phase_map.py`
- `geometry_phase_map_full.csv`
- `geometry_phase_map_summary.csv`
- `geometry_phase_map_compact.csv`
- `geometry_phase_map_overview.png`

### Neue Geometriefolge
- cubic
- tetragonal: \(1.05, 1.10, 1.20, 1.35, 1.50\)
- orthorhombisch: \((1,1.05,1.20)\), \((1,1.10,1.30)\), \((1,1.20,1.50)\)

### Zwei Auswahlkriterien
1. `score_q = mean_abs_xyz × q-Kontrast`
2. `score_iso = mean_abs_xyz × Isotropie`

### Historische Konsequenz
Gerade dieser Schritt macht den Nebenast stark:
- die Kriterien sind nicht identisch,
- aber erzählen im Großen dieselbe Geschichte,
- insbesondere für \(eta=0\) und \(eta=1\) dominiert außerhalb des kubischen Punkts fast sofort `Z1`,
- während bei \(eta=5\) wieder gemischtere Familien wie `X1+Y1+Z1`, `X1+Y1+Z2`, `X1+Y2+Z2` oder `X1+Y2+Z1` auftreten.

Damit wird aus der bloßen Vermutung erstmals eine **Umschaltstruktur**.

## Phase OP4 – Kompakte Ordnungsphasenkarte
Danach wird die Kontinuation explizit in die Form einer **kompakten Ordnungsphasenkarte** gebracht:
- Geometrieindex,
- \(eta\),
- dominante Ordnungsfamilie,
- Stabilität zwischen beiden Kriterien.

Historisch ist das der Übergang von „sehr interessante numerische Tendenz“ zu einem **lesbaren Ergebnis** mit stabilen Zonen und Übergangsbereichen.

Die Kurzlesung dieses Schritts lautet:
- **kubisch** hält gemischtere volumetrische Familien offen,
- schon **milde tetragonale Anisotropie** kippt für kleine und moderate \(eta\) fast sofort in achsengebundene Ordnung,
- **größeres Robin** lockert diese Bindung teilweise wieder auf.

## Phase OP5 – Begriffliche Fixierung: Ordnungsphasen-Notiz
Auf Basis dieser Kontinuation wird anschließend eine explizite **Ordnungsphasen-Notiz** formuliert.

Dort wird der Kern erstmals begrifflich festgehalten:

> Kubische Richtungs-Gleichwertigkeit hält gemischte volumetrische Ordnungsfamilien offen; geometrische Anisotropie bevorzugt achsengebundene Ordnung; Robin selektiert und moduliert den Umschlag zwischen beiden.

Wichtig ist der methodische Zusatz:
- kein endgültiger Phasensatz,
- sondern ein belastbarer Arbeitsstand,
- gestützt durch zwei unterschiedliche Auswahlkriterien.

## Phase OP6 – Prinzipien-Notiz als übergreifende Verdichtung
Am Ende steht die **Prinzipien-Notiz**. Sie ist keine neue numerische Phase, sondern die erste echte **Synthese** des ganzen Paper-1-Kerns.

Dort werden in kondensierter Form zusammengezogen:
1. 3D-Helmholtz als Strukturträger,
2. Robin als Selektionsparameter,
3. Volumenprimat gegenüber natürlichen Flächen-Observablen,
4. hexaedrisch/kubische Klasse als stärkster Träger,
5. primitive kubisch als modellnäheste Lesung,
6. Richtungs-Gleichwertigkeit als Ordnungsfaktor,
7. eine erste physikalisch lesbare Signatur als Mehrkomponenten-Diagnostik.

Damit wird aus vielen Einzelbefunden erstmals ein **explizites wiederkehrendes Ordnungsprinzip**.

## Phase OP7 – Einordnung in Paper 1
Für Paper 1 ist dieser Nebenast kein neuer technischer Hauptstrang, sondern eine **explizite Syntheseschicht**.

Seine Funktion ist:
- die vielen Diagnosen nicht nur nebeneinander stehen zu lassen,
- sondern ihre gemeinsame Logik sichtbar zu machen,
- ohne den Haupttext des Papers mit allen Kontinuationsdetails zu überladen.

Gerade deshalb eignet sich dieser Ast im Repo besonders gut als:
- Forschungshistorie,
- Statusnotiz,
- später eventuell als Appendix-/Supplement-Modul,
- aber noch nicht als vollständig restaurierter lauffähiger Kern.

---

## Was für das Repository daraus folgt
Dieser Nebenast sollte im Repo vorerst **dokumentationsstark und restaurationsoffen** geführt werden:

1. **bereits vorhanden als Quellenbasis**
   - geordnete Chat-Exporte
   - Prinzipien-Notiz
   - Ordnungsphasen-Notiz

2. **historisch belegt, aber noch nicht dateiseitig restauriert**
   - Kontinuationsskript
   - Pilot-CSV
   - Phase-Map-Skript
   - volle / kompakte CSV
   - Übersichtsbild

3. **Dokumentationsschicht**
   - Forschungshistorie
   - Repo-Unterstruktur
   - Status-/Provenienznotiz

## Offene Wiederherstellungen
Für einen dateiseitig vollständigen Nebenast fehlen aktuell diese historisch belegten Artefakte:
- `geometry_continuation_study.py`
- `geometry_continuation_study_first_test.csv`
- `geometry_phase_map.py`
- `geometry_phase_map_full.csv`
- `geometry_phase_map_summary.csv`
- `geometry_phase_map_compact.csv`
- `geometry_phase_map_overview.png`

## Empfehlung
Für die öffentliche GitHub-Fassung sollte der Nebenast jetzt so erscheinen:
- **klare dokumentierte Forschungshistorie**,
- **ehrlich markierte noch fehlende numerische Artefakte**,
- **keine vorgetäuschte Reproduzierbarkeit**, solange die historischen Skripte und CSV/PNG-Dateien nicht wiedergefunden oder kuratiert rekonstruiert sind.

So bleibt das Repo sauber und vollständig zugleich.
