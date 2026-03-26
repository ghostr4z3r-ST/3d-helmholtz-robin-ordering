# Flächendiagnose-Nebenast – Forschungshistorie v1

## Zweck
Diese Notiz rekonstruiert den **Flächendiagnose-Nebenast** als geordneten Forschungsverlauf innerhalb von Paper 1. Sie dient nicht als fertiges Paper-Manuskript, sondern als interne Provenienz- und Strukturhilfe für das Repository.

## Quellenbasis
Verwendet wurden die geordneten Chat-Exporte `chatexport_f_1` bis `chatexport_f_9` sowie die Kurznotiz zur Umformulierung der Arbeitshypothese. Zusätzlich wurde der bereits wiederhergestellte Code-/Datenstand (`face_information_sidebranch.py`, `face_information_first_test.csv`) berücksichtigt.

## Historischer Kernbefund in einem Satz
Die Fragestellung verschob sich von einer starken Erwartung der Form `eine Fläche < alle Flächen < Volumen` zu der robusteren Aussage:

> Die relevante Ordnung ist im untersuchten Modellrahmen **volumetrisch primär**; natürliche Flächen-Observablen tragen sie nur **partiell** und bleiben deutlich hinter Volumen-/Gitter-Observablen zurück.

---

## Phase F0 – Von alter Intuition zu sauberer Frage
Am Anfang steht nicht mehr die große Behauptung „Fläche trägt die Wirklichkeit“, sondern eine mathematisch prüfbare Frage: **Wann und in welchem Maß ist Volumenordnung auf natürlichen Flächen-Observablen lesbar?**

Wichtig ist die methodische Verschiebung:
- weg von einer ontologischen Aussage,
- hin zu einer testbaren Diagnosefrage innerhalb des 3D-Helmholtz-Robin-Rahmens.

Damit wird der Nebenast direkt anschlussfähig an den Hauptstrang: Volumenmoden, Robin-Selektion, Würfel-/Superzellen-Readouts und bereits vorhandene Volumenlabels werden als Referenz genutzt.

## Phase F1 – Projektkarte und erste Arbeitshypothese
Danach wird der Nebenast als eigene kleine Projektkarte formuliert.

### Leitfrage
Wie viel Information über die **Volumenklasse** einer Mode steckt in
- einer einzelnen Fläche,
- mehreren Flächen gemeinsam,
- dem gesamten Volumen?

### Erste Arbeitshypothese
Die Rekonstruierbarkeit steigt mit der Zahl der Flächen, bleibt aber hinter der Volumenlesbarkeit zurück.

### Kennzahlen
- `A_1F`: Klassifikationsgüte aus **einer** Fläche
- `A_6F`: Klassifikationsgüte aus **allen sechs** Flächen
- `A_V`: Klassifikationsgüte aus Volumen-/Gitterfeatures
- `D_1F`, `D_6F`: Defizite der Flächenlesbarkeit relativ zum Volumen

### Zielklassen
Die ersten Zielklassen werden ausdrücklich aus dem Hauptstrang übernommen:
- Achsenmode
- Paarmode
- volle 3D-Mode
- später primitive finite-`q`-Ordnung
- später fcc-flächengetragene Ordnung

Damit ist der Nebenast **nicht isoliert**, sondern klar als Diagnosemodul über einem bereits aufgebauten Volumenkern angelegt.

## Phase F2 – Erster Würfeltest mit Face-Features
Anschließend wird der erste echte Testlauf auf dem **Würfel** aufgesetzt. Historisch dazu gehören der wiederhergestellte Code `face_information_sidebranch.py` und die Referenzdatei `face_information_first_test.csv`.

### Inhalt des ersten Tests
- Betas: `0, 1, 5`
- Klassen: Achsen-, Paar- und volle 3D-Moden
- Vergleich von `A_1F`, `A_6F`, `A_V`
- Flächenfeatures auf Basis einfacher lokaler Signaturen/Hadamard-artiger 2D-Projektionen

### Historische Lesung
Der erste Lauf trägt bereits die Grundrichtung:
- Flächen sind **schwächer** als Volumen,
- mehrere Flächen können hilfreich sein,
- aber der Volumenreadout bleibt die Referenz.

Dieser Schritt ist die eigentliche Geburt des Nebenasts als numerischer Diagnose.

## Phase F3 – Erweiterte Features, mehr Moden, erste Robustheit
Im nächsten Schritt wird der Nebenast methodisch verschärft:
- reichere Flächenfeatures,
- mehr Moden pro Klasse,
- Würfeltests über mehrere `beta`-Werte,
- erste Robustheitsprüfung über die Auflösung,
- Übertragung auf primitive Superzelle und fcc-Gegenlesung.

### Historische Konsequenz
Hier kippt die Interpretation von einer starken zu einer schwächeren, aber robusteren Form:
- **nicht robust:** `A_1F < A_6F` als allgemeines Gesetz,
- **robust:** Flächen-Observablen bleiben insgesamt deutlich schwächer als Volumen-/Gitter-Observablen.

Dieser Schritt ist wissenschaftlich wichtig, weil der Nebenast dadurch **härter** statt bloß größer wird.

### Backup-Status
Die in der Chat-Historie genannten Artefakte zu diesem Schritt sind bisher **nicht** als Originaldateien wiedergefunden:
- `face_information_sidebranch_extended.py`
- `face_information_extended_test.csv`
- `face_information_extended_cube_robustness.csv`

Sie sind daher momentan **historisch belegt**, aber noch **nicht dateiseitig restauriert**.

## Phase F4 – Verschärfung zur Patch-Scan-Frage
Danach wird die Fragestellung nochmals verfeinert:

> Vielleicht ist die Ordnung nicht vollständig auf Flächen lesbar, aber lokal-statistisch über Patch-Scans wenigstens als Familienstruktur erkennbar.

Dafür werden neue Größen eingeführt:
- `A_1F^scan`
- `A_6F^scan`

und lokale Patch-Features wie
- Mittelwert,
- Streuung,
- Vorzeichenbalance,
- kleine 2D-Hadamard-Projektionen,
- einfache Gradienten-/Checkerboard-Proxys.

Methodisch ist das der Übergang von „globale Flächenlesung“ zu „statistischer Flächenfingerabdruck“.

## Phase F5 – Erster Patch-Scan-Test auf dem Würfel
Der erste Patch-Scan-Test läuft wieder auf dem Würfel. Das Ergebnis ist historisch gerade deshalb wichtig, weil es **negativ** ausfällt:
- die Scan-Varianten verbessern die Flächenlesbarkeit nicht zuverlässig,
- teils fallen sie sogar hinter die grobe Flächenlesung zurück,
- das Volumen bleibt klar überlegen.

Damit wird die Hypothese erneut verschärft:

> Selbst einfaches lokales Patch-Scanning rekonstruiert die Volumenordnung nicht robust.

### Backup-Status
Die dazu genannten Dateien sind derzeit nur in der Chat-Historie belegt, aber nicht als Originaldateien vorhanden:
- `face_patch_scan_sidebranch.py`
- `face_patch_scan_first_test.csv`

## Phase F6 – Patch-Scan auf Superzelle und fcc
Der nächste Schritt prüft dieselbe Frage auf der **Gitterebene**:
- primitive kubische Superzelle,
- fcc-Lesung,
- erneut Vergleich von grober Flächenlesung, Patch-Scan und Volumen-/Gitterreadout.

Die historische Schlusslesung lautet:
- Einzelflächen können schwache statistische Spuren tragen,
- aber der Patch-Scan „rettet“ die Flächen nicht,
- auch auf der Gitterebene bleibt die relevante Ordnung volumetrisch bzw. gitterseitig primär.

### Backup-Status
Noch nicht dateiseitig wiedergefunden:
- `face_patch_scan_supercell.py`
- `face_patch_scan_supercell_test.csv`

## Phase F7 – Selbsttest und methodische Einordnung
Danach erfolgt ein wichtiger Selbsttest des ganzen Nebenasts.

Die saubere Einordnung lautet:
- **kein bloßer Zirkelschluss**, weil die Tests anders hätten ausfallen können und mehrere starke Erwartungen gerade **nicht** stabil bestätigt wurden,
- aber auch **kein modellunabhängiger Beweis**, weil alles innerhalb eines konkret definierten Modellrahmens gilt.

Diese Phase ist wichtig für die spätere Sprache im Paper:
- keine metaphysische Weltbehauptung,
- sondern ein **numerisch-mathematisches Ergebnis innerhalb des 3D-Helmholtz-Robin-Rahmens**.

## Phase F8 – Endformulierung des Nebenasts
Am Ende steht die geschärfte Kurzfassung:

> Die relevante Ordnung ist volumetrisch organisiert; natürliche Flächen-Observablen tragen sie nur partiell und erlauben keine vollständige Rekonstruktion der Volumenklasse.

Oder noch knapper:

> **Die Ordnung ist volumetrisch, flächenweise nur partiell lesbar.**

Diese Formulierung ersetzt ausdrücklich die frühere stärkere Erwartung, dass die Gesamtheit aller Flächen stets klar gegen Einzelflächen gewinnt.

---

## Was für das Repository daraus folgt
Der Nebenast sollte im Repo **nicht** als zweiter Hauptstrang auftreten, sondern als klar markierte Unterkategorie innerhalb von Paper 1:

1. **restaurierter Kern**
   - erster Face-Information-Test auf dem Würfel
   - historisches Referenz-CSV

2. **historisch belegt, noch nicht restauriert**
   - erweiterte Features / Robustheit
   - Patch-Scan auf Würfel
   - Patch-Scan auf Superzelle/fcc

3. **Dokumentationsschicht**
   - Forschungshistorie
   - Substruktur-Plan
   - Manifest fehlender historischer Artefakte

## Offene Wiederherstellungen
Für einen dateiseitig vollständigen Nebenast fehlen aktuell noch diese historischen Artefakte:
- `face_information_sidebranch_extended.py`
- `face_information_extended_test.csv`
- `face_information_extended_cube_robustness.csv`
- `face_patch_scan_sidebranch.py`
- `face_patch_scan_first_test.csv`
- `face_patch_scan_supercell.py`
- `face_patch_scan_supercell_test.csv`

## Empfehlung
Für die öffentliche GitHub-Fassung sollte der Nebenast jetzt so erscheinen:
- **ein restaurierter, direkt lauffähiger erster Test**,
- **eine klare historische Einordnung**,
- **ehrlich markierte noch fehlende historische Zwischenartefakte**.

So bleibt das Repo sauber, ohne die tatsächliche Forschungsgeschichte zu verlieren.
