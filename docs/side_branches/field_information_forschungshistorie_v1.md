# Feldinformationsdiagnose-Nebenast – Forschungshistorie v1

## Zweck
Diese Notiz rekonstruiert den **Feldinformationsdiagnose-Nebenast** als geordneten Forschungsverlauf innerhalb von Paper 1. Sie dient nicht als fertiges Manuskript, sondern als interne Provenienz-, Struktur- und Restaurationshilfe für das Repository.

## Quellenbasis
Verwendet wurden die geordneten Chat-Exporte
- `chatexport_field_1_erste Fragestellung und Vorschläge.txt`
- `chatexport_field_projektkarte.txt`
- `chatexport_field_2_erster Testlauf und befunde.txt`
- `chatexport_field_3_.txt`

Zusätzlich wurde der bereits wiederhergestellte Code-/Datenstand `field_information_distribution.py` und `field_information_distribution_first_test.csv` berücksichtigt.

## Historischer Kernbefund in einem Satz
Die Fragestellung verschob sich von der naiven Erwartung „volle 3D = zugleich isotropste **und** breiteste Verteilung" zu der präziseren Aussage:

> Verschiedene Moden- und Gitterklassen unterscheiden sich reproduzierbar in ihrer **Informationsgeometrie**; volle 3D erscheint dabei als **isotropste**, aber nicht notwendigerweise als **delokalisierteste** Klasse, und primitive vs. fcc trennen sich eher über **Subgitter-/Ordnungslesung** als über grobe Vollvolumenverteilung.

---

## Phase FI0 – Formulierung des letzten Diagnoseasts
Der Nebenast entsteht aus einer klaren Anschlussfrage an den Hauptstrang und an die Flächendiagnose:

- Welche Ordnung liegt vor?
- Auf welchem Subgitter bzw. Träger sitzt sie?
- Wie unvollständig tragen Flächen sie?
- **Und wie verteilt sich diese Ordnung überhaupt im 3D-Volumen?**

Damit verschiebt sich der Fokus von der bloßen Klassenfrage hin zu einer **Geometrie der Informationsverteilung**.

Wichtig ist die methodische Einordnung:
- keine neue Theorie,
- keine ontologische Zusatzbehauptung,
- sondern ein messbarer Zusatzast innerhalb des bestehenden 3D-Helmholtz-Robin-Rahmens.

## Phase FI1 – Erste Arbeitshypothese und Kennzahlenvorschlag
Im nächsten Schritt wird der Ast offen und testbar formuliert.

### Leitfrage
Wie ist die relevante Feldinformation im Volumen verteilt?

Nicht nur:
- achsial,
- paarartig,
- volle 3D,
- primitive kubisch,
- fcc,

sondern:
- breit volumetrisch,
- auf wenige Hotspots konzentriert,
- linienartig,
- flächenartig,
- oder als echte 3D-Gitterverteilung.

### Arbeitshypothese
Die relevante Ordnung ist volumetrisch primär, aber nicht notwendigerweise homogen verteilt. Unterschiedliche Moden- und Gitterklassen sollten sich durch
- Delokalisierung,
- räumliche Konzentration,
- geometrische Anisotropie

unterscheiden.

### Vorgeschlagene Kennzahlen
Der Nebenast wird bewusst über Standardgrößen aufgesetzt:
- `f_eff` als effektive Volumenbeteiligung / Participation-Ratio-Lesung,
- `H_norm` als normierte Shannon-Entropie,
- Eigenwerte des zweiten Momenttensors als Form-/Anisotropie-Lesung,
- optional `B` bzw. Randanteil als Anschluss an den Flächennachast.

Diese Wahl ist methodisch wichtig, weil der Ast damit **misst statt behauptet**.

## Phase FI2 – Projektkarte und erste Zielklassen
Danach wird der Nebenast in einer kleinen Projektkarte fixiert.

### Phase 1 – Einzelzelle
- Achsenmode
- Paarmode
- volle 3D-Mode

### Phase 2 – Gitterebene
- primitive finite-`q`-Ordnung
- fcc-flächengetragene Ordnung

### Erwartungsrichtung
Historisch wird anfangs noch erwartet:
- Achse: anisotroper, eher weniger volumetrisch,
- Paar: mittlere Verteilung,
- volle 3D: isotroper und eventuell volumetrisch breiter,
- primitive finite-`q`: eher volumetrisch,
- fcc: eher subgitter-/flächengebundener.

Gerade weil diese Erwartung später **nur teilweise** bestätigt wird, ist die Projektkarte ein wichtiger Ausgangspunkt.

## Phase FI3 – Erster Würfeltest auf Achse / Paar / volle 3D
Danach wird der erste reale Testlauf auf dem **Würfel** gefahren. Historisch dazu gehören das wiedergefundene Skript `field_information_distribution.py` und die Referenzdatei `field_information_distribution_first_test.csv`.

### Testaufbau
- `beta = 0, 1, 5`
- je eine repräsentative Achsen-, Paar- und volle 3D-Mode
- Auswertung von `f_eff`, `H_norm`, Anisotropie und Randanteil

### Tragender erster Befund
Der erste Test ist gerade deshalb stark, weil er der naiven Erwartung **nicht einfach folgt**:
- volle 3D ist durchgehend die **isotropste** Klasse,
- aber **nicht** die breiteste im Sinn von `f_eff`/Entropie,
- Achse und Paar bleiben geometrisch klar unterscheidbar,
- mit wachsendem `beta` sinkt zusätzlich der Randanteil.

Historisch kippt die Lesung hier von
- „volle 3D ist wohl am stärksten überall verteilt"
zu
- „volle 3D ist vor allem **isotrop**, nicht automatisch **maximal delokalisiert**“.

Das ist der eigentliche Startpunkt des Nebenasts als **ehrliche Informationsgeometrie**.

## Phase FI4 – Übergang zur Superzelle
Im nächsten Schritt wird explizit gefragt, ob dieselben Kennzahlen auch auf die Gitterebene übertragen werden können:
- primitive finite-`q`-Superzelle,
- fcc-flächengetragene Ordnung,
- Vergleich auf Basis derselben Verteilungsmaße.

Damit verlässt der Ast die Einzelzelle und wird an die bereits vorher aufgebaute Superzellenlogik des Hauptstrangs angeschlossen.

## Phase FI5 – Superzellentest primitive vs. fcc
Die nächste historische Stufe ist der 3×3×3-Superzellentest.

### Historische Lesung
Das Ergebnis ist in zweierlei Hinsicht aufschlussreich:

1. Für `beta = 1` und `beta = 5` liegen primitive und fcc in dieser Lesung auf **derselben zugrunde liegenden Eigenmode**. Die globale Vollvolumenverteilung fällt darum identisch aus.

2. Für `beta = 0` unterscheiden sich primitive und fcc zwar, aber nur **schwach** in `f_eff`, Entropie, Anisotropie und Randanteil.

### Konsequenz
Genau hier wird die Erwartung nochmals geschärft:
- primitive vs. fcc trennen sich **nicht zwingend** über grobe Vollvolumen-Dichtegeometrie,
- sondern offenbar eher über
  - Subgitter-Bindung,
  - Ordnungsvektor,
  - Peakfamilie,
  - und die konkrete Ordnungslesung.

Das ist eine wichtige negative Präzisierung, keine Schwäche.

## Phase FI6 – Endformulierung des Nebenasts
Am Ende steht eine deutlich robustere Kurzfassung als am Anfang.

### Für Einzelzellen
Achsen-, Paar- und volle 3D-Moden sind nicht nur topologisch verschieden, sondern auch in ihrer räumlichen Informationsverteilung.

### Für volle 3D
Volle 3D ist die **isotropste**, aber nicht die **delokalisierteste** Verteilung.

### Für primitive vs. fcc
Primitive finite-`q`-Ordnung und fcc-Ordnung unterscheiden sich in diesem Ast **nicht primär** über grobe Vollvolumenverteilung, sondern eher über feinere Ordnungs- und Sublattice-Struktur.

Die geschärfte Endform lautet daher:

> Nicht nur die Art der Ordnung, sondern auch der Grad und die Geometrie ihrer räumlichen Verteilung sind im 3D-Helmholtz-Robin-System strukturierte Größen; diese Verteilungsgeometrie trennt Klassen reproduzierbar, aber nicht immer in der anfangs erwarteten Richtung.

## Phase FI7 – Einordnung in Paper 1
Für Paper 1 bleibt dieser Ast ein **Nebenast**. Er ist diagnostisch wertvoll, aber nicht der Hauptträger der Kernargumentation.

Seine Stärke liegt darin, dass er
- Volumen-vs.-Fläche,
- Modenklassifikation,
- primitive vs. fcc,
- und die Frage nach lokalisierter vs. breit verteilter Information

in einer einzigen Zusatzdiagnostik zusammenführt.

---

## Was für das Repository daraus folgt
Der Nebenast sollte im Repo **nicht** als neuer Hauptstrang auftreten, sondern als klar markiertes Zusatzmodul:

1. **restaurierter Kern**
   - erster Würfeltest mit `field_information_distribution.py`
   - historisches Referenz-CSV `field_information_distribution_first_test.csv`

2. **historisch belegt, aber noch nicht dateiseitig restauriert**
   - `field_information_distribution_supercell.py`
   - `field_information_distribution_supercell_test.csv`

3. **Dokumentationsschicht**
   - Forschungshistorie
   - Repo-Unterstruktur
   - Status-/Provenienznotiz

## Offene Wiederherstellungen
Für einen dateiseitig vollständigen Nebenast fehlen aktuell noch diese historisch belegten Artefakte:
- `field_information_distribution_supercell.py`
- `field_information_distribution_supercell_test.csv`

## Empfehlung
Für die öffentliche GitHub-Fassung sollte der Nebenast jetzt so erscheinen:
- **ein restaurierter, direkt lauffähiger erster Würfeltest**,
- **eine klare historische Einordnung**,
- **ehrlich markierte noch fehlende Superzellen-Artefakte**.

So bleibt das Repo sauber, ohne die tatsächliche Forschungsgeschichte zu verkürzen.
