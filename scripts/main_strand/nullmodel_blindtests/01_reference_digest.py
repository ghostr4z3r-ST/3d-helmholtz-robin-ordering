from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
REF = ROOT / "results" / "main_strand" / "nullmodel_blindtests" / "reference"
OUT = ROOT / "results" / "main_strand" / "nullmodel_blindtests" / "reference_digest.md"


def read_csv(name: str):
    path = REF / name
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def fmt(x: str) -> str:
    try:
        return f"{float(x):.3f}"
    except Exception:
        return str(x)


def main() -> None:
    fa = read_csv("feature_ablation_summary_reduced.csv")
    cube = read_csv("field_shuffle_null_tests_50.csv")
    sc = read_csv("field_shuffle_supercell_null_tests.csv")

    lines = []
    lines.append("# Nullmodell-/Blindtest-Referenzdigest")
    lines.append("")
    lines.append("Diese Datei fasst **mitgelieferte Referenz-CSV** zusammen. Sie reproduziert die historischen Tests nicht neu, sondern bietet einen kuratierten Leseeinstieg in Phase J des Hauptstrangs.")
    lines.append("")
    lines.append("## Feature-Ablation (reduzierter erster Lauf)")
    for row in fa:
        lines.append(f"- {row['dataset']}: voll={fmt(row['full_accuracy'])}, wichtigste Gruppe={row['most_informative_group']}, beste Einzelgruppe={row['best_single_group']} ({fmt(row['best_single_accuracy'])})")
    lines.append("")
    lines.append("## Feld-Shuffle auf Würfel-Aufgaben")
    for row in cube:
        lines.append(f"- {row['task']}: echt={fmt(row['A_real'])}, Shuffle={fmt(row['shuffle_mean'])} ± {fmt(row['shuffle_std'])}, p={fmt(row['p_ge'])}")
    lines.append("")
    lines.append("## Feld-Shuffle auf Superzelle / q-Diagnostik")
    for row in sc:
        if row['section'] == 'classification':
            lines.append(f"- {row['task']}: echt={fmt(row['real'])}, Shuffle={fmt(row['shuffle_mean'])} ± {fmt(row['shuffle_std'])}, p={fmt(row['p_ge'])}")
    lines.append("")
    lines.append("## q-Metriken mit positiver Shuffle-Robustheit")
    for row in sc:
        if row['section'] == 'q_metrics':
            z = float(row['z_score'])
            if z > 0:
                lines.append(f"- {row['task']}: echt={fmt(row['real'])}, Shuffle={fmt(row['shuffle_mean'])} ± {fmt(row['shuffle_std'])}, z={fmt(row['z_score'])}, p={fmt(row['p_ge'])}")
    lines.append("")
    OUT.write_text("\n".join(lines), encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
