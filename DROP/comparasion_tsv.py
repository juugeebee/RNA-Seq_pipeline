import csv
import sys

file1 = './results_fraser2_run1_drop1.5.tsv'
file2 = './Moabi_fraser2mRNA-D_NOVASEQ-Run1.tsv'

key_col = 0  # colonne clé (toujours comparée)

out_header_diff = "header_diff.tsv"
out_only_1 = "only_in_fichier_drop1.5.tsv"
out_only_2 = "only_in_fichier_drop1.4tsv"
out_diff = "differences.tsv"


def read_tsv(filename):
    with open(filename, newline="", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        rows = {row[key_col]: row for row in reader}
    return header, rows


# --- Lecture ---
header1, data1 = read_tsv(file1)
header2, data2 = read_tsv(file2)

# --- Détection des colonnes différentes ---
ignored_cols = set()

with open(out_header_diff, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["index", "fichier1", "fichier2"])

    max_len = max(len(header1), len(header2))
    for i in range(max_len):
        h1 = header1[i] if i < len(header1) else ""
        h2 = header2[i] if i < len(header2) else ""
        if h1 != h2:
            writer.writerow([i, h1, h2])
            ignored_cols.add(i)

# ⚠️ la colonne clé ne doit jamais être ignorée
ignored_cols.discard(key_col)

print(f"Colonnes ignorées : {sorted(ignored_cols)}")

# --- Fonction de comparaison sans les colonnes ignorées ---
def filtered_row(row):
    return [v for i, v in enumerate(row) if i not in ignored_cols]


# --- Comparaison des données ---
keys1 = set(data1)
keys2 = set(data2)

only_1 = keys1 - keys2
only_2 = keys2 - keys1
common = keys1 & keys2

# --- Sorties ---
with open(out_only_1, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header1)
    for k in sorted(only_1):
        writer.writerow(data1[k])

with open(out_only_2, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(header1)
    for k in sorted(only_2):
        writer.writerow(data2[k])

with open(out_diff, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["clé", "source", *header1])

    for k in sorted(common):
        r1 = filtered_row(data1[k])
        r2 = filtered_row(data2[k])
        if r1 != r2:
            writer.writerow([k, "fichier1.5", *data1[k]])
            writer.writerow([k, "fichier1.4", *data2[k]])

print("✅ Comparaison terminée (colonnes ignorées prises en compte)")