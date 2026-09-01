import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
# ── Configuration ─────────────────────────────────────────────────────────────
import glob as _glob
_matches = _glob.glob("./OUTRIDER_*_annote.xlsx")
if not _matches:
    raise FileNotFoundError("Aucun fichier OUTRIDER_*_annote.xlsx trouvé dans le répertoire courant.")
if len(_matches) > 1:
    raise ValueError(f"Plusieurs fichiers correspondants : {_matches}. Précisez lequel utiliser.")
INPUT_FILE = _matches[0]
OUTPUT_DIR   = "./outrider_volcano_plots"    # dossier de sortie pour les PNG
ZSCORE_COL   = "zScore"
PVALUE_COL   = "pValue"
GENE_COL     = "gene_symbol"
PATIENT_COL  = "sampleID"
# Seuils pour la mise en couleur
ZSCORE_THRESH  = 2     # |zScore| >= seuil → significatif
#PVALUE_THRESH  = 1.3     # -log10(pValue) >= seuil → significatif (ex: 1.3 ≈ p=0.05)
PVALUE_THRESH  = 3.5     # -log10(pValue) >= seuil → significatif (ex: 1.3 ≈ p=0.05)
# ──────────────────────────────────────────────────────────────────────────────


def classify(row, neg_log10_p):
    sig_z = abs(row[ZSCORE_COL]) >= ZSCORE_THRESH
    sig_p = neg_log10_p >= PVALUE_THRESH

    if sig_z and sig_p:
        return "sig"   # rouge
    else:
        return "ns"    # gris


COLOR_MAP = {
    "sig": "red",
    "ns":  "grey",
}

LABEL_MAP = {
    "sig": f"|Z|>={ZSCORE_THRESH} & -log10(p)>={PVALUE_THRESH}",
    "ns":  "Non significatif",
}


def volcano_for_patient(df_patient, patient_id, out_dir):
    df = df_patient.copy()

    # Supprimer les lignes sans p-value ou zScore valides
    df = df.dropna(subset=[ZSCORE_COL, PVALUE_COL, GENE_COL])
    df = df[df[PVALUE_COL] > 0]  # éviter log(0)

    df["-log10p"] = -np.log10(df[PVALUE_COL])
    df["category"] = df.apply(lambda row: classify(row, row["-log10p"]), axis=1)

    fig, ax = plt.subplots(figsize=(10, 7))

    # Tracer les points par catégorie
    for cat, grp in df.groupby("category"):
        ax.scatter(
            grp[ZSCORE_COL],
            grp["-log10p"],
            c=COLOR_MAP[cat],
            alpha=0.75,
            edgecolors="none",
            s=40,
            label=LABEL_MAP[cat],
            zorder=2,
        )

    # Étiquettes des gènes (tous les points colorés)
    labeled = df[df["category"] != "ns"]
    for _, row in labeled.iterrows():
        ax.annotate(
            row[GENE_COL],
            xy=(row[ZSCORE_COL], row["-log10p"]),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=6,
            color="#222222",
            ha="left",
        )

    # Lignes de seuil
    ax.axvline( ZSCORE_THRESH, color="#888", lw=0.8, ls="--", zorder=1)
    ax.axvline(-ZSCORE_THRESH, color="#888", lw=0.8, ls="--", zorder=1)
    ax.axhline(PVALUE_THRESH, color="#888", lw=0.8, ls="--", zorder=1)

    # Échelles centrées sur les données avec marge de 10 %
    pad_x = (df[ZSCORE_COL].max() - df[ZSCORE_COL].min()) * 0.10 or 0.5
    pad_y = (df["-log10p"].max() - df["-log10p"].min()) * 0.10 or 0.5
    ax.set_xlim(df[ZSCORE_COL].min() - pad_x, df[ZSCORE_COL].max() + pad_x)
    #ax.set_ylim(df["-log10p"].min() - pad_y, df["-log10p"].max() + pad_y)
    ax.set_ylim(0.5)

    ax.set_xlabel("Z-score", fontsize=12)
    ax.set_ylabel("-log₁₀(p-value)", fontsize=12)
    ax.set_title(f"Volcano plot — Patient {patient_id}", fontsize=13, fontweight="bold")
    ax.legend(fontsize=8, framealpha=0.8, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    safe_id = str(patient_id).replace("/", "_").replace(" ", "_")
    out_path = os.path.join(out_dir, f"volcano_{safe_id}.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  Sauvegardé : {out_path}")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print(f"Lecture de {INPUT_FILE} …")
    df = pd.read_excel(INPUT_FILE)

    patients = df[PATIENT_COL].unique()
    print(f"{len(patients)} patient(s) trouvé(s) : {list(patients)}\n")

    for pid in patients:
        print(f"→ Patient : {pid}")
        volcano_for_patient(df[df[PATIENT_COL] == pid], pid, OUTPUT_DIR)

    print(f"\nTerminé. Les plots sont dans le dossier « {OUTPUT_DIR}/ ».")


if __name__ == "__main__":
    main()
