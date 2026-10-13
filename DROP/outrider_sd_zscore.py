#!/usr/bin/env python3
"""
Ajoute les colonnes mean_l2fc, sd_l2fc et z_l2fc à chacun des fichiers OUTRIDER annotés.
Les originaux ne sont pas modifiés : une copie *_zscore.xlsx est écrite dans OUT_DIR
(mêmes feuilles, mêmes colonnes + les nouvelles colonnes).

Usage : renseigner FILES ci-dessous, puis  python outrider_sd_zscore.py
"""
import os
import sys

import pandas as pd

# ---------------------------------------------------------------- À ÉDITER
FILES = [
    "/run/user/1000/gvfs/smb-share:server=wprod.ds.aphp.fr,share=dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev5/_RNASEQ/RNASEQ_2025/Il-str-mRNA-D_NOVASEQ-Run1_17072025/03_Analyse/OUTRIDER_annote 2.xlsx",
    "/run/user/1000/gvfs/smb-share:server=wprod.ds.aphp.fr,share=dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev5/_RNASEQ/RNASEQ_2025/Il-str-mRNA-D_NOVASEQ-Run2_05092025/03_Analyse/PSL/OUTRIDER_Il-str-mRNA-D_NOVASEQ-Run2_annote_01.xlsx",
    "/run/user/1000/gvfs/smb-share:server=wprod.ds.aphp.fr,share=dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev5/_RNASEQ/RNASEQ_2025/Il-str-mRNA-D_NOVASEQ-Run2_05092025/03_Analyse/PSL/OUTRIDER_Il-str-mRNA-D_NOVASEQ-Run2_annote_02.xlsx",
    "/run/user/1000/gvfs/smb-share:server=wprod.ds.aphp.fr,share=dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev5/_RNASEQ/RNASEQ_2025/Il-str-mRNA-D_NOVASEQ-Run3_10112025/03_Analyse/PSL/OUTRIDER_Il-str-mRNA-D_NOVASEQ-Run3_annote_01_PSL.xlsx",
    "/run/user/1000/gvfs/smb-share:server=wprod.ds.aphp.fr,share=dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev5/_RNASEQ/RNASEQ_2025/Il-str-mRNA-D_NOVASEQ-Run3_10112025/03_Analyse/PSL/OUTRIDER_Il-str-mRNA-D_NOVASEQ-Run3_annote_02_PSL.xlsx",
    "/run/user/1000/gvfs/smb-share:server=wprod.ds.aphp.fr,share=dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev5/_RNASEQ/RNASEQ_2025/Il-str-mRNA-D_NOVASEQ-Run4_13012026/03_Analyse/outrider/PSL/OUTRIDER_PAXonly_Il-str-mRNA-D_NOVASEQ-Run4_annote.xlsx",
    "/run/user/1000/gvfs/smb-share:server=wprod.ds.aphp.fr,share=dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev5/_RNASEQ/RNASEQ_2026/Il-str-mRNA-D_NOVASEQ-Run5_15042026/03_Analyse/outrider/OUTRIDER_Il-str-mRNA-D_NOVASEQ-Run5_annote.xlsx"
]
OUT_DIR = "/run/user/1000/gvfs/smb-share:server=wprod.ds.aphp.fr,share=dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev5/_RNASEQ/OUTRIDER_sortie_zscore"
# "global" : une SD unique sur toutes les lignes des fichiers annotes
# "gene"   : SD par gène (ensg) sur les 5 fichiers
# "source" : SD calculée séparément au sein de chaque fichier
MODE = "global"
# --------------------------------------------------------------------------

VALUE = "l2fc"  # déjà en log2 ; foldChange est arrondi (0.00 -> log2 = -inf)


def read_sheets(path, cols=None):
    """Toutes les feuilles du fichier -> {nom_feuille: DataFrame}."""
    return pd.read_excel(path, sheet_name=None, usecols=cols)


def compute_stats(files):
    """Passe 1 : statistiques calculées sur l'ensemble des fichiers (modes global / gene)."""
    frames = []
    for f in files:
        print(f"Lecture de {os.path.basename(f)} ...")
        frames += list(read_sheets(f, ["ensg", VALUE]).values())
    
    #frames = [d for f in files for d in read_sheets(f, ["ensg", VALUE]).values()]
    
    lean = pd.concat(frames, ignore_index=True).dropna(subset=[VALUE])
    if MODE == "global":
        mu, sd = lean[VALUE].mean(), lean[VALUE].std()
        print(f"Moyenne globale = {mu:.4f} ; SD globale = {sd:.4f} ({len(lean)} lignes)")
        return mu, sd
    return (lean.groupby("ensg")[VALUE]
                .agg(mean_l2fc="mean", sd_l2fc="std", n_obs="count")
                .reset_index())


def add_columns(sheets, stats):
    """Passe 2 : ajoute mean_l2fc, sd_l2fc, z_l2fc à chaque feuille."""
    if MODE == "source":
        allv = pd.concat([d[VALUE] for d in sheets.values()])
        mu, sd = allv.mean(), allv.std()
    elif MODE == "global":
        mu, sd = stats
    out = {}
    for name, d in sheets.items():
        if MODE == "gene":
            d = d.merge(stats, on="ensg", how="left")  # how="left" conserve l'ordre des lignes
        else:
            d = d.assign(mean_l2fc=mu, sd_l2fc=sd)
        d["z_l2fc"] = (d[VALUE] - d["mean_l2fc"]) / d["sd_l2fc"]
        out[name] = d
    return out


def write(sheets, path):
    try:
        import xlsxwriter  # noqa: F401  (plus rapide qu'openpyxl en écriture)
        engine = "xlsxwriter"
    except ImportError:
        engine = "openpyxl"
    with pd.ExcelWriter(path, engine=engine) as w:
        for name, d in sheets.items():
            d.to_excel(w, sheet_name=name, index=False)


def main():
    missing = [f for f in FILES if not os.path.exists(f)]
    if missing:
        sys.exit("Fichiers introuvables :\n" + "\n".join(missing))
    os.makedirs(OUT_DIR, exist_ok=True)

    stats = compute_stats(FILES) if MODE != "source" else None

    for f in FILES:
        print(f"Traitement de {os.path.basename(f)} ...")
        sheets = add_columns(read_sheets(f), stats)
        out = os.path.join(OUT_DIR, os.path.basename(f).replace(".xlsx", "_zscore.xlsx"))
        write(sheets, out)
        print(f"  -> {out}")


if __name__ == "__main__":
    main()
