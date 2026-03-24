#!/usr/bin/env python3

import pandas as pd
import sys
import os
from itertools import combinations

def detect_gene_column(df):
    possible_cols = ["geneID", "hgncSymbol", "gene_id", "Gene", "ENSEMBL"]
    for col in possible_cols:
        if col in df.columns:
            return col
    return None

def load_genes(file):
    df = pd.read_csv(file, sep="\t")
    gene_col = detect_gene_column(df)
    
    if gene_col is None:
        raise ValueError(f"Aucune colonne gène détectée dans {file}")
    
    genes = set(df[gene_col].dropna().unique())
    return genes, gene_col

def main(files):
    gene_sets = {}
    
    for file in files:
        genes, gene_col = load_genes(file)
        gene_sets[file] = genes
        print(f"\n📄 {file}")
        print(f"   Colonne gène : {gene_col}")
        print(f"   Nombre de gènes uniques : {len(genes)}")

    print("\n🔬 Comparaisons pair-à-pair :")
    for f1, f2 in combinations(files, 2):
        common = gene_sets[f1].intersection(gene_sets[f2])
        only_f1 = gene_sets[f1] - gene_sets[f2]
        only_f2 = gene_sets[f2] - gene_sets[f1]

        print(f"\n🆚 {os.path.basename(f1)} vs {os.path.basename(f2)}")
        print(f"   Gènes communs : {len(common)}")
        print(f"   Spécifiques à {os.path.basename(f1)} : {len(only_f1)}")
        print(f"   Spécifiques à {os.path.basename(f2)} : {len(only_f2)}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python compare_genes.py fichier1.tsv fichier2.tsv [fichier3.tsv ...]")
        sys.exit(1)

    main(sys.argv[1:])