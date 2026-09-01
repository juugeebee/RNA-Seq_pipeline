#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================
Automatic clustering of cellular deconvolution profiles
==============================================================

Compatible with:
    - CIBERSORTx output (samples in rows / "Mixture" column as
      index, cell types in columns, plus P-value / Correlation /
      RMSE stat columns which are automatically dropped)
    - immunedeconv-style output (cell types in rows, samples in
      columns, e.g. deconvolution_abis.csv / deconvolution_quantiseq.csv)

The script auto-detects the orientation and the presence of
CIBERSORTx stat columns, so no manual transposition is needed.

Input:
    deconvolution_abis.csv (or any deconvolution csv, see -i)

Output folder:
    clustering/

Generated files:
    clusters.csv
    silhouette_scores.csv
    PCA_clusters.png
    Dendrogram.png
    Heatmap.png
    CellType_Mean_ByCluster.csv
    Summary.txt

Author:
==============================================================
"""

import os
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")

###############################################################
# PARAMETERS
###############################################################

DEFAULT_INPUT_FILE = "deconvolution_abis.csv"

OUTPUT_DIR = "clustering"

MIN_CLUSTERS = 2
MAX_CLUSTERS = 10

FIG_DPI = 300

# Columns appended by CIBERSORTx that are NOT cell types and must
# be excluded from the clustering (case-insensitive match).
CIBERSORTX_STAT_COLUMNS = {"p-value", "correlation", "rmse", "pvalue"}


###############################################################
# ARGUMENTS
###############################################################

parser = argparse.ArgumentParser(
    description="Automatic clustering of cellular deconvolution profiles "
                "(CIBERSORTx or immunedeconv output)."
)
parser.add_argument(
    "-i", "--input",
    default=DEFAULT_INPUT_FILE,
    help=f"Input CSV file (default: {DEFAULT_INPUT_FILE})"
)
parser.add_argument(
    "-o", "--output",
    default=OUTPUT_DIR,
    help=f"Output directory (default: {OUTPUT_DIR})"
)
args = parser.parse_args()

INPUT_FILE = args.input
OUTPUT_DIR = args.output

###############################################################
# CREATE OUTPUT DIRECTORY
###############################################################

os.makedirs(OUTPUT_DIR, exist_ok=True)

###############################################################
# LOAD DATA (auto-detect orientation: CIBERSORTx vs immunedeconv)
###############################################################

print("Loading data...")

raw = pd.read_csv(INPUT_FILE, index_col=0)

stat_cols_present = [c for c in raw.columns if c.strip().lower() in CIBERSORTX_STAT_COLUMNS]

if stat_cols_present:
    # CIBERSORTx format: samples already in rows ("Mixture" as index),
    # cell types in columns, plus stat columns to drop.
    print(f"CIBERSORTx format detected (dropping columns: {stat_cols_present})")
    X = raw.drop(columns=stat_cols_present)
else:
    # immunedeconv format (e.g. ABIS/quanTIseq via immunedeconv):
    # cell types in rows, samples in columns -> transpose.
    print("immunedeconv-style format detected (cell types in rows), transposing...")
    X = raw.T

# Make sure everything is numeric (CIBERSORTx sometimes writes numbers as strings)
X = X.apply(pd.to_numeric, errors="coerce")

# Drop any cell type / sample that ended up entirely NaN (safety net)
X = X.dropna(axis=1, how="all").dropna(axis=0, how="all")

print(f"{X.shape[0]} samples")
print(f"{X.shape[1]} cell types")

###############################################################
# STANDARDIZATION
###############################################################

print("Standardizing data...")

scaler = StandardScaler()

X_scaled = scaler.fit_transform(X)

###############################################################
# HIERARCHICAL CLUSTERING
###############################################################

print("Running hierarchical clustering...")

Z = linkage(X_scaled, method="ward")

###############################################################
# AUTOMATIC NUMBER OF CLUSTERS
###############################################################

print("Searching optimal number of clusters...")

scores = {}

max_k = min(MAX_CLUSTERS, X.shape[0] - 1)

for k in range(MIN_CLUSTERS, max_k + 1):

    labels = fcluster(Z, k, criterion="maxclust")

    score = silhouette_score(X_scaled, labels)

    scores[k] = score

best_k = max(scores, key=scores.get)

print(f"Optimal number of clusters : {best_k}")
print(f"Silhouette score : {scores[best_k]:.3f}")

###############################################################
# FINAL CLUSTERS
###############################################################

clusters = fcluster(Z, best_k, criterion="maxclust")

cluster_table = pd.DataFrame({
    "Sample": X.index,
    "Cluster": clusters
})

cluster_table.to_csv(
    os.path.join(OUTPUT_DIR, "clusters.csv"),
    index=False
)

###############################################################
# SILHOUETTE PLOT
###############################################################

score_df = pd.DataFrame({
    "Clusters": list(scores.keys()),
    "Silhouette": list(scores.values())
})

score_df.to_csv(
    os.path.join(OUTPUT_DIR, "silhouette_scores.csv"),
    index=False
)

plt.figure(figsize=(7,5))

plt.plot(
    score_df["Clusters"],
    score_df["Silhouette"],
    "-o",
    linewidth=2,
    markersize=8
)

plt.grid(alpha=0.3)

plt.title("Automatic Selection of the Number of Clusters")

plt.xlabel("Number of clusters")

plt.ylabel("Mean silhouette score")

plt.tight_layout()

plt.savefig(
    os.path.join(OUTPUT_DIR, "Silhouette_scores.png"),
    dpi=FIG_DPI
)

plt.close()

###############################################################
# PCA
###############################################################

print("Generating PCA...")

pca = PCA(n_components=2)

coord = pca.fit_transform(X_scaled)

plt.figure(figsize=(9,7))

scatter = plt.scatter(
    coord[:,0],
    coord[:,1],
    c=clusters,
    cmap="tab10",
    edgecolor="black",
    s=90
)

for i, sample in enumerate(X.index):

    plt.text(
        coord[i,0],
        coord[i,1],
        sample,
        fontsize=6
    )

plt.xlabel(
    f"Principal Component 1 ({100*pca.explained_variance_ratio_[0]:.1f}% variance explained)"
)

plt.ylabel(
    f"Principal Component 2 ({100*pca.explained_variance_ratio_[1]:.1f}% variance explained)"
)

plt.title("Principal Component Analysis of Cellular Profiles")

legend = plt.legend(
    *scatter.legend_elements(),
    title="Cluster",
    loc="best"
)

plt.gca().add_artist(legend)

plt.grid(alpha=0.25)

plt.tight_layout()

plt.savefig(
    os.path.join(OUTPUT_DIR, "PCA_clusters.png"),
    dpi=FIG_DPI
)

plt.close()

###############################################################
# DENDROGRAM
###############################################################

print("Generating dendrogram...")

plt.figure(figsize=(15,7))

dendrogram(
    Z,
    labels=X.index,
    leaf_rotation=90,
    leaf_font_size=7
)

plt.title("Hierarchical Clustering of Samples")

plt.xlabel("Samples")

plt.ylabel("Ward distance")

plt.tight_layout()

plt.savefig(
    os.path.join(OUTPUT_DIR, "Dendrogram.png"),
    dpi=FIG_DPI
)

plt.close()

###############################################################
# HEATMAP
###############################################################

print("Generating heatmap...")

ordered_samples = cluster_table.sort_values("Cluster")["Sample"]

ordered = X.loc[ordered_samples]

cluster_colors = pd.Series(clusters, index=X.index).loc[ordered_samples]

palette = sns.color_palette("tab10", best_k)

lut = {
    i + 1: palette[i]
    for i in range(best_k)
}

col_colors = cluster_colors.map(lut)

g = sns.clustermap(
    ordered.T,
    cmap="viridis",
    figsize=(16,9),
    col_cluster=False,
    row_cluster=True,
    col_colors=col_colors,
    linewidths=0.2,
    cbar_kws={
        "label": "Standardized cell proportion"
    }
)

g.fig.suptitle(
    "Heatmap of Cellular Composition",
    fontsize=16,
    y=1.02
)

g.savefig(
    os.path.join(OUTPUT_DIR, "Heatmap.png"),
    dpi=FIG_DPI
)

plt.close()

###############################################################
# MEAN CELL TYPE PER CLUSTER
###############################################################

print("Calculating cluster averages...")

means = ordered.copy()

means["Cluster"] = cluster_colors.values

cluster_mean = means.groupby("Cluster").mean()

cluster_mean.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "CellType_Mean_ByCluster.csv"
    )
)

###############################################################
# SUMMARY REPORT
###############################################################

print("Writing report...")

with open(
    os.path.join(OUTPUT_DIR, "Summary.txt"),
    "w"
) as f:

    f.write("=====================================================\n")
    f.write("Automatic clustering report\n")
    f.write("=====================================================\n\n")

    f.write(f"Input file : {INPUT_FILE}\n")
    f.write(f"Number of samples : {X.shape[0]}\n")
    f.write(f"Number of cell types : {X.shape[1]}\n\n")

    f.write(f"Optimal number of clusters : {best_k}\n")
    f.write(f"Silhouette score : {scores[best_k]:.3f}\n\n")

    f.write("Samples per cluster\n")
    f.write("-------------------\n")

    for c in sorted(np.unique(clusters)):
        n = np.sum(clusters == c)
        f.write(f"Cluster {c}: {n} samples\n")

    f.write("\n")

    f.write("Average cell proportions by cluster saved in:\n")
    f.write("CellType_Mean_ByCluster.csv\n")

###############################################################
# FINISHED
###############################################################

print("\n========================================")
print("Analysis completed successfully")
print("========================================")
print(f"Results saved in : {OUTPUT_DIR}")
print(f"Optimal number of clusters : {best_k}")
print(f"Silhouette score : {scores[best_k]:.3f}")
print("========================================")
