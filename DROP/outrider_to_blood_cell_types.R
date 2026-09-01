#!/usr/bin/env Rscript
# ==============================================================================
# outrider_to_blood_cell_types.R
#
# Extrait l'expression normalisee (TPM) depuis un OutriderDataSet issu du
# pipeline DROP, puis lance une deconvolution cellulaire du sang total
# (quanTIseq + ABIS) via le package immunedeconv.
#
# Usage:
#   Rscript outrider_to_blood_cell_types.R <chemin_vers_ods.Rds> <dossier_sortie>
#
# Exemple:
#   Rscript outrider_to_blood_cell_types.R \
#       Output/processed_results/aberrant_expression/outrider/RNASeq_Run5/ods.Rds \
#       Fichiers_annotes/blood_cell_types
# ==============================================================================

suppressPackageStartupMessages({
  library(OUTRIDER)      # deja installe si tu fais tourner DROP
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript outrider_to_blood_cell_types.R <ods.Rds> <output_dir>")
}
ods_path   <- args[1]
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------
# 1. Charger l'OutriderDataSet et extraire les comptages bruts + longueurs
# ------------------------------------------------------------------------
message("Chargement de l'ods : ", ods_path)
ods <- readRDS(ods_path)

# Comptages bruts (genes x echantillons)
raw_counts <- counts(ods, normalized = FALSE)
message("Dimensions comptages bruts : ", paste(dim(raw_counts), collapse = " x "))

# ------------------------------------------------------------------------
# 2. Conversion en TPM
# ------------------------------------------------------------------------
# OUTRIDER/DROP stocke generalement les largeurs de genes (basepairs) dans
# les rowRanges de l'objet si l'annotation GTF a ete fournie. On les recupere
# pour un calcul TPM correct. Si absent, on bascule sur biomaRt / GTF externe.

if (!is.null(rowRanges(ods)) && "basepairs" %in% colnames(mcols(rowRanges(ods)))) {
  gene_lengths <- mcols(rowRanges(ods))$basepairs
  names(gene_lengths) <- rownames(ods)
} else {
  stop(
    "Impossible de trouver les longueurs de genes (basepairs) dans l'ods.\n",
    "Fournis un fichier de longueurs de genes separe, ou recalcule-les depuis ",
    "le GTF utilise pour l'annotation DROP (colonne 'basepairs' attendue)."
  )
}

# Garder uniquement les genes communs entre comptages et longueurs
common_genes <- intersect(rownames(raw_counts), names(gene_lengths))
raw_counts   <- raw_counts[common_genes, , drop = FALSE]
gene_lengths <- gene_lengths[common_genes]

# TPM = (counts / length_kb) puis normalisation par million
rate <- raw_counts / (gene_lengths / 1000)
tpm  <- t(t(rate) / colSums(rate)) * 1e6

# ------------------------------------------------------------------------
# 3. Conversion Ensembl ID -> gene symbol (immunedeconv attend des symbols)
# ------------------------------------------------------------------------
ensembl_ids <- sub("\\..*$", "", rownames(tpm))  # enlever suffixe de version

# Essai biomaRt avec plusieurs mirrors + retry, puis repli local sur
# org.Hs.eg.db si tous les mirrors echouent (serveur Ensembl instable/405).
get_mapping_biomart <- function(ids, mirrors = c("useast", "uswest", "asia", "www"),
                                 n_retry = 2, pause = 3) {
  suppressPackageStartupMessages(library(biomaRt))
  for (m in mirrors) {
    for (attempt in seq_len(n_retry)) {
      message(sprintf("biomaRt: essai mirror '%s' (tentative %d/%d)...", m, attempt, n_retry))
      mart <- tryCatch(
        useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = m),
        error = function(e) NULL
      )
      if (is.null(mart)) { Sys.sleep(pause); next }
      mp <- tryCatch(
        getBM(
          attributes = c("ensembl_gene_id", "hgnc_symbol"),
          filters    = "ensembl_gene_id",
          values     = ids,
          mart       = mart
        ),
        error = function(e) NULL
      )
      if (!is.null(mp)) {
        message("biomaRt OK via mirror: ", m)
        return(mp)
      }
      Sys.sleep(pause)
    }
  }
  NULL
}

get_mapping_orgdb <- function(ids) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop(
      "org.Hs.eg.db n'est pas installe. Installe-le avec :\n",
      "  BiocManager::install('org.Hs.eg.db')\n",
      "pour permettre le fallback local sans reseau."
    )
  }
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(AnnotationDbi))
  symbols <- suppressMessages(
    mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL",
           keytype = "ENSEMBL", multiVals = "first")
  )
  data.frame(
    ensembl_gene_id = names(symbols),
    hgnc_symbol     = unname(symbols),
    stringsAsFactors = FALSE
  )
}

message("Conversion Ensembl ID -> gene symbol via biomaRt (avec repli local)...")
mapping <- get_mapping_biomart(ensembl_ids)
if (is.null(mapping)) {
  message("Tous les mirrors biomaRt ont echoue. Repli sur org.Hs.eg.db (local, sans reseau)...")
  mapping <- get_mapping_orgdb(ensembl_ids)
}

mapping <- mapping[!is.na(mapping$hgnc_symbol) & mapping$hgnc_symbol != "", ]
mapping <- mapping[!duplicated(mapping$ensembl_gene_id), ]

rownames(tpm) <- ensembl_ids
tpm <- tpm[rownames(tpm) %in% mapping$ensembl_gene_id, , drop = FALSE]
symbol_map <- setNames(mapping$hgnc_symbol, mapping$ensembl_gene_id)
rownames(tpm) <- symbol_map[rownames(tpm)]

# En cas de doublons de symbols (plusieurs Ensembl ID -> meme gene), on garde
# la ligne avec la plus forte expression moyenne
if (any(duplicated(rownames(tpm)))) {
  message("Fusion des doublons de gene symbols (on garde le plus exprime)...")
  avg_expr <- rowMeans(tpm)
  keep <- !duplicated(rownames(tpm)) |
    ave(avg_expr, rownames(tpm), FUN = function(x) x == max(x))
  tpm <- tpm[order(-avg_expr), ]
  tpm <- tpm[!duplicated(rownames(tpm)), ]
}

write.csv(tpm, file.path(output_dir, "expression_tpm.csv"))
message("Matrice TPM sauvegardee : ", file.path(output_dir, "expression_tpm.csv"))

# ------------------------------------------------------------------------
# 4. Deconvolution cellulaire
# ------------------------------------------------------------------------
suppressPackageStartupMessages(library(immunedeconv))

message("Lancement quanTIseq...")
res_quantiseq <- deconvolute(tpm, method = "quantiseq")

message("Lancement ABIS...")
res_abis <- deconvolute(tpm, method = "abis")

write.csv(res_quantiseq, file.path(output_dir, "deconvolution_quantiseq.csv"), row.names = FALSE)
write.csv(res_abis,      file.path(output_dir, "deconvolution_abis.csv"),      row.names = FALSE)

message("Termine. Resultats dans : ", output_dir)
message("  - expression_tpm.csv")
message("  - deconvolution_quantiseq.csv")
message("  - deconvolution_abis.csv")
