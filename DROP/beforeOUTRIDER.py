#A ready-to-use Python snippet to integrate into your DROP pipeline before OUTRIDER, 
#which filters for protein-coding genes and high-expression genes, 
#ensuring numerical stability for your 96-sample cohort with Gencode v48.


# ------------------------------
# DROP pre-filter for OUTRIDER
# ------------------------------

import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

# Path to your unfitted OUTRIDER object
ods_path = "output/processed_results/aberrant_expression/v48/outrider/outrider/ods_unfitted.Rds"
ods_filtered_path = "output/processed_results/aberrant_expression/v48/outrider/outrider/ods_unfitted_filtered.Rds"

# Load R library
r('library(OUTRIDER)')

# Load ods object
r(f'ods <- readRDS("{ods_path}")')

# --- Filter high-expression genes ---
# keep only genes with mean counts > threshold (e.g., 10)
r('keep_expr <- rowMeans(assay(ods)) > 10')
r('ods <- ods[keep_expr, ]')

# --- Filter protein-coding genes only ---
# assumes rowData(ods)$gene_type is present
r('if("gene_type" %in% colnames(rowData(ods))) { '
  'keep_pc <- rowData(ods)$gene_type == "protein_coding"; '
  'ods <- ods[keep_pc, ] }')

# Recalculate size factors
r('ods <- estimateSizeFactors(ods)')

# Save filtered ods
r(f'saveRDS(ods, "{ods_filtered_path}")')

print("✅ ods object filtered and saved successfully for OUTRIDER")