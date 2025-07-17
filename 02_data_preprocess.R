# scripts/02_data_preprocess.R
# ---------------------------------------------------
# Load, clean, and harmonize TCGA/GTEx data for downstream ceRNA analysis
# ---------------------------------------------------

# 1. Load libraries
source("scripts/00_libraries.R")

# 2. Load downloaded RDS and GTF files
tcga_data   <- readRDS("TCGA_GBM_data.rds")
mir_data    <- readRDS("TCGA_GBM_miRNA_data.rds")
rse_brain   <- readRDS("GTEx_Brain_rse.rds")
gtf         <- rtracklayer::import("gencode.v38.annotation.gtf")

# 3. Extract count matrices, clean gene IDs, etc.
# (Insert your pipeline’s code for normalization, filtering, harmonization...)

# 4. Save processed data for downstream scripts
saveRDS(processed_tcga, "processed_TCGA.rds")
saveRDS(processed_gtex, "processed_GTEx.rds")
# etc.

cat("Data preprocessing complete!\n")
