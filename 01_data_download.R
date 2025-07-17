# scripts/01_data_download.R
# ---------------------------------------------------
# Download TCGA-GBM and GTEx brain gene expression data
# plus GENCODE annotation GTF for gene info.
# ---------------------------------------------------

# Load required packages (assumes Bioconductor for TCGAbiolinks, recount3, etc.)
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("TCGAbiolinks")
}
if (!requireNamespace("recount3", quietly = TRUE)) BiocManager::install("recount3")
if (!requireNamespace("rtracklayer", quietly = TRUE)) BiocManager::install("rtracklayer")
if (!requireNamespace("R.utils", quietly = TRUE)) install.packages("R.utils")

library(TCGAbiolinks)
library(recount3)
library(rtracklayer)
library(R.utils)

# --------------------------------------
# 1. Download TCGA GBM gene expression
# --------------------------------------
query_tcga <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
GDCdownload(query_tcga)
saveRDS(query_tcga, file = "TCGA_GBM_query.rds")  # Save query for reproducibility

tcga_data <- GDCprepare(query_tcga)
saveRDS(tcga_data, file = "TCGA_GBM_data.rds")    # Save SummarizedExperiment object

# --------------------------------------
# 2. Download TCGA GBM miRNA expression
# --------------------------------------
query_mir <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "miRNA Expression Quantification",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
GDCdownload(query_mir)
saveRDS(query_mir, file = "TCGA_GBM_miRNA_query.rds")

mir_data <- GDCprepare(query_mir)
saveRDS(mir_data, file = "TCGA_GBM_miRNA_data.rds")

# --------------------------------------
# 3. Download GTEx BRAIN expression data (via recount3)
# --------------------------------------
gtex_proj <- subset(available_projects("human"), project == "BRAIN" & file_source == "gtex")
rse_brain <- create_rse(gtex_proj)
saveRDS(rse_brain, file = "GTEx_Brain_rse.rds")

# --------------------------------------
# 4. Download GENCODE v38 annotation GTF (for gene types/info)
# --------------------------------------
gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"
gtf_file <- "gencode.v38.annotation.gtf.gz"
gtf_unzipped <- "gencode.v38.annotation.gtf"
options(timeout = 300)
if (!file.exists(gtf_file)) {
  tryCatch({
    download.file(gtf_url, destfile = gtf_file, method = "curl")
  }, error = function(e) {
    download.file(gtf_url, destfile = gtf_file)
  })
}
if (!file.exists(gtf_unzipped)) {
  R.utils::gunzip(gtf_file, overwrite = TRUE)
}
# You can now load this GTF in downstream steps

cat("All downloads complete!\n")
