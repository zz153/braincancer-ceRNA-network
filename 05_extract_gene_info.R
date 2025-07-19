#!/usr/bin/env Rscript

# --------------------------------------------------------
# Script: extract_gene_info.R
# Purpose: Download GENCODE v38 GTF, extract gene annotation,
#          and save as gene_info.csv for ceRNA pipeline.
# --------------------------------------------------------

# 1. Load required libraries
if (!requireNamespace("rtracklayer", quietly=TRUE)) install.packages("rtracklayer")
if (!requireNamespace("R.utils", quietly=TRUE)) install.packages("R.utils")
library(rtracklayer)
library(R.utils)

# 2. Download and unzip GTF if not already present
gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz"
gtf_file <- "gencode.v38.annotation.gtf.gz"
gtf_unzipped <- "gencode.v38.annotation.gtf"

if (!file.exists(gtf_file)) {
  cat("Downloading GTF file...\n")
  download.file(gtf_url, destfile = gtf_file, mode = "wb")
}
if (!file.exists(gtf_unzipped)) {
  cat("Unzipping GTF file...\n")
  R.utils::gunzip(gtf_file, overwrite = TRUE)
}

# 3. Import GTF and extract gene info
cat("Parsing GTF and extracting gene info...\n")
gtf <- rtracklayer::import(gtf_unzipped)
genes_gtf <- gtf[gtf$type == "gene"]

gene_info <- data.frame(
  ensembl_gene_id = gsub("\\..*", "", genes_gtf$gene_id),
  gene_biotype = genes_gtf$gene_type,
  external_gene_name = genes_gtf$gene_name,
  stringsAsFactors = FALSE
)

# 4. Save to CSV for ceRNA scripts
write.csv(gene_info, "gene_info.csv", row.names = FALSE)
cat("gene_info.csv saved in working directory!\n")
