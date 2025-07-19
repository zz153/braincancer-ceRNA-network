#!/usr/bin/env Rscript

# --------------------------------------------------------
# Script: process_encori_ceRNA.R
# Purpose: Process ENCORI CLIP-supported lncRNA–miRNA file for ceRNA analysis
#          Outputs: validated_lnc_mir.csv
# --------------------------------------------------------

# 1. (Optional) Download the ENCORI file if not present
encori_url <- "https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=lncRNA&miRNA=all&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=0&program=None&target=all&cellType=all"
encori_file <- "ENCORI_lncRNA_miRNA_CLIP.txt"
if (!file.exists(encori_file)) {
  download.file(encori_url, destfile = encori_file, mode = "wb")
}

# 2. Load libraries
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
if (!requireNamespace("miRBaseConverter", quietly=TRUE)) install.packages("miRBaseConverter")
library(dplyr)
library(miRBaseConverter)

# 3. Read ENCORI file (skip comment lines if any)
# (If you know exactly how many lines to skip, set skip=3 or another value. Otherwise, try skip=0.)
# If error: try increasing skip to match your header row.
encori <- read.delim(encori_file, stringsAsFactors = FALSE, header = TRUE)

# 4. Select relevant columns
encori_pairs <- encori %>%
  dplyr::select(miRNA = miRNAname, lncRNA = geneName)

# 5. Map mature miRNA to precursor
prec_map <- as.data.frame(
  miRNA_MatureToPrecursor(unique(encori_pairs$miRNA)),
  stringsAsFactors = FALSE
)
colnames(prec_map) <- c("Mature", "Precursor")
encori_pairs <- merge(encori_pairs, prec_map, by.x = "miRNA", by.y = "Mature", all.x = TRUE)

# 6. Load gene_info data.frame (must be prepared by user)
# gene_info should have: external_gene_name, ensembl_gene_id, gene_biotype
if (!file.exists("gene_info.csv")) {
  stop("Please provide gene_info.csv with columns external_gene_name, ensembl_gene_id, and gene_biotype!")
}
gene_info <- read.csv("gene_info.csv", stringsAsFactors = FALSE)

# 7. Map lncRNA symbol to Ensembl gene ID
lnc_lookup <- gene_info %>%
  filter(gene_biotype == "lncRNA") %>%
  dplyr::select(external_gene_name, ensembl_gene_id)
encori_annot <- merge(
  encori_pairs,
  lnc_lookup,
  by.x = "lncRNA",
  by.y = "external_gene_name",
  all.x = TRUE
)

# 8. Filter for validated, mapped pairs
validated_lnc_mir <- encori_annot %>%
  filter(!is.na(Precursor) & !is.na(ensembl_gene_id)) %>%
  dplyr::select(miRNA = Precursor, lncRNA_ensembl = ensembl_gene_id) %>%
  distinct()

# 9. Output for use in ceRNA analysis
cat("Validated ENCORI lncRNA–miRNA pairs:", nrow(validated_lnc_mir), "\n")
write.csv(validated_lnc_mir, "validated_lnc_mir.csv", row.names = FALSE)

cat("Done! File saved: validated_lnc_mir.csv\n")
