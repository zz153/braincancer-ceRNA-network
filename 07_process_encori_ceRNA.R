#!/usr/bin/env Rscript

# --------------------------------------------------------
# Script: process_encori_ceRNA.R
# Purpose: Process ENCORI CLIP-supported lncRNA–miRNA file for ceRNA analysis
#          Output: validated_lnc_mir.csv
# --------------------------------------------------------

# 1. Set ENCORI file name (assume already downloaded)
encori_file <- "ENCORI_lncRNA_miRNA_CLIP.txt"

# 2. Load libraries
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
if (!requireNamespace("miRBaseConverter", quietly=TRUE)) install.packages("miRBaseConverter")
library(dplyr)
library(miRBaseConverter)

# 3. Detect and skip comment lines in ENCORI file
lines <- readLines(encori_file, n = 10)
skip_n <- sum(grepl("^#", lines))
cat("Auto-detected", skip_n, "comment lines to skip in ENCORI file\n")

# 4. Read ENCORI file with correct skip value
encori <- read.delim(
  encori_file,
  stringsAsFactors = FALSE,
  header = TRUE,
  skip = skip_n,
  comment.char = "#"
)

# 5. Select relevant columns
encori_pairs <- encori %>%
  dplyr::select(miRNA = miRNAname, lncRNA = geneName)

# 6. Map mature miRNA to precursor
prec_map <- as.data.frame(
  miRNA_MatureToPrecursor(unique(encori_pairs$miRNA)),
  stringsAsFactors = FALSE
)
colnames(prec_map) <- c("Mature", "Precursor")
encori_pairs <- merge(encori_pairs, prec_map, by.x = "miRNA", by.y = "Mature", all.x = TRUE)

# 7. Load gene_info data.frame (must be prepared by user)
# gene_info should have: external_gene_name, ensembl_gene_id, gene_biotype
if (!file.exists("gene_info.csv")) {
  stop("Please provide gene_info.csv with columns external_gene_name, ensembl_gene_id, and gene_biotype!")
}
gene_info <- read.csv("gene_info.csv", stringsAsFactors = FALSE)

# 8. Map lncRNA symbol to Ensembl gene ID
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

# 9. Filter for validated, mapped pairs
validated_lnc_mir <- encori_annot %>%
  filter(!is.na(Precursor) & !is.na(ensembl_gene_id)) %>%
  dplyr::select(miRNA = Precursor, lncRNA_ensembl = ensembl_gene_id) %>%
  distinct()

# 10. Output for use in ceRNA analysis
cat("Validated ENCORI lncRNA–miRNA pairs:", nrow(validated_lnc_mir), "\n")
write.csv(validated_lnc_mir, "validated_lnc_mir.csv", row.names = FALSE)
cat("Done! File saved: validated_lnc_mir.csv\n")

