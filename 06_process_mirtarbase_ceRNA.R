#!/usr/bin/env Rscript

# --------------------------------------------------------
# Script: process_mirtarbase_ceRNA.R
# Purpose: Process miRTarBase 2025 file for ceRNA analysis
#          Outputs: validated_mir_mrna_db.csv
# --------------------------------------------------------

# 1. (Optional) Download the miRTarBase file if not present
if (!file.exists("miRTarBase_MTI_2025.xlsx")) {
  download.file(
    "https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/data/miRTarBase_MTI.xlsx",
    destfile = "miRTarBase_MTI_2025.xlsx",
    mode = "wb"
  )
}

# 2. Load libraries
if (!requireNamespace("readxl", quietly=TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
if (!requireNamespace("miRBaseConverter", quietly=TRUE)) install.packages("miRBaseConverter")
library(readxl)
library(dplyr)
library(miRBaseConverter)

# 3. Read in miRTarBase Excel file
mirtar <- readxl::read_xlsx("miRTarBase_MTI_2025.xlsx")

# 4. Subset for human miRNA-mRNA interactions
mirtar_hsa <- mirtar %>%
  filter(
    `Species (miRNA)` == "Homo sapiens",
    `Species (Target Gene)` == "Homo sapiens"
  ) %>%
  dplyr::select(miRNA, `Target Gene`)

# 5. Clean whitespace and keep mature human miRNAs only
mirtar_hsa$miRNA <- trimws(as.character(mirtar_hsa$miRNA))
mirtar_hsa$`Target Gene` <- trimws(as.character(mirtar_hsa$`Target Gene`))
mirtar_hsa <- mirtar_hsa %>% filter(grepl("^hsa-", miRNA), !is.na(miRNA))

# 6. Map mature miRNA → precursor
prec_map <- as.data.frame(
  miRNA_MatureToPrecursor(unique(mirtar_hsa$miRNA)),
  stringsAsFactors = FALSE
)
colnames(prec_map) <- c("Mature", "Precursor")
mirtar_hsa <- merge(mirtar_hsa, prec_map, by.x = "miRNA", by.y = "Mature", all.x = TRUE)

# 7. Load gene_info data.frame (must be prepared by user)
# gene_info should have: external_gene_name, ensembl_gene_id, gene_biotype
if (!file.exists("gene_info.csv")) {
  stop("Please provide gene_info.csv with columns external_gene_name and ensembl_gene_id!")
}
gene_info <- read.csv("gene_info.csv", stringsAsFactors = FALSE)

# 8. Map gene symbols to Ensembl gene IDs
mirtar_annot <- merge(
  mirtar_hsa,
  gene_info[, c("external_gene_name", "ensembl_gene_id")],
  by.x = "Target Gene",
  by.y = "external_gene_name",
  all.x = TRUE
)

# 9. Keep only pairs with both miRNA precursor and Ensembl gene ID
validated_mir_mrna_db <- mirtar_annot %>%
  filter(!is.na(Precursor) & !is.na(ensembl_gene_id)) %>%
  dplyr::select(miRNA_precursor = Precursor, mRNA_ensembl = ensembl_gene_id) %>%
  distinct()

# 10. Output for use in ceRNA analysis
cat("Validated miRNA–mRNA pairs (for ceRNA):", nrow(validated_mir_mrna_db), "\n")
write.csv(validated_mir_mrna_db, "validated_mir_mrna_db.csv", row.names = FALSE)

cat("Done! File saved: validated_mir_mrna_db.csv\n")
