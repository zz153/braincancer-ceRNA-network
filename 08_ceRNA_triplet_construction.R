#!/usr/bin/env Rscript

# -----------------------------------------------------------
# Script: 08_ceRNA_triplet_construction.R
# Purpose: Build ceRNA triplets (lncRNA–miRNA–mRNA) network 
#          using DE results and validated ENCORI & miRTarBase data
# -----------------------------------------------------------

# --------------------- #
#    LOAD DATA FILES    #
# --------------------- #
de_file     <- "DEG_all_Tumor_vs_AllNormals_batchCorrected.csv"
mir_file    <- "DE_miRNA_Tumor_vs_Normal.csv"
annot_file  <- "gene_info.csv"
lncmir_file <- "validated_lnc_mir.csv"
mirmrna_file<- "validated_mir_mrna_db.csv"

# Read files
res_annot         <- read_csv(de_file, show_col_types = FALSE)
gene_info         <- read_csv(annot_file, show_col_types = FALSE)
res_mir           <- read_csv(mir_file, show_col_types = FALSE)
validated_lnc_mir <- read_csv(lncmir_file, show_col_types = FALSE)
validated_mir_mrna <- read_csv(mirmrna_file, show_col_types = FALSE)

library(dplyr)

#Rename gene ID column
res_annot <- res_annot %>%
  rename(ensembl_gene_id = ...1)

#Join with gene_info to add biotype and gene symbol
res_annot <- res_annot %>%
  left_join(gene_info, by = "ensembl_gene_id")


#Get upregulated lncRNAs (Ensembl IDs)
up_lncRNAs <- res_annot %>%
  filter(
    gene_biotype == "lncRNA",
    adj.P.Val < 0.05,
    logFC > 1
  ) %>%
  pull(ensembl_gene_id)

#Subset validated lncRNA–miRNA pairs for those upregulated lncRNAs
validated_lnc_mir_up <- validated_lnc_mir %>%
  filter(lncRNA_ensembl %in% up_lncRNAs)

cat("Upregulated validated lncRNA–miRNA pairs:", nrow(validated_lnc_mir_up), "\n")
head(validated_lnc_mir_up)
# 3. Check the first few rows and biotype distribution
head(res_annot)
table(res_annot$gene_biotype, useNA = "ifany")

library(dplyr)

#Rename miRNA ID column for clarity
res_mir <- res_mir %>%
  rename(miRNA = ...1)

#Quick check
head(res_mir)

#Define upregulated mRNAs
up_mRNAs <- res_annot %>%
  filter(
    gene_biotype == "protein_coding",
    adj.P.Val < 0.05,
    logFC > 1
  ) %>%
  pull(ensembl_gene_id)

#Find which validated miRTarBase mRNAs are upregulated
validated_mir_mrna_up <- validated_mir_mrna %>%
  filter(mRNA_ensembl %in% up_mRNAs)

cat("Upregulated validated miRNA–mRNA pairs:", nrow(validated_mir_mrna_up), "\n")
head(validated_mir_mrna_up)

#Get upregulated lncRNAs (Ensembl IDs)
up_lncRNAs <- res_annot %>%
  filter(
    gene_biotype == "lncRNA",
    adj.P.Val < 0.05,
    logFC > 1
  ) %>%
  pull(ensembl_gene_id)

#Subset validated lncRNA–miRNA pairs for those upregulated lncRNAs
validated_lnc_mir_up <- validated_lnc_mir %>%
  filter(lncRNA_ensembl %in% up_lncRNAs)

cat("Upregulated validated lncRNA–miRNA pairs:", nrow(validated_lnc_mir_up), "\n")
head(validated_lnc_mir_up)

library(dplyr)

#Make sure the miRNA column is named the same in both dataframes
validated_mir_mrna_up <- validated_mir_mrna_up %>%
  rename(miRNA = miRNA_precursor)

#Find common miRNAs between lncRNA–miRNA and miRNA–mRNA
common_miRNAs <- intersect(validated_lnc_mir_up$miRNA, validated_mir_mrna_up$miRNA)
cat("Number of bridge miRNAs:", length(common_miRNAs), "\n")
print(head(common_miRNAs))

down_mirnas <- res_mir %>%
  filter(logFC < -1, adj.P.Val < 0.05) %>%
  pull(miRNA)

cat("Number of downregulated miRNAs in DE set:", length(down_mirnas), "\n")

# Which common miRNAs are downregulated?
bridge_down_mirnas <- intersect(common_miRNAs, down_mirnas)

cat("Number of bridge miRNAs that are downregulated:", length(bridge_down_mirnas), "\n")
print(bridge_down_mirnas)

head(validated_lnc_mir_up)

# Assuming validated_lnc_mir_up and validated_mrna_mir_up both have a 'miRNA' column

# lncRNA–miRNA pairs for downregulated bridge miRNAs
lnc_mir_ce <- validated_lnc_mir_up %>%
  filter(miRNA %in% bridge_down_mirnas)

# miRNA–mRNA pairs for downregulated bridge miRNAs
mir_mrna_ce <- validated_mir_mrna_up %>%
  filter(miRNA %in% bridge_down_mirnas)

lnc_mir_ce <- as.data.frame(lnc_mir_ce)
mir_mrna_ce <- as.data.frame(mir_mrna_ce)

ceRNA_triplets <- inner_join(
  lnc_mir_ce %>% rename(lncRNA = lncRNA_ensembl),
  mir_mrna_ce %>% rename(mRNA = mRNA_ensembl),
  by = "miRNA",
  relationship = "many-to-many"
) %>%
  dplyr::select(lncRNA, miRNA, mRNA) %>%
  distinct()


cat("Unique lncRNAs: ", length(unique(ceRNA_triplets$lncRNA)), "\n")
cat("Unique miRNAs:  ", length(unique(ceRNA_triplets$miRNA)), "\n")
cat("Unique mRNAs:   ", length(unique(ceRNA_triplets$mRNA)), "\n")

unique_lncRNAs <- unique(ceRNA_triplets$lncRNA)
unique_miRNAs  <- unique(ceRNA_triplets$miRNA)
unique_mRNAs   <- unique(ceRNA_triplets$mRNA)

# To view them:
head(unique_lncRNAs)
head(unique_miRNAs)
head(unique_mRNAs)

# Save as one-column CSVs
write.csv(data.frame(lncRNA = unique_lncRNAs), "unique_lncRNAs.csv", row.names = FALSE)
write.csv(data.frame(miRNA  = unique_miRNAs),  "unique_miRNAs.csv",  row.names = FALSE)
write.csv(data.frame(mRNA   = unique_mRNAs),   "unique_mRNAs.csv",   row.names = FALSE)







