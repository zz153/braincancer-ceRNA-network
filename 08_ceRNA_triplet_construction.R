#!/usr/bin/env Rscript

# -----------------------------------------------------------
# Script: 08_ceRNA_triplet_construction.R
# Purpose: Build ceRNA triplets (lncRNA–miRNA–mRNA) network 
#          using DE results and validated ENCORI & miRTarBase data
# -----------------------------------------------------------

# 1. Load Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# 2. Define Input File Paths (adjust as needed)
de_file        <- "DEG_all_Tumor_vs_AllNormals_batchCorrected.csv"
mir_file       <- "DE_miRNA_Tumor_vs_Normal.csv"
annot_file     <- "gene_info.csv"
lncmir_file    <- "validated_lnc_mir.csv"
mirmrna_file   <- "validated_mir_mrna_db.csv"

# 3. Load Data
res_annot          <- read_csv(de_file, show_col_types = FALSE)
gene_info          <- read_csv(annot_file, show_col_types = FALSE)
res_mir            <- read_csv(mir_file, show_col_types = FALSE)
validated_lnc_mir  <- read_csv(lncmir_file, show_col_types = FALSE)
validated_mir_mrna <- read_csv(mirmrna_file, show_col_types = FALSE)

# 4. Harmonize and Annotate Gene Table
res_annot <- res_annot %>%
  rename(ensembl_gene_id = ...1) %>%
  left_join(gene_info, by = "ensembl_gene_id")

# 5. Upregulated lncRNAs and mRNAs
up_lncRNAs <- res_annot %>%
  filter(gene_biotype == "lncRNA", adj.P.Val < 0.05, logFC > 1) %>%
  pull(ensembl_gene_id)

up_mRNAs <- res_annot %>%
  filter(gene_biotype == "protein_coding", adj.P.Val < 0.05, logFC > 1) %>%
  pull(ensembl_gene_id)

# 6. Subset validated interactions for upregulated nodes
validated_lnc_mir_up <- validated_lnc_mir %>%
  filter(lncRNA_ensembl %in% up_lncRNAs)

validated_mir_mrna_up <- validated_mir_mrna %>%
  filter(mRNA_ensembl %in% up_mRNAs) %>%
  rename(miRNA = miRNA_precursor)

# 7. Find common miRNAs ("bridges") between validated lncRNA–miRNA and miRNA–mRNA
common_miRNAs <- intersect(validated_lnc_mir_up$miRNA, validated_mir_mrna_up$miRNA)
cat("Number of bridge miRNAs:", length(common_miRNAs), "\n")

# 8. Identify downregulated miRNAs from DE results
res_mir <- res_mir %>% rename(miRNA = ...1)
down_mirnas <- res_mir %>%
  filter(logFC < -1, adj.P.Val < 0.05) %>%
  pull(miRNA)
cat("Number of downregulated miRNAs in DE set:", length(down_mirnas), "\n")

# 9. Find bridge miRNAs that are also downregulated
bridge_down_mirnas <- intersect(common_miRNAs, down_mirnas)
cat("Number of bridge miRNAs that are downregulated:", length(bridge_down_mirnas), "\n")
print(bridge_down_mirnas)

# 10. Subset interactions for these "bridge" miRNAs
lnc_mir_ce <- validated_lnc_mir_up %>%
  filter(miRNA %in% bridge_down_mirnas) %>%
  as.data.frame()

mir_mrna_ce <- validated_mir_mrna_up %>%
  filter(miRNA %in% bridge_down_mirnas) %>%
  as.data.frame()

# 11. Build ceRNA triplets (lncRNA–miRNA–mRNA)
ceRNA_triplets <- inner_join(
  lnc_mir_ce %>% rename(lncRNA = lncRNA_ensembl),
  mir_mrna_ce %>% rename(mRNA = mRNA_ensembl),
  by = "miRNA",
  relationship = "many-to-many"
) %>%
  select(lncRNA, miRNA, mRNA) %>%
  distinct()

cat("Number of ceRNA triplets constructed:", nrow(ceRNA_triplets), "\n")
print(head(ceRNA_triplets, 10))

# 12. Output for downstream analysis or visualization
write_csv(ceRNA_triplets, "ceRNA_triplets.csv")
cat("ceRNA triplets written to ceRNA_triplets.csv\n")

# ------------------- END OF SCRIPT ------------------- #
