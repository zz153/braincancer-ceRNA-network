#!/usr/bin/env Rscript

# -----------------------------------------------------------
# Script: 10_unique_lncRNAs_survivalplots.R
# Purpose: Build ceRNA triplets (lncRNA–miRNA–mRNA) network 
#          using DE results and validated ENCORI & miRTarBase data
# -----------------------------------------------------------

# --- Required Libraries (assume already loaded in your env) ---
library(survival)
library(survminer)
library(dplyr)
library(tibble)

# --- 1. Extract clinical survival info from TCGA data ---
clin_df <- as.data.frame(colData(tcga_data))
# Use "Primary Tumor" samples only for survival analysis
clin_df <- clin_df[clin_df$sample_type == "Primary Tumor", ]

# Required columns: patient barcode, vital_status, days_to_death, days_to_last_follow_up
clin_df$deceased <- clin_df$vital_status == "Dead"
clin_df$overall_survival <- ifelse(
  clin_df$deceased,
  clin_df$days_to_death,
  clin_df$days_to_last_follow_up
)

# --- 2. Get the normalized expression matrix ---
# Use voom-normalized matrix from previous step (v$E)
expr_mat <- v$E # rows: genes, columns: samples
# Subset for "Primary Tumor" columns to match clin_df
tumor_samples <- clin_df$barcode
expr_mat_tumor <- expr_mat[, tumor_samples, drop = FALSE]

# --- 3. Get gene symbols for titles ---
gene_info_tbl <- as_tibble(gene_info)
lnc_annot <- gene_info %>%
  dplyr::filter(ensembl_gene_id %in% unique_lncRNAs) %>%
  dplyr::select(ensembl_gene_id, external_gene_name)

# --- 4. Output folder for plots ---
dir.create("survival_plots", showWarnings = FALSE)

# --- 5. Loop through each lncRNA ---
for(i in seq_along(unique_lncRNAs)) {
  lnc_id <- unique_lncRNAs[i]
  # Expression for this gene
  expr <- expr_mat_tumor[lnc_id, ]
  # Merge with clinical for the right patients
  df <- clin_df
  df$lnc_expr <- expr[df$barcode]
  
  # Group by median expression (High/Low)
  median_expr <- median(df$lnc_expr, na.rm = TRUE)
  df$expr_group <- ifelse(df$lnc_expr >= median_expr, "High", "Low")
  
  # Get gene symbol for the plot title
  gene_symbol <- lnc_annot$external_gene_name[lnc_annot$ensembl_gene_id == lnc_id]
  if(length(gene_symbol) == 0 || is.na(gene_symbol)) gene_symbol <- lnc_id
  
  # --- Kaplan–Meier survival analysis ---
  surv_obj <- Surv(time = df$overall_survival, event = df$deceased)
  fit <- survfit(surv_obj ~ expr_group, data = df)
  
  # --- Plot and save ---
  surv_plot <- ggsurvplot(
    fit, data = df, pval = TRUE, risk.table = TRUE,
    ggtheme = theme_bw(base_size = 14),
    legend.title = "",
    legend.labs = c("High Expression", "Low Expression"),
    title = paste0("Survival: ", gene_symbol, " (", lnc_id, ")")
  )
  ggsave(
    filename = file.path("survival_plots", paste0("surv_", gene_symbol, "_", lnc_id, ".pdf")),
    plot = print(surv_plot), width = 6, height = 6
  )
}

cat("Kaplan–Meier survival plots saved in /survival_plots!\n")
