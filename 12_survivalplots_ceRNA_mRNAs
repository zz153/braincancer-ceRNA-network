# -----------------------------------------------------------
# Script: survival_ceRNA_mRNAs.R
# Purpose: Kaplan–Meier survival plots for all mRNAs in ceRNA
#          networks for CYTOR and MIR4435-2HG in TCGA GBM
# -----------------------------------------------------------

library(survival)
library(survminer)
library(dplyr)
library(tibble)

# --- 1. Prepare Clinical Data (Primary Tumors Only) ---
clin_df <- as.data.frame(colData(tcga_data))
clin_df <- clin_df[clin_df$sample_type == "Primary Tumor", ]
clin_df$deceased <- clin_df$vital_status == "Dead"
clin_df$overall_survival <- ifelse(
  clin_df$deceased,
  clin_df$days_to_death,
  clin_df$days_to_last_follow_up
)

# --- 2. Prepare Expression Data ---
expr_mat <- v$E # voom-normalized expression matrix
tumor_samples <- clin_df$barcode
expr_mat_tumor <- expr_mat[, tumor_samples, drop = FALSE]

# --- 3. Prepare mRNA List ---
# (Example: use your actual Ensembl ID vectors)
mrnas_cytor <- unique(ceRNA_CYTOR$mRNA)
mrnas_mir4435hg <- unique(ceRNA_MIR44352HG$mRNA)
mrna_ceRNAs <- unique(c(mrnas_cytor, mrnas_mir4435hg))

# --- 4. Annotate Gene Symbols ---
mrna_annot <- gene_info %>%
  dplyr::filter(ensembl_gene_id %in% mrna_ceRNAs) %>%
  dplyr::select(ensembl_gene_id, external_gene_name)


# --- 5. Output Folder ---
dir.create("survival_plots_mRNAs", showWarnings = FALSE)

# --- 6. Survival Analysis Loop ---
for (i in seq_along(mrna_ceRNAs)) {
  mrna_id <- mrna_ceRNAs[i]
  # Expression for this gene
  expr <- expr_mat_tumor[mrna_id, ]
  df <- clin_df
  df$gene_expr <- expr[df$barcode]
  
  # Median split
  median_expr <- median(df$gene_expr, na.rm = TRUE)
  df$expr_group <- ifelse(df$gene_expr >= median_expr, "High", "Low")
  
  # Get gene symbol for title
  gene_symbol <- mrna_annot$external_gene_name[mrna_annot$ensembl_gene_id == mrna_id]
  if (length(gene_symbol) == 0 || is.na(gene_symbol)) gene_symbol <- mrna_id
  
  # Kaplan–Meier
  surv_obj <- Surv(time = df$overall_survival, event = df$deceased)
  fit <- survfit(surv_obj ~ expr_group, data = df)
  
  surv_plot <- ggsurvplot(
    fit, data = df, pval = TRUE, risk.table = TRUE,
    ggtheme = theme_bw(base_size = 14),
    legend.title = "",
    legend.labs = c("High Expression", "Low Expression"),
    title = paste0("Survival: ", gene_symbol, " (", mrna_id, ")")
  )
  ggsave(
    filename = file.path("survival_plots_mRNAs", paste0("surv_", gene_symbol, "_", mrna_id, ".pdf")),
    plot = print(surv_plot), width = 6, height = 6
  )
}

cat("Kaplan–Meier survival plots for ceRNA mRNAs saved in /survival_plots_mRNAs!\n")

cytor_id <- "ENSG00000222041"
mir4435hg_id <- "ENSG00000172965"

mRNAs_cytor <- ceRNA_triplets %>% filter(lncRNA == cytor_id) %>% pull(mRNA)
mRNAs_mir4435hg <- ceRNA_triplets %>% filter(lncRNA == mir4435hg_id) %>% pull(mRNA)

shared_mRNAs <- intersect(mRNAs_cytor, mRNAs_mir4435hg)
length(shared_mRNAs)  # Should be 25

other_lncRNAs <- ceRNA_triplets %>%
  filter(mRNA %in% shared_mRNAs & !(lncRNA %in% c(cytor_id, mir4435hg_id))) %>%
  distinct(lncRNA, mRNA)

# Get counts: For each shared mRNA, how many lncRNAs (besides CYTOR & MIR4435-2HG) regulate it?
lncRNA_counts <- other_lncRNAs %>%
  group_by(mRNA) %>%
  summarize(num_other_lncRNAs = n(), lncRNA_list = paste(unique(lncRNA), collapse = "; "))

print(lncRNA_counts)

lncRNA_counts <- other_lncRNAs %>%
  dplyr::group_by(mRNA) %>%
  dplyr::summarise(
    num_other_lncRNAs = n(),
    lncRNA_list = paste(unique(lncRNA), collapse = "; ")
  ) %>%
  dplyr::ungroup()

count(lncRNA_counts)
