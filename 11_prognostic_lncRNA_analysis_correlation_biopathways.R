head(ceRNA_triplets)

# Using Ensembl ID for CYTOR/MIR4435-2HG
ceRNA_CYTOR <- ceRNA_triplets %>% filter(lncRNA == "ENSG00000222041")
ceRNA_MIR44352HG <- ceRNA_triplets %>% filter(lncRNA == "ENSG00000172965")

head(ceRNA_MIR44352HG)
head(ceRNA_CYTOR)

# Get unique mRNAs (Ensembl IDs)
mRNAs_CYTOR <- unique(ceRNA_CYTOR$mRNA)
mRNAs_mir44352HG <- unique(ceRNA_MIR44352HG$mRNA)

lnc_ids <- c("ENSG00000222041", "ENSG00000172965")

library(dplyr)
# This makes sure you're using dplyr's select:
all_partners <- ceRNA_triplets %>%
  filter(lncRNA %in% c("ENSG00000222041", "ENSG00000172965")) %>%
  dplyr::select(lncRNA, miRNA, mRNA) %>%
  distinct() %>%
  arrange(lncRNA, miRNA, mRNA)

expr_cytor     <- as.numeric(expr_mat["ENSG00000222041", ])
expr_mir4435hg <- as.numeric(expr_mat["ENSG00000172965", ])

# Pearson correlation (linear)
cor_test <- cor.test(expr_cytor, expr_mir4435hg, method = "pearson")
print(cor_test)

# For non-normal data, Spearman can be used:
# cor.test(expr_cytor, expr_mir4435hg, method = "spearman")

library(ggplot2)
df <- data.frame(CYTOR = expr_cytor, MIR4435_1HG = expr_mir4435hg)
ggplot(df, aes(x = CYTOR, y = MIR4435_1HG)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, col = "red") +
  labs(
    x = "CYTOR (ENSG00000222041) Expression",
    y = "MIR4435-1HG (ENSG00000172965) Expression",
    title = "Correlation of lncRNA Expression (GBM Tumor Samples)"
  ) +
  theme_minimal(base_size = 14)

# Assuming ceRNA_triplets contains your filtered results
shared_mRNAs <- ceRNA_triplets %>%
  filter(lncRNA %in% c("ENSG00000222041", "ENSG00000172965")) %>%
  pull(mRNA) %>%
  unique()
length(shared_mRNAs)
head(shared_mRNAs)

library(clusterProfiler)
library(org.Hs.eg.db)

# Convert Ensembl IDs to Entrez IDs
entrez_ids <- mapIds(org.Hs.eg.db, keys=shared_mRNAs, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
entrez_ids <- na.omit(entrez_ids)

# KEGG pathway enrichment
kegg_res <- enrichKEGG(gene=entrez_ids, organism="hsa", pvalueCutoff=0.05)
head(kegg_res)

mrna_annot <- gene_info %>%
  filter(ensembl_gene_id %in% shared_mRNAs, gene_biotype == "protein_coding")
gene_symbols <- mrna_annot$external_gene_name
print(gene_symbols)

# GO Biological Process
go_bp <- enrichGO(gene=entrez_ids, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)
head(go_bp)
head(gene_info)
library(clusterProfiler)
gene_symbols <- mRNA_annot$external_gene_name
gene_desc <- bitr(gene_symbols, fromType="SYMBOL", toType=c("GENENAME"), OrgDb="org.Hs.eg.db")
head(gene_desc)

# Set Ensembl IDs
cytor_id <- "ENSG00000222041"
mir4435hg_id <- "ENSG00000172965"

# --------- Correlation Plot: CYTOR vs MIR4435-2HG ---------
# Required: expr_mat (normalized matrix, Ensembl IDs as rownames)

library(ggplot2)
library(tibble)

# 1. Set Ensembl IDs
cytor_id <- "ENSG00000222041"
mir4435hg_id <- "ENSG00000172965"

# 2. Extract expression values (all samples or tumor subset as needed)
expr_cytor <- expr_mat[cytor_id, ]
expr_mir4435hg <- expr_mat[mir4435hg_id, ]

# 3. Build a data frame
df_corr <- tibble(
  CYTOR = as.numeric(expr_cytor),
  MIR4435_2HG = as.numeric(expr_mir4435hg)
)

# 4. Pearson correlation test
cor_test <- cor.test(df_corr$CYTOR, df_corr$MIR4435_2HG, method = "pearson")
cor_val <- round(cor_test$estimate, 3)
p_lab <- ifelse(cor_test$p.value < 2e-16, "p < 2e-16", paste0("p = ", signif(cor_test$p.value, 3)))

# 5. Plot with ggplot2
p <- ggplot(df_corr, aes(x = CYTOR, y = MIR4435_2HG)) +
  geom_point(color = "#3182bd", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "#e6550d", fill = "#e6550d", alpha = 0.2) +
  theme_bw(base_size = 15) +
  labs(
    x = "CYTOR expression (logCPM, normalized)",
    y = "MIR4435-2HG expression (logCPM, normalized)",
    title = "Correlation between CYTOR and MIR4435-2HG expression"
  ) +
  annotate(
    "text",
    x = min(df_corr$CYTOR, na.rm = TRUE),
    y = max(df_corr$MIR4435_2HG, na.rm = TRUE),
    hjust = 0, vjust = 1,
    label = paste0("Pearson's r = ", cor_val, "\n", p_lab),
    size = 5,
    color = "black"
  )

# 6. Save and show plot
ggsave("correlation_CYTOR_MIR44352HG.pdf", p, width = 6, height = 5)
print(p)


length(unique(ceRNA_triplets$mRNA))

# --- Correlation: CYTOR vs MIR4435-2HG in GBM Tumors ---

# Ensembl IDs
cytor_id <- "ENSG00000222041"
mir4435hg_id <- "ENSG00000172965"

# Check both genes exist
if (!(cytor_id %in% rownames(expr_mat)) | !(mir4435hg_id %in% rownames(expr_mat))) {
  stop("One or both lncRNAs not found in expression matrix.")
}

# Expression in tumor samples (already subset to tumor columns via expr_mat_tumor)
expr_cytor <- expr_mat[cytor_id, tumor_samples]
expr_mir4435hg <- expr_mat[mir4435hg_id, tumor_samples]

# Correlation test (Pearson)
cor_test <- cor.test(expr_cytor, expr_mir4435hg, method = "pearson")
cat("Pearson correlation between CYTOR and MIR4435-2HG:\n")
print(cor_test)

# Create plot
library(ggplot2)
library(tibble)

df_corr <- tibble(
  CYTOR = as.numeric(expr_cytor),
  MIR4435_2HG = as.numeric(expr_mir4435hg)
)

p_corr <- ggplot(df_corr, aes(x = CYTOR, y = MIR4435_2HG)) +
  geom_point(color = "#2c7fb8", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "#d95f0e", se = FALSE) +
  theme_bw(base_size = 14) +
  labs(
    title = "CYTOR vs MIR4435-2HG in GBM Tumor Samples",
    x = "CYTOR expression (logCPM)",
    y = "MIR4435-2HG expression (logCPM)"
  ) +
  annotate("text",
           x = min(df_corr$CYTOR, na.rm = TRUE),
           y = max(df_corr$MIR4435_2HG, na.rm = TRUE),
           label = paste0("r = ", round(cor_test$estimate, 3), "\n",
                          "p = ", signif(cor_test$p.value, 3)),
           hjust = 0, vjust = 1,
           size = 5)

ggsave("correlation_CYTOR_vs_MIR4435-2HG.pdf", plot = p_corr, width = 6, height = 5)


# --- Correlation: CYTOR vs MIR4435-2HG in GBM Tumors (with R²) ---

cytor_id <- "ENSG00000222041"
mir4435hg_id <- "ENSG00000172965"

if (!(cytor_id %in% rownames(expr_mat)) | !(mir4435hg_id %in% rownames(expr_mat))) {
  stop("One or both lncRNAs not found in expression matrix.")
}

expr_cytor <- expr_mat[cytor_id, tumor_samples]
expr_mir4435hg <- expr_mat[mir4435hg_id, tumor_samples]

cor_test <- cor.test(expr_cytor, expr_mir4435hg, method = "pearson")
r_val <- round(as.numeric(cor_test$estimate), 3)
r2_val <- round(r_val^2, 3)  # Calculate R²
p_lab <- ifelse(cor_test$p.value < 2e-16, "p < 2e-16", paste0("p = ", signif(cor_test$p.value, 3)))

df_corr <- tibble(
  CYTOR = as.numeric(expr_cytor),
  MIR4435_2HG = as.numeric(expr_mir4435hg)
)

# Create annotation label
corr_label <- paste0("r = ", r_val, 
                     "\nR² = ", r2_val, 
                     "\n", p_lab)

p_corr <- ggplot(df_corr, aes(x = CYTOR, y = MIR4435_2HG)) +
  geom_point(color = "#2c7fb8", alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", color = "#d95f0e", se = FALSE) +
  theme_bw(base_size = 14) +
  labs(
    title = "CYTOR vs MIR4435-2HG in GBM Tumor Samples",
    x = "CYTOR expression (logCPM)",
    y = "MIR4435-2HG expression (logCPM)"
  ) +
  annotate("text",
           x = min(df_corr$CYTOR, na.rm = TRUE),
           y = max(df_corr$MIR4435_2HG, na.rm = TRUE),
           label = corr_label,
           hjust = 0, vjust = 1,
           size = 5)

ggsave("correlation_CYTOR_vs_MIR4435-2HG.pdf", plot = p_corr, width = 6, height = 5)

