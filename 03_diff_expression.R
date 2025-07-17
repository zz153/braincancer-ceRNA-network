# scripts/03_diff_expression.R
# -----------------------------------------------------------
# Differential Expression Analysis for mRNA, lncRNA, miRNA
# in TCGA-GBM vs. Normal (TCGA + GTEx)
# -----------------------------------------------------------

# --- Load Required Libraries ---
libs <- c(
  "TCGAbiolinks", "SummarizedExperiment", "edgeR", "limma", "biomaRt",
  "clusterProfiler", "org.Hs.eg.db", "tidyverse", "rtracklayer",
  "R.utils", "miRBaseConverter", "recount3", "ggpubr", "pheatmap",
  "survival", "survminer", "Hmisc", "STRINGdb", "tibble", "ggrepel"
)
invisible(lapply(libs, require, character.only = TRUE))

# --- Load Preprocessed Data ---
tcga_data <- readRDS("TCGA_GBM_data.rds")
rse_brain <- readRDS("GTEx_Brain_rse.rds")
gene_info <- readRDS("gene_info.rds") # Save this in preprocessing

# --- Expression Matrices ---
# Prepare gene expression matrices (simplified, you may load these from RDS if already preprocessed)
tcga_counts <- assay(tcga_data)
rownames(tcga_counts) <- gsub("\\..*", "", rownames(tcga_counts))
gtex_expr <- scale_counts(rse_brain) # Assuming you have a function scale_counts
rownames(gtex_expr) <- gsub("\\..*", "", rownames(gtex_expr))

# --- Harmonize genes across datasets ---
common_genes <- Reduce(intersect, list(
  rownames(tcga_counts), rownames(gtex_expr)
))
tcga_counts_h <- tcga_counts[common_genes, ]
gtex_expr_h <- gtex_expr[common_genes, ]

# --- Build group and batch factors ---
tcga_meta <- colData(tcga_data)
tcga_tumor_samples <- rownames(tcga_meta[tcga_meta$sample_type == "Primary Tumor", ])
tcga_normal_samples <- rownames(tcga_meta[tcga_meta$sample_type == "Solid Tissue Normal", ])
tcga_tumor_expr <- tcga_counts_h[, tcga_tumor_samples]
tcga_normal_expr <- tcga_counts_h[, tcga_normal_samples]
gtex_normal_expr <- gtex_expr_h

# --- Combine data for DE analysis ---
counts_combined <- cbind(tcga_tumor_expr, tcga_normal_expr, gtex_normal_expr)
group <- factor(c(
  rep("Tumor",       ncol(tcga_tumor_expr)),
  rep("TCGA_Normal", ncol(tcga_normal_expr)),
  rep("GTEx_Normal", ncol(gtex_normal_expr))
), levels = c("TCGA_Normal", "GTEx_Normal", "Tumor"))
batch <- factor(c(
  rep("TCGA", ncol(tcga_tumor_expr) + ncol(tcga_normal_expr)),
  rep("GTEx", ncol(gtex_normal_expr))
), levels = c("TCGA", "GTEx"))

# --- edgeR DGEList, filtering, normalization ---
dge <- DGEList(counts = counts_combined, group = group)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# --- Design matrix & voom ---
design <- model.matrix(~ group + batch)
v <- voom(dge, design, plot = TRUE) # mean-variance plot

# --- PCA (Quality Control) ---
pca <- prcomp(t(v$E))
cols <- c("TCGA_Normal" = "#0072B2", "GTEx_Normal" = "#009E73", "Tumor" = "#D55E00")
plot(
  pca$x[,1], pca$x[,2],
  col = cols[as.character(group)],
  pch = 19,
  xlab = "PC1", ylab = "PC2",
  main = "PCA: All Samples (after normalization)"
)
legend("bottomright", legend = names(cols), fill = cols)

# --- Boxplots by sample ---
boxplot(v$E, col = cols[as.character(group)], las = 2, main = "logCPM by Sample", outline = FALSE)

# --- Differential Expression (limma) ---
fit <- lmFit(v, design)
fit <- eBayes(fit)
res_all <- topTable(fit, coef = "groupTumor", number = Inf, adjust.method = "BH")
res_sig <- res_all %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)

write.csv(res_all, "DEG_all_Tumor_vs_AllNormals_batchCorrected.csv")
write.csv(res_sig, "DEG_sig_Tumor_vs_AllNormals_batchCorrected.csv")

# --- Volcano Plot ---
df_mrna <- res_all %>%
  rownames_to_column("ensembl_gene_id") %>%
  left_join(gene_info[, c("ensembl_gene_id", "external_gene_name")], by = "ensembl_gene_id") %>%
  mutate(
    category = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    label = external_gene_name
  )
pal <- c(Up = "#D55E00", Down = "#0072B2", NS = "grey70")
top_hits <- df_mrna %>% filter(category != "NS") %>% arrange(adj.P.Val) %>% slice_head(n = 10)
p_mrna <- ggplot(df_mrna, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = category), alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(data = top_hits, aes(label = label), size = 3, max.overlaps = 10) +
  scale_color_manual(values = pal, name = NULL) +
  labs(
    title = "Volcano: mRNAs (Tumor vs All Normals)",
    x = expression(log[2]~fold-change),
    y = expression(-log[10]~"(adj. P-value)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
ggsave("volcano_mRNA_labeled.pdf", p_mrna, width = 6, height = 5)

# --- Heatmap of significant DEGs ---
sig_genes <- rownames(res_sig)
expr_mat <- v$E[sig_genes, , drop = FALSE]
expr_mat_z <- t(scale(t(expr_mat)))
annotation_col <- data.frame(Group = group[colnames(expr_mat_z)])
rownames(annotation_col) <- colnames(expr_mat_z)
pheatmap(
  expr_mat_z,
  annotation_col = annotation_col,
  color = colorRampPalette(c("blue", "black", "yellow"))(100),
  show_colnames = FALSE, show_rownames = FALSE,
  cluster_cols = TRUE, cluster_rows = TRUE,
  fontsize = 6,
  main = "DE genes (z-scored)",
  breaks = seq(-2, 2, length.out = 101)
)

# --- Repeat for miRNA if available ---
# (Load/process miRNA as you did for mRNA/lncRNA: normalization, limma, volcano, etc.)

# --- Save workspace and session info ---
save.image("03_diff_expression_workspace.RData")
writeLines(capture.output(sessionInfo()), "sessionInfo_03_diff_expression.txt")

cat("Differential expression analysis complete!\n")
