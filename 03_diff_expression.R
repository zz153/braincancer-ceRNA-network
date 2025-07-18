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
tcga_data   <- readRDS("TCGA_GBM_data.rds")
rse_brain   <- readRDS("GTEx_Brain_rse.rds")
gene_info   <- readRDS("gene_info.rds")
mir_data    <- readRDS("TCGA_GBM_miRNA_data.rds")

# --- Make function for GTEx counts normalization ---
scale_counts <- function(rse) {
  round(t(t(assay(rse)) / colSums(assay(rse))) * 1e6)
}

# --- Build Expression Matrices ---
tcga_counts <- assay(tcga_data)
rownames(tcga_counts) <- gsub("\\..*", "", rownames(tcga_counts))
gtex_expr <- scale_counts(rse_brain)
rownames(gtex_expr) <- gsub("\\..*", "", rownames(gtex_expr))

# --- Harmonize genes across datasets ---
common_genes <- Reduce(intersect, list(
  rownames(tcga_counts), rownames(gtex_expr)
))
tcga_counts_h <- tcga_counts[common_genes, ]
gtex_expr_h   <- gtex_expr[common_genes, ]

# --- Identify Tumor, TCGA-Normal, and GTEx-Normal columns ---
tcga_meta <- colData(tcga_data)
tcga_tumor_samples  <- rownames(tcga_meta[tcga_meta$sample_type == "Primary Tumor", ])
tcga_normal_samples <- rownames(tcga_meta[tcga_meta$sample_type == "Solid Tissue Normal", ])
tcga_tumor_expr     <- tcga_counts_h[, tcga_tumor_samples, drop=FALSE]
tcga_normal_expr    <- tcga_counts_h[, tcga_normal_samples, drop=FALSE]
gtex_normal_expr    <- gtex_expr_h

# --- Combine all columns for DE analysis ---
counts_combined <- cbind(tcga_tumor_expr, tcga_normal_expr, gtex_normal_expr)
sample_names    <- colnames(counts_combined)
group_labels <- c(
  rep("Tumor",       ncol(tcga_tumor_expr)),
  rep("TCGA_Normal", ncol(tcga_normal_expr)),
  rep("GTEx_Normal", ncol(gtex_normal_expr))
)
group <- factor(group_labels, levels = c("TCGA_Normal", "GTEx_Normal", "Tumor"))
names(group) <- sample_names
batch <- factor(c(
  rep("TCGA", ncol(tcga_tumor_expr) + ncol(tcga_normal_expr)),
  rep("GTEx", ncol(gtex_normal_expr))
), levels = c("TCGA", "GTEx"))
names(batch) <- sample_names

# --- edgeR DGEList, filtering, normalization ---
dge <- DGEList(counts = counts_combined, group = group)
keep <- filterByExpr(dge)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge)

# --- Design matrix & voom ---
design <- model.matrix(~ group + batch)
v      <- voom(dge, design, plot = TRUE)

# --- PCA (Quality Control) ---
cols <- c("TCGA_Normal" = "#0072B2", "GTEx_Normal" = "#009E73", "Tumor" = "#D55E00")
col_vector <- cols[as.character(group[colnames(v$E)])]
if (any(is.na(col_vector))) {
  warning("Some samples missing group/color! Using grey for those.")
  col_vector[is.na(col_vector)] <- "grey"
}
pca <- prcomp(t(v$E))
plot(
  pca$x[,1], pca$x[,2],
  col = col_vector,
  pch = 19,
  xlab = "PC1", ylab = "PC2",
  main = "PCA: All Samples (after normalization)"
)
legend("bottomright", legend = names(cols), fill = cols)

# --- Boxplots by sample ---
box_cols <- col_vector
if (length(box_cols) != ncol(v$E)) {
  box_cols <- rep("grey", ncol(v$E))
}
boxplot(v$E, col = box_cols, las = 2, main = "logCPM by Sample", outline = FALSE)

# --- Differential Expression (limma) ---
fit    <- lmFit(v, design)
fit    <- eBayes(fit)
res_all <- topTable(fit, coef = "groupTumor", number = Inf, adjust.method = "BH")
res_sig <- res_all %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)

write.csv(res_all, "DEG_all_Tumor_vs_AllNormals_batchCorrected.csv")
write.csv(res_sig, "DEG_sig_Tumor_vs_AllNormals_batchCorrected.csv")

# --- Volcano Plot for mRNA/lncRNA ---
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
ann_colors <- list(Group = cols)
if (nrow(annotation_col) == 0 | ncol(expr_mat_z) == 0 | nrow(expr_mat_z) == 0) {
  cat("Skipping heatmap: no significant genes or groups\n")
} else {
  pheatmap(
    expr_mat_z,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("blue", "black", "yellow"))(100),
    show_colnames = FALSE, show_rownames = FALSE,
    cluster_cols = TRUE, cluster_rows = TRUE,
    fontsize = 6,
    main = "DE genes (z-scored)",
    breaks = seq(-2, 2, length.out = 101)
  )
}

# ==================================================
# ========== Repeat for miRNA ======================
# ==================================================

# --- Extract and prepare miRNA matrix ---
# Handles SummarizedExperiment or data.frame!
if ("SummarizedExperiment" %in% class(mir_data)) {
  read_count_cols <- grep("^read_count_", colnames(assay(mir_data)), value = TRUE)
  mir_expr <- as.matrix(assay(mir_data)[read_count_cols, ])
  if (nrow(mir_expr) == 0) mir_expr <- t(assay(mir_data)[read_count_cols, ])
  colnames(mir_expr) <- gsub("^read_count_", "", colnames(mir_expr))
  rownames(mir_expr) <- rowData(mir_data)$miRNA_ID
  # Try to extract groups from colData if present
  mir_coldata <- as.data.frame(colData(mir_data))
  if ("sample_type" %in% colnames(mir_coldata)) {
    mir_sample_types <- mir_coldata$sample_type
    names(mir_sample_types) <- rownames(mir_coldata)
    mir_group <- mir_sample_types[colnames(mir_expr)]
    mir_group <- factor(ifelse(mir_group == "Primary Tumor", "Tumor", "Normal"))
  } else {
    mir_group <- ifelse(grepl("01A", colnames(mir_expr)), "Tumor", "Normal")
    mir_group <- factor(mir_group)
  }
} else {
  # Assume mir_data is a data.frame with read_count columns
  read_count_cols <- grep("^read_count_", colnames(mir_data), value = TRUE)
  mir_expr <- as.matrix(mir_data[, read_count_cols])
  colnames(mir_expr) <- gsub("^read_count_", "", colnames(mir_expr))
  rownames(mir_expr) <- mir_data$miRNA_ID
  # Guess group from column names
  mir_group <- ifelse(grepl("01A", colnames(mir_expr)), "Tumor", "Normal")
  mir_group <- factor(mir_group)
}

# --- DGE for miRNA ---
dge_mir <- DGEList(counts = mir_expr)
keep_mir <- filterByExpr(dge_mir)
dge_mir <- dge_mir[keep_mir, , keep.lib.sizes = FALSE]
dge_mir <- calcNormFactors(dge_mir)
design_mir <- model.matrix(~ mir_group)
v_mir <- voom(dge_mir, design_mir, plot = TRUE)
fit_mir <- lmFit(v_mir, design_mir)
fit_mir <- eBayes(fit_mir)
res_mir <- topTable(fit_mir, coef = "mir_groupTumor", number = Inf, adjust.method = "BH")
write.csv(res_mir, "DE_miRNA_Tumor_vs_Normal.csv")

# --- Volcano plot for miRNA ---
df_mir <- res_mir %>%
  rownames_to_column("miRNA") %>%
  mutate(
    category = case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Up",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    label = miRNA
  )
pal_mir <- c(Up = "#D55E00", Down = "#0072B2", NS = "grey70")
top_hits_mir <- df_mir %>% filter(category != "NS") %>% arrange(adj.P.Val) %>% slice_head(n = 10)
p_mir <- ggplot(df_mir, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = category), alpha = 0.7, size = 1.6) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text_repel(data = top_hits_mir, aes(label = label), size = 3, max.overlaps = 10) +
  scale_color_manual(values = pal_mir, name = NULL) +
  labs(
    title = "Volcano: miRNAs (Tumor vs Normal)",
    x = expression(log[2]~fold-change),
    y = expression(-log[10]~"(adj. P-value)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
ggsave("volcano_miRNA_labeled.pdf", p_mir, width = 6, height = 5)

# --- Save workspace and session info ---
save.image("03_diff_expression_workspace.RData")
writeLines(capture.output(sessionInfo()), "sessionInfo_03_diff_expression.txt")

cat("Differential expression analysis complete!\n")
