# --- Load harmonized expression matrices ---
tcga_data   <- readRDS("TCGA_GBM_data.rds")
rse_brain   <- readRDS("GTEx_Brain_rse.rds")

# --- Get expression matrices and harmonize ---
tcga_counts <- assay(tcga_data)
rownames(tcga_counts) <- gsub("\\..*", "", rownames(tcga_counts))
scale_counts <- function(rse) round(t(t(assay(rse)) / colSums(assay(rse))) * 1e6)
gtex_expr <- scale_counts(rse_brain)
rownames(gtex_expr) <- gsub("\\..*", "", rownames(gtex_expr))

# --- Harmonize genes across datasets ---
common_genes <- intersect(rownames(tcga_counts), rownames(gtex_expr))
tcga_counts_h <- tcga_counts[common_genes, ]
gtex_expr_h   <- gtex_expr[common_genes, ]

# --- Identify sample columns ---
tcga_meta <- colData(tcga_data)
tcga_tumor_samples  <- rownames(tcga_meta[tcga_meta$sample_type == "Primary Tumor", ])
tcga_normal_samples <- rownames(tcga_meta[tcga_meta$sample_type == "Solid Tissue Normal", ])
tcga_tumor_expr     <- tcga_counts_h[, tcga_tumor_samples, drop=FALSE]
tcga_normal_expr    <- tcga_counts_h[, tcga_normal_samples, drop=FALSE]
gtex_normal_expr    <- gtex_expr_h

# --- Combine all for DE analysis ---
counts_combined <- cbind(tcga_tumor_expr, tcga_normal_expr, gtex_normal_expr)

# --- Print summary ---
cat("=== Harmonized/Combined matrix used for mRNA/lncRNA DE analysis ===\n")
cat("Genes (rows):    ", nrow(counts_combined), "\n")
cat("Samples (cols):  ", ncol(counts_combined), "\n")
cat(" - TCGA Tumor:   ", ncol(tcga_tumor_expr), "\n")
cat(" - TCGA Normal:  ", ncol(tcga_normal_expr), "\n")
cat(" - GTEx Normal:  ", ncol(gtex_normal_expr), "\n")

# --- Load TCGA miRNA data ---
mir_data <- readRDS("TCGA_GBM_miRNA_data.rds")

# --- Extract miRNA count matrix and group info (handles SummarizedExperiment and data.frame) ---
if ("SummarizedExperiment" %in% class(mir_data)) {
  read_count_cols <- grep("^read_count_", colnames(assay(mir_data)), value = TRUE)
  mir_expr <- as.matrix(assay(mir_data)[read_count_cols, ])
  if (nrow(mir_expr) == 0) mir_expr <- t(assay(mir_data)[read_count_cols, ])
  colnames(mir_expr) <- gsub("^read_count_", "", colnames(mir_expr))
  rownames(mir_expr) <- rowData(mir_data)$miRNA_ID
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
  read_count_cols <- grep("^read_count_", colnames(mir_data), value = TRUE)
  mir_expr <- as.matrix(mir_data[, read_count_cols])
  colnames(mir_expr) <- gsub("^read_count_", "", colnames(mir_expr))
  rownames(mir_expr) <- mir_data$miRNA_ID
  mir_group <- ifelse(grepl("01A", colnames(mir_expr)), "Tumor", "Normal")
  mir_group <- factor(mir_group)
}

# --- Print summary ---
cat("=== TCGA miRNA expression matrix used for DE analysis ===\n")
cat("miRNAs (rows):     ", nrow(mir_expr), "\n")
cat("Samples (columns): ", ncol(mir_expr), "\n")
cat(" - Tumor samples:  ", sum(mir_group == "Tumor"), "\n")
cat(" - Normal samples: ", sum(mir_group == "Normal"), "\n")

# Make a data.frame summarizing your sample/gene/miRNA counts
summary_table <- data.frame(
  DataType     = c("mRNA/lncRNA", "mRNA/lncRNA", "mRNA/lncRNA", "mRNA/lncRNA", "miRNA", "miRNA", "miRNA"),
  Subset       = c("All Samples", "TCGA Tumor", "TCGA Normal", "GTEx Normal", "All Samples", "TCGA Tumor", "TCGA Normal"),
  Features     = c(57562, 57562, 57562, 57562, 1881, 1881, 1881),
  Samples      = c(3308, 372, 5, 2931, 277, 240, 37)
)
install.packages("kableExtra")
# Clean up for a nicer table display
library(knitr)
library(kableExtra)

kable(summary_table, align = 'c', caption = "Summary of harmonized and combined datasets for differential expression analysis") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

# 1. Make per-sample summary for your expression matrix (genes x samples)
summary_df <- data.frame(
  sample = colnames(v$E),  # or colnames(log2_cpm_mir) for miRNA
  median = apply(v$E, 2, median, na.rm = TRUE),
  Q1     = apply(v$E, 2, function(x) quantile(x, 0.25, na.rm = TRUE)),
  Q3     = apply(v$E, 2, function(x) quantile(x, 0.75, na.rm = TRUE))
)
summary_df$group <- group[summary_df$sample]

# 2. Keep only TCGA GBM & TCGA Normal
summary_df_tcga <- summary_df[summary_df$group %in% c("TCGA_Normal", "Tumor"), ]
summary_df_tcga$group <- factor(summary_df_tcga$group, levels = c("TCGA_Normal", "Tumor"))
summary_df_tcga <- summary_df_tcga[order(summary_df_tcga$group), ]

# 3. Plot
library(ggplot2)
ggplot(summary_df_tcga, aes(x = factor(sample, levels = sample), y = median, color = group)) +
  geom_errorbar(aes(ymin = Q1, ymax = Q3), width = 0.2) +
  geom_point(size = 0.7) +
  scale_color_manual(values = c("TCGA_Normal" = "#0072B2", "Tumor" = "#D55E00")) +
  labs(
    title = expression("Normalized log"[2]*"-CPM: TCGA GBM vs Normal"),
    x     = NULL,
    y     = expression(log[2]*"-CPM")
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.grid.major = element_line(linewidth = 0.1),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5)
  )

# Your v$E: genes x samples
summary_df <- data.frame(
  sample = colnames(v$E),
  median = apply(v$E, 2, median, na.rm=TRUE),
  group  = group[colnames(v$E)]
)

library(ggplot2)
ggplot(summary_df, aes(x = group, y = median, fill = group)) +
  geom_boxplot(outlier.size=0.5, lwd=0.1, alpha=0.7) +
  scale_fill_manual(values = c("TCGA_Normal"="#0072B2", "GTEx_Normal"="#009E73", "Tumor"="#D55E00")) +
  labs(title="Sample medians by group", y="Sample median log2CPM", x=NULL) +
  theme_bw(base_size=12)

library(ggplot2)
ggplot(summary_df, aes(x = group, y = median, fill = group)) +
  geom_boxplot(outlier.size = 0.8, lwd = 0.2, alpha = 0.75, width = 0.65) +
  scale_fill_manual(values = c("TCGA_Normal"="#0072B2", "GTEx_Normal"="#009E73", "Tumor"="#D55E00")) +
  labs(
    title = "Sample medians by group",
    y = "Sample median log2CPM",
    x = NULL
  ) +
  theme_bw(base_size = 18) +  # This increases base font size
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 22),
    axis.title.y    = element_text(size = 20, face = "bold"),
    axis.text.x     = element_text(size = 18, face = "bold"),
    axis.text.y     = element_text(size = 18),
    legend.title    = element_text(size = 16, face = "bold"),
    legend.text     = element_text(size = 15),
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"),
    panel.grid.minor = element_blank(),
    panel.border    = element_rect(size = 1),
    legend.position = "none"
  )

summary_df <- data.frame(
  sample = colnames(v$E),
  median = apply(v$E, 2, median, na.rm = TRUE),
  Q1     = apply(v$E, 2, function(x) quantile(x, 0.25, na.rm = TRUE)),
  Q3     = apply(v$E, 2, function(x) quantile(x, 0.75, na.rm = TRUE))
)
summary_df$group <- group[summary_df$sample]

str(summary_df)

library(dplyr)
iqr_table <- summary_df %>%
  group_by(group) %>%
  summarise(
    median_Q1 = median(Q1, na.rm=TRUE),
    median_Q3 = median(Q3, na.rm=TRUE),
    IQR = median_Q3 - median_Q1
  )

# For annotation
iqr_table$label <- paste0("IQR = ", round(iqr_table$IQR, 2))

library(ggplot2)
ggplot(summary_df, aes(x = group, y = median, fill = group)) +
  geom_boxplot(outlier.size = 0.8, lwd = 0.2, alpha = 0.75, width = 0.65) +
  scale_fill_manual(values = c("TCGA_Normal"="#0072B2", "GTEx_Normal"="#009E73", "Tumor"="#D55E00")) +
  labs(
    title = "Sample medians by group",
    y = "Sample median log2CPM",
    x = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 22),
    axis.title.y    = element_text(size = 20, face = "bold"),
    axis.text.x     = element_text(size = 18, face = "bold"),
    axis.text.y     = element_text(size = 18),
    legend.position = "none"
  ) +
  geom_text(
    data = iqr_table,
    aes(x = group, y = median_Q1 - 0.2, label = label),
    inherit.aes = FALSE,
    color = "black",
    size = 5,
    fontface = "italic"
  )

library(ggplot2)

ggplot(summary_df, aes(x = group, y = median, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.6, width = 0.9, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.8, outlier.alpha = 0.4, 
               fill = "white", color = "black", lwd = 0.4, alpha = 0.8) +
  scale_fill_manual(values = c("TCGA_Normal" = "#0072B2",
                               "GTEx_Normal" = "#009E73",
                               "Tumor" = "#D55E00")) +
  labs(
    title = "Per-sample expression median distributions",
    subtitle = "Sample median log2CPM by group (with violin + IQR boxplot)",
    x = NULL,
    y = "Sample median log2CPM"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 19),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "gray30"),
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 15),
    legend.position = "none"
  )

# 1. Calculate per-sample IQR
sample_iqr <- apply(v$E, 2, function(x) quantile(x, 0.75) - quantile(x, 0.25))

# 2. Create a dataframe with group info
sample_iqr_df <- data.frame(
  sample = names(sample_iqr),
  IQR = sample_iqr,
  group = group[names(sample_iqr)]
)

library(ggplot2)

ggplot(sample_iqr_df, aes(x = group, y = IQR, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.9, color = NA) +
  geom_boxplot(width = 0.15, outlier.size = 0.8, fill = "white", color = "black", lwd = 0.4, alpha = 0.8) +
  scale_fill_manual(values = c("TCGA_Normal"="#0072B2", "GTEx_Normal"="#009E73", "Tumor"="#D55E00")) +
  labs(
    title = "Per-sample IQRs by group",
    y = "Sample IQR (log2-CPM)",
    x = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 15),
    legend.position = "none"
  )

# Calculate per-sample IQR
sample_iqr <- apply(v$E, 2, function(x) quantile(x, 0.75) - quantile(x, 0.25))
sample_iqr_df <- data.frame(
  sample = names(sample_iqr),
  IQR = sample_iqr,
  group = group[names(sample_iqr)]
)

# Plot as boxplot
library(ggplot2)
ggplot(sample_iqr_df, aes(x = group, y = IQR, fill = group)) +
  geom_boxplot(outlier.size = 0.8, lwd = 0.2, alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c("TCGA_Normal"="#0072B2", "GTEx_Normal"="#009E73", "Tumor"="#D55E00")) +
  scale_y_continuous(limits = c(0, 10)) +  # y-axis from 0 to 10
  labs(
    title = "Sample-wise IQR (log2-CPM) distributions by group",
    y = "Sample IQR (log2-CPM)",
    x = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
    axis.title.y = element_text(size = 17),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 15),
    legend.position = "none"
  )
