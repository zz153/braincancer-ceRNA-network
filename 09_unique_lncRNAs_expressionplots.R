## Title:plot_ceRNA_unique_lncRNA_expression.R
## Purpose: # Boxplots of ceRNA Hub lncRNAs in GBM and Normal Tissues

sample_annot <- data.frame(
  Sample = colnames(lnc_expr),
  Group = case_when(
    grepl("^TCGA.*-01", colnames(lnc_expr)) ~ "Tumor",
    grepl("^TCGA.*-11", colnames(lnc_expr)) ~ "Normal_TCGA",
    grepl("^GTEX", colnames(lnc_expr)) ~ "Normal_GTEx",
    TRUE ~ "Other"
  )
)

library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Annotate lncRNAs with gene symbols
lnc_annot <- gene_info %>%
  dplyr::filter(ensembl_gene_id %in% unique_lncRNAs) %>%
  dplyr::select(ensembl_gene_id, external_gene_name)


# 2. Subset and tidy
lnc_expr_long <- as.data.frame(lnc_expr) %>%
  mutate(ensembl_gene_id = rownames(.)) %>%
  pivot_longer(
    cols = -ensembl_gene_id,
    names_to = "Sample",
    values_to = "Expression"
  ) %>%
  left_join(sample_annot, by = "Sample") %>%
  left_join(lnc_annot, by = "ensembl_gene_id")

# Drop "Other" if you want only those three
lnc_expr_long <- lnc_expr_long %>% filter(Group %in% c("Tumor", "Normal_TCGA", "Normal_GTEx"))

p <- ggplot(lnc_expr_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ external_gene_name, scales = "free_y") +
  labs(
    title = "Tumor vs Normal (TCGA) vs Normal (GTEx): ceRNA lncRNAs",
    x = "",
    y = "Expression (logCPM)"
  ) +
  theme_bw(base_size = 12) +
  scale_fill_manual(values = c("Tumor" = "#e41a1c", "Normal_TCGA" = "#377eb8", "Normal_GTEx" = "#4daf4a")) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )
# PDF (vector, great for publication)
ggsave("lncRNA_ceRNA_boxplots.pdf", plot = p, width = 12, height = 8)

# PNG (raster, good for slides)
ggsave("lncRNA_ceRNA_boxplots.png", plot = p, width = 12, height = 8, dpi = 300)

# TIFF (high-res for journals)
ggsave("lncRNA_ceRNA_boxplots.tiff", plot = p, width = 12, height = 8, dpi = 600)

# ---- Boxplots of ceRNA lncRNAs with p-values ----

# 1. Required libraries
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
library(ggplot2)
library(ggpubr)
library(dplyr)

# 2. (Optional) Ensure Group is a factor with correct order
lnc_expr_long$Group <- factor(
  lnc_expr_long$Group,
  levels = c("Tumor", "Normal_TCGA", "Normal_GTEx")
)

# 3. Define pairwise comparisons
comparisons <- list(
  c("Tumor", "Normal_TCGA"),
  c("Tumor", "Normal_GTEx"),
  c("Normal_TCGA", "Normal_GTEx")
)

# 4. Generate the plot
p <- ggboxplot(
  lnc_expr_long,
  x = "Group", y = "Expression", fill = "Group",
  facet.by = "external_gene_name",
  palette = c("Tumor" = "#e41a1c", "Normal_TCGA" = "#377eb8", "Normal_GTEx" = "#4daf4a"),
  outlier.shape = NA
) +
  stat_compare_means(
    method = "kruskal.test", # global p-value per gene
    label = "p.format",
    label.y = 1.1 * max(lnc_expr_long$Expression, na.rm = TRUE)
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.signif",
    hide.ns = TRUE
  ) +
  labs(
    title = "Expression of ceRNA lncRNAs (boxplots + p-values)",
    x = "",
    y = "Expression (logCPM)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

# 5. Save the plot
ggsave("ceRNA_lncRNAs_boxplots_with_pvals.pdf", plot = p, width = 12, height = 8)


# ---- Save individual boxplots for each ceRNA lncRNA ----

library(ggplot2)
library(ggpubr)
library(dplyr)

# Ensure Group is a factor with the desired order
lnc_expr_long$Group <- factor(
  lnc_expr_long$Group,
  levels = c("Tumor", "Normal_TCGA", "Normal_GTEx")
)

# Get unique lncRNAs and gene symbols
genes <- unique(lnc_expr_long$external_gene_name)

comparisons <- list(
  c("Tumor", "Normal_TCGA"),
  c("Tumor", "Normal_GTEx"),
  c("Normal_TCGA", "Normal_GTEx")
)

for (g in genes) {
  # Subset data for the current gene
  gene_data <- lnc_expr_long %>% filter(external_gene_name == g)
  
  # Build plot
  p <- ggboxplot(
    gene_data,
    x = "Group", y = "Expression", fill = "Group",
    palette = c("Tumor" = "#e41a1c", "Normal_TCGA" = "#377eb8", "Normal_GTEx" = "#4daf4a"),
    outlier.shape = NA
  ) +
    stat_compare_means(method = "kruskal.test", label = "p.format") +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif", hide.ns = TRUE) +
    labs(
      title = paste0("Expression of ", g, " (ceRNA lncRNA)"),
      x = "",
      y = "Expression (logCPM)"
    ) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top")
  
  # Save as PDF (change to ".png" for PNGs)
  out_name <- paste0(gsub("[^A-Za-z0-9_]", "_", g), "_boxplot.pdf")
  ggsave(out_name, plot = p, width = 5, height = 6)
}

cat("Saved individual boxplots for each lncRNA!\n")
