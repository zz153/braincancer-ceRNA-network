# scripts/02_data_preprocess.R
# ---------------------------------------------------
# Harmonize and extract mRNA/lncRNA matrices for TCGA-GBM & GTEx
# ---------------------------------------------------

# --- Load libraries ---
libs <- c("SummarizedExperiment", "rtracklayer", "dplyr")
invisible(lapply(libs, require, character.only = TRUE))

# --- Load Data ---
tcga_data <- readRDS("TCGA_GBM_data.rds")
rse_brain <- readRDS("GTEx_Brain_rse.rds")

# --- Load & Process Annotation ---
gtf <- rtracklayer::import("gencode.v38.annotation.gtf")
genes_gtf <- gtf[gtf$type == "gene"]
gene_info <- data.frame(
  ensembl_gene_id   = gsub("\\..*", "", genes_gtf$gene_id),
  gene_biotype      = genes_gtf$gene_type,
  external_gene_name = genes_gtf$gene_name,
  stringsAsFactors = FALSE
)
saveRDS(gene_info, file = "gene_info.rds")

# --- Get counts ---
tcga_counts <- assay(tcga_data)
rownames(tcga_counts) <- gsub("\\..*", "", rownames(tcga_counts))

# Simple function for GTEx normalization to CPM
scale_counts <- function(rse) {
  round(t(t(assay(rse)) / colSums(assay(rse))) * 1e6)
}
gtex_expr <- scale_counts(rse_brain)
rownames(gtex_expr) <- gsub("\\..*", "", rownames(gtex_expr))

# --- Harmonize gene sets ---
common_genes <- intersect(rownames(tcga_counts), rownames(gtex_expr))

tcga_counts_h <- tcga_counts[common_genes, ]
gtex_expr_h   <- gtex_expr[common_genes, ]

# --- Identify mRNA and lncRNA gene IDs using GENCODE annotation ---
mRNA_ids   <- gene_info$ensembl_gene_id[gene_info$gene_biotype == "protein_coding"]
lncRNA_ids <- gene_info$ensembl_gene_id[gene_info$gene_biotype == "lncRNA"]

# Only keep IDs present in the matrices
mRNA_ids_in_counts   <- intersect(mRNA_ids, rownames(tcga_counts_h))
lncRNA_ids_in_counts <- intersect(lncRNA_ids, rownames(tcga_counts_h))

# --- Extract & save the matrices (TCGA) ---
tcga_mrna   <- tcga_counts_h[mRNA_ids_in_counts, , drop = FALSE]
tcga_lncrna <- tcga_counts_h[lncRNA_ids_in_counts, , drop = FALSE]

saveRDS(tcga_mrna,   "TCGA_GBM_mRNA_counts.rds")
saveRDS(tcga_lncrna, "TCGA_GBM_lncRNA_counts.rds")

# --- Extract & save the matrices (GTEx) ---
gtex_mrna   <- gtex_expr_h[mRNA_ids_in_counts, , drop = FALSE]
gtex_lncrna <- gtex_expr_h[lncRNA_ids_in_counts, , drop = FALSE]

saveRDS(gtex_mrna,   "GTEx_Brain_mRNA_counts.rds")
saveRDS(gtex_lncrna, "GTEx_Brain_lncRNA_counts.rds")

# --- Quick sanity summary ---
cat("Genes in TCGA-GBM:", nrow(tcga_counts), "\n")
cat("Genes in GTEx:", nrow(gtex_expr), "\n")
cat("Common genes:", length(common_genes), "\n")
cat("Protein-coding genes extracted:", length(mRNA_ids_in_counts), "\n")
cat("lncRNAs extracted:", length(lncRNA_ids_in_counts), "\n")
cat("All harmonized matrices saved!\n")

