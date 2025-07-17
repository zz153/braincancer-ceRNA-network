# scripts/02_preprocessing.R
# -----------------------------------------------------------
# Preprocessing: Extract, clean, harmonize TCGA/GTEx expression
# Save ready-to-analyze count matrices and annotation.
# -----------------------------------------------------------

# Load libraries
libs <- c("SummarizedExperiment", "dplyr", "rtracklayer", "tidyverse")
invisible(lapply(libs, require, character.only = TRUE))

# --- Load raw SummarizedExperiment objects ---
tcga_data   <- readRDS("TCGA_GBM_data.rds")
mir_data    <- readRDS("TCGA_GBM_miRNA_data.rds")
rse_brain   <- readRDS("GTEx_Brain_rse.rds")

# --- Clean TCGA gene counts (remove version numbers) ---
tcga_counts <- assay(tcga_data)
rownames(tcga_counts) <- gsub("\\..*", "", rownames(tcga_counts))

# --- Extract TCGA sample metadata ---
tcga_meta <- colData(tcga_data)
tcga_tumor_samples  <- rownames(tcga_meta[tcga_meta$sample_type == "Primary Tumor", ])
tcga_normal_samples <- rownames(tcga_meta[tcga_meta$sample_type == "Solid Tissue Normal", ])

# --- Prepare GTEx gene expression (e.g., scale counts to CPM) ---
scale_counts <- function(rse) {
  round(t(t(assay(rse)) / colSums(assay(rse))) * 1e6)
}
gtex_expr <- scale_counts(rse_brain)
rownames(gtex_expr) <- gsub("\\..*", "", rownames(gtex_expr))

# --- Harmonize gene IDs across datasets ---
common_genes <- Reduce(intersect, list(
  rownames(tcga_counts), rownames(gtex_expr)
))
tcga_counts_h <- tcga_counts[common_genes, ]
gtex_expr_h   <- gtex_expr[common_genes, ]

# --- Extract GENCODE annotation, save gene info as .rds for easy loading ---
gtf_file <- "gencode.v38.annotation.gtf"
gtf <- rtracklayer::import(gtf_file)
genes_gtf <- gtf[gtf$type == "gene"]
gene_info <- data.frame(
  ensembl_gene_id    = gsub("\\..*", "", genes_gtf$gene_id),
  gene_biotype       = genes_gtf$gene_type,
  external_gene_name = genes_gtf$gene_name
)
saveRDS(gene_info, "gene_info.rds")

# --- Save harmonized, clean matrices for downstream DE ---
saveRDS(tcga_counts_h, "TCGA_GBM_counts_harmonized.rds")
saveRDS(gtex_expr_h,   "GTEx_Brain_counts_harmonized.rds")

# --- Optional: Save lncRNA/mRNA split for quick analysis ---
mRNA_ids  <- gene_info$ensembl_gene_id[gene_info$gene_biotype == "protein_coding"]
lncRNA_ids<- gene_info$ensembl_gene_id[gene_info$gene_biotype == "lncRNA"]

saveRDS(tcga_counts_h[mRNA_ids, , drop=FALSE],  "TCGA_GBM_mRNA_counts.rds")
saveRDS(tcga_counts_h[lncRNA_ids, , drop=FALSE],"TCGA_GBM_lncRNA_counts.rds")
saveRDS(gtex_expr_h[mRNA_ids, , drop=FALSE],    "GTEx_Brain_mRNA_counts.rds")
saveRDS(gtex_expr_h[lncRNA_ids, , drop=FALSE],  "GTEx_Brain_lncRNA_counts.rds")

# --- Process miRNA matrix (optional: adapt to your exact format) ---
# You might want to harmonize barcodes, remove lowly expressed miRNAs, etc.
# Save miRNA count matrix for DE
saveRDS(mir_data, "TCGA_GBM_miRNA_data_cleaned.rds")

cat("Preprocessing and harmonization complete!\n")
