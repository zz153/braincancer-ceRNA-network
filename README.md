# 🧬 braincancer-ceRNA-network

A reproducible pipeline for constructing and analyzing competing endogenous RNA (ceRNA) networks in glioblastoma (GBM) using RNA-seq data from TCGA and GTEx.

> 📄 **Published in Bioinformatics and Biology Insights**
> Rana Z & Grans J. *Comprehensive ceRNA Profiling Uncovers Clinically Relevant Hub lncRNAs in Glioblastoma.*
> https://doi.org/10.1177/11779322251411169

---

## 🔍 Overview

Competing endogenous RNAs (ceRNAs) regulate gene expression by sequestering shared microRNAs. This pipeline identifies clinically relevant hub lncRNAs in GBM through multi-cohort ceRNA network construction and survival analysis.

---

## ⚙️ Pipeline Structure

```
01_data_download.R                          # TCGA and GTEx data acquisition
02_datapreprocessing.R                      # Normalisation and harmonization
03_diff_expression.R                        # Differential expression analysis
04_download_interaction_datasets.sh         # miRNA-target interaction data download
05_extract_gene_info.R                      # Gene annotation and filtering
06_process_mirtarbase_ceRNA.R               # miRTarBase interaction processing
07_process_encori_ceRNA.R                   # ENCORI CLIP-seq interaction processing
08_ceRNA_triplet_construction.R             # ceRNA triplet (lncRNA-miRNA-mRNA) construction
09_unique_lncRNAs_expressionplots.R         # lncRNA expression visualisation
10_unique_lncRNAs_survivalplots.R           # Survival analysis for hub lncRNAs
11_prognostic_lncRNA_analysis_correlation_biopathways.R  # Pathway enrichment and correlation
12_survivalplots_ceRNA_mRNAs.R              # Survival plots for ceRNA-associated mRNAs
13_harmonizationsummary_QC.R               # Harmonization QC summary
```

---

## 🧠 Methods Summary

| Step | Description |
|---|---|
| Data | TCGA-GBM tumor + GTEx normal brain RNA-seq |
| Differential Expression | DESeq2-based tumor vs normal comparison |
| Interaction Data | miRTarBase + ENCORI CLIP-seq |
| ceRNA Network | lncRNA-miRNA-mRNA triplet construction |
| Survival Analysis | Cox regression and Kaplan–Meier stratification |

---

## 🚀 Getting Started

```r
# Step 1: Download data
source("01_data_download.R")

# Step 2: Preprocess
source("02_datapreprocessing.R")

# Continue sequentially through scripts 03-13
```

---

## 📬 Contact

**Zohaib Rana**
Postdoctoral Fellow, Department of Biochemistry, University of Otago
[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--5171--501X-A6CE39?style=flat-square&logo=orcid&logoColor=white)](https://orcid.org/0000-0002-5171-501X)
📧 zohaib.rana@otago.ac.nz

---

## 📜 Citation

If you use this pipeline, please cite:

> Rana Z & Grans J. Comprehensive ceRNA Profiling Uncovers Clinically Relevant Hub lncRNAs in Glioblastoma. *Bioinformatics and Biology Insights.* 2025. https://doi.org/10.1177/11779322251411169

---

*University of Otago · Department of Biochemistry · Dunedin, NZ*
