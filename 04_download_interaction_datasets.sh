#!/bin/bash

# ---- miRTarBase 2025 ----
echo "Downloading miRTarBase 2025..."
curl -L 'https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/data/miRTarBase_MTI.xlsx' -o miRTarBase_MTI_2025.xlsx

# ---- ENCORI lncRNA–miRNA CLIP-supported ----
echo "Downloading ENCORI lncRNA–miRNA CLIP-supported pairs..."
curl 'https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=lncRNA&miRNA=all&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=0&program=None&target=all&cellType=all' -o ENCORI_lncRNA_miRNA_CLIP.txt

echo "All files downloaded!"
