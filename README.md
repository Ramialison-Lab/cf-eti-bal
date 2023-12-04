# CF-ETI-scRNA

## Description
Computational analyses of single-cell RNA-seq data from bronchoalveolar lavage (BAL) of preschool cystic fibrosis receiving Elexacaftor-Texacaftor-Ivacaftor (ETI) therapy.
Single-cell RNA-seq libraries were generated using the Cell Hashing with 10x Single Cell 3' Reagent Kit v3.1 (Dual Index) Protocol (Stoeckius et al. 2018).  
  
There are two batches in this dataset.  
For batch 1, ~110,000 cells from 10 hashed samples (each with a distinct Hashtag Antibody, or HTO) were pooled into one tube to generate four captures.  
For batch 2, ~60,000 cells from 6 hashed samples were pooled to generate two captures.  
QC and demultiplexing were carried out separately for each capture. Two batches were integrated into one dataset for clustering and downstream analysis.
  
Code included in the "src" directory:  
1. Shell script for cellranger multi  
2. R script for QC, HTO demultiplexing and preprocessing with reference to Maksimovic et al. (2022)  
3. R script for reference-based Human Lung Cell Atlas v2 annotation using Azimuth
4. R script for separate sub-clustering of Monocytes and Macrophages, T and NK cells, and other cells.
5. R script for differential expression and pathway enrichment of recruited lung monocytes & macrophages
6. R script for pseudotime analysis of recruited lung monocytes & macrophages  

Developers: Anson Wong (main), Hieu T Nim (contributor), Ramialison Laboratory, Australian Regenerative Medicine Institute & Murdoch Children's Research Institute, Australia  

## References
Stoeckius, M., Zheng, S., Houck-Loomis, B. et al. Cell Hashing with barcoded antibodies enables multiplexing and doublet detection for single cell genomics. Genome Biol 19, 224 (2018).  
Maksimovic, J. et al. Single-cell atlas of bronchoalveolar lavage from preschool cystic fibrosis reveals new cell phenotypes. bioRxiv 2022.2006.2017.496207 (2022).

## Last updated
04-Dec-2023
