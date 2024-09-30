# SumRank

These are the relevant scripts for the following manuscript:

**Citation:** 
<br/>
Nakatsuka, N.; Adler, D.; Jiang, L.; Hartman, A.; Cheng, E.; Klann, E.; Satija, R. “A Reproducibility Focused Meta-Analysis Method for Single-Cell Transcriptomic Case-Control Studies Uncovers Robust Differentially Expressed Genes.” In revision.

**Contact:** Nathan Nakatsuka: 08nanaka@gmail.com


Notes:
<br/>
The SumRank_code.R document has code for the SumRank algorithm as well as the sex specific analyses.

The UCell_AUC_code.R document has code for obtaining UCell scores from gene sets and using them to test case-control statuses in held out datasets.



## <p>Steps for use:</p>

Preprocessing steps:
Perform QC and determine cell type of each cell (e.g. by mapping to an atlas).
Perform differential expression on each cell type in each dataset.:
a) Perform pseudobulking per individual and per cell type if using a bulk RNA-sequencing differential expression method (e.g. DESeq2). 
Alternatively, b) use a mixed model differential expression method to control for within individual correlation of gene expression.

Combine the results of the datasets with SumRank.
Step 1) Run SumRank on the raw data.

Step 2) Calibrate p-values empirically. Perform permutations of case-control status on datasets and do differential 


















Vignette Example:
1) Download the following COVID-19 datasets in Seurat object form from: https://atlas.fredhutch.org/fredhutch/covid/
-Wilk, Arunachalam, Lee, Wen

2) Perform differential expression
3) 











   


