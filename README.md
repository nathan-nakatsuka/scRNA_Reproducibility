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

**Preprocessing steps:**
1) Perform QC and determine cell type of each cell (e.g. by mapping to an atlas).<br/>
2) Perform differential expression on each cell type in each dataset.:<br/>
Alternatives: <br/>
a) Perform pseudobulking per individual and per cell type if using a bulk RNA-sequencing differential expression method (e.g. DESeq2). 
Alternatively, b) use a mixed model differential expression method to control for within individual correlation of gene expression.
<br/>
Combine the results of the datasets with SumRank.
Step 1) Run SumRank on your real data.

Step 2) Calibrate p-values empirically.
2a) Perform permutations of case-control status and do differential expression on each of your datasets. Then run SumRank on the results of these differential expression analyses as in Step 1.
-Note: The number of permutations to perform will depend on the user's desire of p-value accuracy (this will only determine the cutoff for which genes are considered significant from the results of Step 1). Usually 500-1,000 permutations are enough to get stable p-values.



2b) Concatenate the p-values from Step 2a and compare them to the p-values from Step 1 to obtain calibrated p-values.


<br/>
<br/>
<br/>


