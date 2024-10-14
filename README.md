# SumRank

These are the relevant scripts for the SumRank software and additional scripts for the following manuscript:

**Citation:** 
<br/>
Nakatsuka, N.; Adler, D.; Jiang, L.; Hartman, A.; Cheng, E.; Klann, E.; Satija, R. “A Reproducibility Focused Meta-Analysis Method for Single-Cell Transcriptomic Case-Control Studies Uncovers Robust Differentially Expressed Genes.” In revision.

**Contact:** Nathan Nakatsuka: 08nanaka@gmail.com

## <p>Overview:</p>
SumRank is a software that allows the assessment of the reproducibility of genes across multiple datasets. The algorithm works by summing up the relative differential expression ranks of genes in each dataset (while also allowing only a certain percentage of datasets to be used for each gene if desired). Genes with high relative rank in many datasets are considered reproducibly differentially expressed. Permutations are used to calibrate p-values. </p>

## <p>Software Requirements:</p>

**OS Requirements:**

This package is supported for macOS and LInux. The package has been tested on the following systems:</p>
* macOS: Sonoma (14.6.1)
* Linux: Ubuntu (23.6.0)

**R Dependencies:**
</p>
R version 4.3.1 or higher

R libraries:</p>
```
Seurat version 4.3.0 or higher
SeuratDisk
DESeq2
unifed version 1.16
ggplot2
```

## <p>Steps for use:</p>

**Preprocessing steps:**
1) Perform QC and determine cell type of each cell (e.g. by mapping to an atlas). The MappingtoAzimuthReference.R code can be used for this.<br/>
2) Perform differential expression on each cell type in each dataset:<br/>
Alternatives: <br/>
a) Perform pseudobulking per individual and per cell type if using a bulk RNA-sequencing differential expression method (e.g. DESeq2). The PseudoBulking.R code can be used for this. <br/>
Alternatively, b) use a mixed model differential expression method to control for within individual correlation of gene expression.<br/>
-Important note: when doing differential expression, you must output ALL genes (i.e. do NOT impose any threshold cutoffs), because SumRank is dependent on the relative ranks of all genes in all datasets.<br/>
<br/>

**Combine the results of the datasets with SumRank.** 
<br/>
1) Run SumRank on your real data using the SumRank function.<br/>
2) Calibrate p-values empirically.<br/>
2a) Perform permutations of case-control status using the PermuteCaseControl function and do differential expression on each of your datasets. Then run SumRank on the results of these differential expression analyses as in Step 1.<br/>
-Note: The number of permutations to perform will depend on the user's desire of p-value accuracy (this will only determine the cutoff for which genes are considered significant from the results of Step 1). Usually 500-1,000 permutations are enough to get stable p-values.<br/>
2b) Concatenate the p-values from Step 2a and compare them to the p-values from Step 1 using the CalibratePValueswithPermutations function to obtain calibrated p-values.
3) Plot Manhattan plot with the MakeManhattanPlot function.

<br/>
See Vignette.md for example of how to run the code.
<br/>
<br/>

Notes:
<br/>
The AdditionalExamplecode.R has example code for the merge, inverse variance, and sex specific analyses.

The UCell_AUC_code.R document has example code for obtaining UCell scores from gene sets and using them to test case-control statuses in held out datasets.

## <p>License:</p>
This project is covered under the Apache 2.0 License.
