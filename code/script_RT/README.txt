README (written by Ryokei Tanaka)

This is a folder to sahre my R codes for other project members. 
Note that these scripts may not be the final version.
I may put R codes under development, or I may change the R codes in the future.

My R codes assume that there is a RAWDATA folder in your working directory, where all the raw data files are stored.
Currently (Dec 10th, 2020), RAWDATA folder has the same subfolder structure as the shared google drive. 
For example, raw expression data is located under "RAWDATA/Seetcorn_TagSeq".

# Update Note: 2/8/2021
Many codes afrer the BLUE calculation are updated. Below is the list of changes:
[1] Genes were removed if more than 90% samples have zero-rlog2 expression proflie (in other words, at least 10% of samples have non-zero rlog2 expression profile for the retained genes)
[2] Outlier removal step was inserted after the PEER calcualtion. An "intercept-only" model is fitted to each gene. We remove a PEER residual of the sample for the gene, if its (absolute value of the) studentized residual is greater than the Bonferroni adjusted P = 0.05 threshold. This is same as the Plant Cell paper (Diepenbrock et al. 2017)
[3] We use BIC to find an optimal model for TWAS. We test whether we should include kernel mutant type, kinship matrix, and PCs from the kinship matrix
[4] For the genomic prediction, we use GRM calculated from all SNPs (for TWAS, we are using kinship matrix calculated from a subset of SNPs)


