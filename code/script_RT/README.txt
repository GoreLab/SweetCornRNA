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


# ---------- quick notes ---------- #
Q. Why do we use rrBLUP::GWAS instead of GAPIT for the TWAS analysis?
A. Because we need to use continuious explanatory variable (= PEER residuals) in TWAS. Further, we have missing values in our PEER data due to the outlier removal. For this type of data, GAPIT does not work well, but rrBLUP does. In our internal comparison, rrBLUP::GWAS and GAPIT are comparable, meaning that their P-values are almost 100% correlated (though their absolute values are not perfectly the same). Thus, we think it is no problem to use rrBLUP instead of GAPIT for TWAS. Note that I also checked that the principal component calcualtion in rrBLUP is equivalent to the one in GAPIT (both calcualte eigenvector from the kinship matrix)

Q. Why there are two outlier removal codes?
A. One code uses base-R functions for the outlier removal, while the other uses asreml. I had confiemd that they are 100% equivalent (i.e., studentized residuals are same), at least for the B73-aligend dataset. As we are fitting an intercept-only model (very simple!), we do not have to use asreml. The reason why we have asreml-version is simply because we had been using it (together with other calcualtion steps such as BLUE/BLUP model fitting). You can use either one (as they are same).


