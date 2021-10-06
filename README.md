#  Hershberger et al. Transcriptome-wide association and prediction for carotenoids and tocochromanols in fresh sweet corn kernels
![language: R](https://img.shields.io/badge/language-R-blue.svg)
[![DOI:10.1101/2021.09.24.461734](http://img.shields.io/badge/DOI-10.1101/2021.09.24.461734-B31B1B.svg)](https://doi.org/10.1101/2021.09.24.461734)

# PREPRINT

This repository and website documents all analyses, summary, tables and figures associated with the following PREPRINT: [Transcriptome-wide association and prediction for carotenoids and tocochromanols in fresh sweet corn kernels](https://doi.org/10.1101/2021.09.24.461734)!

# Abstract

Sweet corn is consistently one of the most highly consumed vegetables in the U.S., providing a valuable opportunity to increase nutrient intake through biofortification. Significant variation for carotenoid (provitamin A, lutein, zeaxanthin) and tocochromanol (vitamin E, antioxidants) levels is present in temperate sweet corn germplasm, yet previous genome-wide association studies (GWAS) of these traits have been limited by low statistical power and mapping resolution. Here, we employed a high-quality transcriptomic dataset collected from fresh sweet corn kernels to conduct transcriptome-wide association studies (TWAS) and transcriptome prediction studies for 39 carotenoid and tocochromanol traits. In agreement with previous GWAS findings, TWAS detected significant associations for four causal genes, _β-carotene hydroxylase_ (_crtRB1_), _lycopene epsilon cyclase_ (_lcyE_), _γ-tocopherol methyltransferase_ (_vte4_), and _homogentisate geranylgeranyltransferase_ (_hggt1_) on a transcriptome-wide level. Pathway-level analysis revealed additional associations for _deoxy-xylulose synthase2_ (_dxs2_), _diphosphocytidyl methyl erythritol synthase2_ (_dmes2_), _cytidine methyl kinase1_ (_cmk1_), and _geranylgeranyl hydrogenase1_ (_ggh1_), of which, _dmes2_, _cmk1_, and _ggh1_ have not previously been identified through maize association studies. Evaluation of prediction models incorporating genome-wide markers and transcriptome-wide abundances revealed a trait-dependent benefit to the inclusion of both genomic and transcriptomic data over solely genomic data, but both transcriptome- and genome-wide datasets outperformed a priori candidate gene-targeted prediction models for most traits. Altogether, this study represents an important step towards understanding the role of regulatory variation in the accumulation of vitamins in fresh sweet corn kernels.

# Data availability
Much of the supporting data and output from the analyses documented here are too large for GitHub.

The `data` and `output` folders for this repository will be mirrored [here](https://datacommons.cyverse.org/browse/iplant/home/shared/GoreLab/dataFromPubs/Hershberger_SweetCornRNA_2021).

# Code

The folder `CODE` contains analysis scripts. Each of the script filenames include numbers that correspond to the order of the analysis. Note that some numbers are skipped (e.g., there is no code starting from "2"). All the scripts were executed the order shown in ExecuteScripts.sh. In total, this folder contains 32 scripts:

1. ExecuteScripts.sh

### Scripts 1.1 to 1.3 are related to the expression data processing before BLUE calculation

- 1.1-LookTagSeq_MakeData.R and 1.1-LookTagSeq_MakeData_Ia453.R curate expression data for each accession but perform no calculations. They remove some control samples for RNAseq that were retained in the raw data files.

- 1.2-GeneFiltering performs gene filtering as described in the manuscript

- 1.3-OutlierRemoval_and_Imputation.R performes outlier removal and imputation based on MAD


### Scripts 3.4 to 3.6 are related to BLUE calculation and outlier removal

- 3.4-BLUE_each_NoHarvDate.R uses asreml-R to calculate expression BLUEs. This code was executed using 3.4-RunMultipleR_B73_NoHarvDate.sh and 3.4-RunMultipleR_Ia453_NoHarvDate.sh for all genes

- 3.5-BLUE_merge_NoHarvDate.R merges the BLUE result files for each gene to the single file of BLUE values

- 3.6-BLUE_RmErr_NoHarvDate.R removes genes for which BLUE calculation did not converge


### Script 4.2 connects expression BLUE data and the phenotypic data files, taking the common subset of accessions etc.

- 4.2-MakeBlueDatasets_NoHarvDate.R takes a common subset of accessions from expression BLUE and phenotypic data files


### Scripts from 5.5 to 5.8 are related to PEER and subsequent outlier removal

- 5.5-Peer_Use25Fact_NoHarvDate.R runs PEER with 25 factors

- 5.6-MakeScreePlot_NoHarvDate.R makes scree plots

- 5.7-Peer_UseOptFact_NoHarvDate.R runs PEER with the given number of factors

- 5.8-OutlierRemoval_NoHarvDate.R performs outlier removal based on Studentized residuals


### Scripts 6.1 and 6.2 are pre-TWAS analyses

- 6.1-CalcKinship_GAPIT calculates genomic relationship matrix (kinship matrix) from the two (11K and 163K) SNP sets

- 6.2-ModelComparison_BIC.R calculates BIC for the model comparison for TWAS.


### Script 7.3 performs TWAS analysis for all vitamin traits.

- All TWAS analyses were performed using 7.3-TWAS_NoHarvDate.R, with the exception of aT using Ia453 expression data. TWAS for this trait was performed using 7.3-TWAS_NoHarvDate_Ia453_aT.R. The "_Ia453_aT" code is equivalent to the 7.3-TWAS_NoHarvDate.R, but there is no for-loop and only that trait was analyzed.


### Scripts 9.6 and 9.7 perform TWAS analysis for the kernel mutant type

- 9.6-ModelComparison_BIC_mutant_type_ReAnalysis.R calculates BIC to choose optimal model

- 9.7-TWAS_Kernel_mutant_type_ReAnalysis.R performs TWAS


### Scripts 13.1 and 13.2 are related to genomic and transcriptome-based prediction

- 13.1-Prediction_stratified.R performs genomic prediction and transcriptome-based prediction with all genes

- 13.2-Prediction_stratified_useCand.R performs transcriptomic prediction using a priori candidate genes


### Scripts 14.1 through 14.6 compile TWAS results and match them with annotations.

- 14.1-TWAS_NAM_JL_QTL.R joins uplifted NAM JL QTL results from Diepenbrock et al 2017, Diepenbrock et al 2021

- 14.2-Match_apriori_with_annotations.R compiles an a priori candidate list and matches it with B73 v4 and Ia453 genome annotations

- 14.3-TWAS_top_hits_apriori_FUNCTIONS.R provides functions for compiling TWAS results and matching with annotations, a priori candidate lists, and NAM JL QTL intervals

- 14.4-match_TWAS_with_annotations_and_apriori.R compiles B73 TWAS results and matches them with annotations, a priori candidate lists, and NAM JL QTL intervals

- 14.5-kernel_TWAS.R compiles kernel type TWAS results

- 14.6-Ia453_TWAS.R sompile TWAS results for Ia453 genome alignment and matches them with annotations and a priori candidate lists


### Script 15.1 compiles prediction results

- 15.1-Prediction_plots.R compiles prediction results and generate a plot (Figure 3)


### Script 16.1 generates a table of kernel endosperm mutation assignments

- 16.1-kernel_mutant_assignments.R assigns kernel mutant types based on line name


---

For questions regarding released code and datasets contact the corresponding authors:

- [Jenna Hershberger, Ph.D.](mailto:jmh579@cornell.edu)

- [Prof. Michael Gore (Cornell University)](mag87@cornell.edu)
