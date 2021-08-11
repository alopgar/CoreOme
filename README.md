[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![Last version](https://img.shields.io/github/tag/alopgar/CoreOme.svg)](https://img.shields.io/github/tag/alopgar/CoreOme.svg)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/CoreOme)](https://cran.r-project.org/package=CoreOme)

# CoreOme
Taxonomy and prevalence filter pipeline for microbiome datasets adapted to SqueezeMeta_Reads output

This pipeline is designed to pre-process SqueezeMeta_reads output by multiple steps:

1) **Samples filter**: Matches phenotypic data (metadata) and microbiome data sample names, listing and removing those with no coincidence in one of both databases. A genotype dataset can be optionally added.
2) **Taxonomy filter**: Used only with taxonomy SQMr files (not KEGG or COG files). Removes those taxa not interesting for post-analysis (normally, animalia, plants, virus and unclassified taxa). At the moment, **this filter must be edited manually**.
3) **Prevalence filter**: This is the main filter. It uses phenotypic data to create a RandomForest model (randomForest package in R) in order to take the microbiome data subcomposition classifying phenotypes with the lowest OOB-error. Data subcompositions are created removing the 5% lowest prevalent features (i.e., those features present in less than 5% of samples), and increasing percentage 5 to 5 to 95% of lowest prevalence (i.e., features present in less than 95% of samples). The subcomposition with the lowest OOB-error will be chosen as core microbiome.
4) **Adding low-prev/high-VI features**: We after perform another RandomForest with the inverse subcomposition of the chosen one, onlly using the discarded features. As an example, if our **core microbiome** contains features present **at least in 40%** of samples (threshold = 40%) our **low-prevalence (LP)** subset will be composed by features present in **less than 40%** of samples, whith which LP-hVI RandomForest is done. We take LP features with a variable importance (VI) higher than 1.0 and a prevalence higher than 5% (can be modified). These LP/hVI features are attached to core microbiome dataset.

## How to use it

a) Use `Rparams.R` to customise paths and other parameters for the pipeline. More info inside the file.  
b) Use `Rphenotmods.R` to modify R classes of phenotypic variables. Useful to convert continuous variables into categorical factors.  
c) Run `SQMr_exe.sh`. It is required to check input, output and bin paths inside the script.
