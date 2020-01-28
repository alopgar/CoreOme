# MicrobFilter
Taxonomy and prevalence filter pipeline for microbiome datasets adapted to SqueezeMeta_Reads output

This pipeline is designed to pre-process SqueezeMeta_reads output by multiple steps:

1) Samples filter: Matches phenotypic data (metadata) and microbiome data sample names, listing and removing those with no coincidence in one of both databases. A genotype dataset can be optionally added.
2) Taxonomy filter: Used only with taxonomy SQMr files (not kegg COG files). Removes those taxa not interesting for post-analysis (normally, animalia, plants, virus and unclassified taxa). At the moment, **this filter must be edited manually**.
3) Prevalence filter: This is the main filter. It uses phenotypic data to create a RandomForest model in order to take the microbiome data subcomposition classifying phenotypes with the lowest OOB-error. Data subcompositions are created removing the lowest prevalent taxa, (...)

## How to use it
