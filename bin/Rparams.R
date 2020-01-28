## R parameter file:

## INITIAL PARAMETERS ---------------------------------------------------------------------------------------
# WDIR: directory with SQMreads outputs stored.
# DbPath: Reference database directory.
# metadata: path and name of metadata phenotype file.
# genotypes: path and name of genotypes file. FOr not providing genotypes, use FALSE.
# datatype: Which data table are we filtering. 
#			Possible values: "taxonomy", "keggs".
# taxrank: If we are using taxonomy tables, which taxonomic rank are we analising. 
#		   Possible values: "superkingdom", "phylum", "class", "order", "family", "genus", "species".

WD <- "/mnt/lustre/scratch/home/otras/ini/alg/Results/Metagenome/METALGEN_SQMReads/"
DbPath <- "/mnt/lustre/scratch/home/otras/ini/alg/db/"

metadata <- paste0(WD, "metadata.csv")
genotypes <- FALSE
#genotypes <- paste0(WD, "METALGENgenotypes_25_SET_2019.txt")

datatype <- "keggs"
rfmethod <- "RandomForest"
taxrank <- "genus"

## Other optional parameters to edit...
if (datatype == "taxonomy"){
  dataselect <- taxrank
  # Names assigned for taxa ranges
  taxdenom <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species") 
  # Level for classification limit (max: 7 = species)
  taxlvl <- which(taxdenom %in% dataselect)
  taxrng <- seq(taxlvl)
  taxDAlvls <- taxdenom[1:which(taxdenom %in% dataselect)]
} else if (datatype == "keggs") {
  keggvars <- c("KEGG", "Function", "Class")
}

## RANDOM FOREST PARAMETERS ---------------------------------------------------------------------------------
# RanForFactor: Phenotype variable used for discrimination in RandomForest.
# RanForVars: Phenotype variables used as effects in RandomForest model.
# trf_*: Parameters for mtry tuning.
# rf_*: Parameters for Random Forest.

RanForFactor <- "ch4.spkpm_adjusted"
RanForVars <- c("NLACTA", "DIASLE")

trf_ntree = 100
trf_nodesize = 3
trf_SFactor = 2
trf_improve = 0.05

rf_ntree = 200
rf_nodesize = 3

## VIP FILTER PARAMETERS ------------------------------------------------------------------------------------
prevfilter = 5
glmfilter = FALSE
