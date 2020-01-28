options(warn=-1)

suppressMessages(require(tidyverse)); suppressMessages(require(randomForest))

## INIT VARIABLES -------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(T)
source(args[1])

# DATA LOADING ---------------------------------------------------------------------------------------------------------------------------------
subfile <- args[2]
#subfile <- paste0(WD, "Filtered/RFSubsets/RanFor_thr95.txt")
cat("\n", "--- Loading data:", subfile, "\n")
ranfordat <- read.table(subfile, sep = " ", dec = ".", stringsAsFactors = F, header = TRUE) %>% column_to_rownames("ID")
sg_phenot <- read.csv(paste0(WD, "Filtered/", "F1.metadata.csv"), sep = ";", dec = ",", stringsAsFactors = F, row.names = 1)

# METADATA MODIFICATION ------------------------------------------------------------------------------------------------------------------------
for(i in RanForVars){
	sg_phenot[,i] <- as.factor(sg_phenot[,i])
	ranfordat[,i] <- as.factor(ranfordat[,i])
}

# FILTER 3: PREVALENCE FILTERING ---------------------------------------------------------------------------------------------------------------
cat("\n", "--- Applying FILTER 3: RandomForest", "\n")
## A) Do RandomForest.
pheno_good <- c(RanForFactor, RanForVars)
pre_mtry <- round(dim(ranfordat)[2]/3)

cat("\n", "- Tuning Random Forest: mtry optimum.", "\n",
    "  Parameters:", "\n",
    "   mtry = ", pre_mtry, "\n",
    "   ntree = ", trf_ntree, "\n",
    "   nodesize = ", trf_nodesize, "\n",
    "   StepFactor = ", trf_SFactor, "\n",
    "   improve = ", trf_improve, "\n\n", sep = "")

set.seed(10)
tunerf <- tuneRF(x = ranfordat[, !names(ranfordat) %in% pheno_good] , y = ranfordat[,RanForFactor], ntreeTry = trf_ntree, mtryStart = pre_mtry, 
                 nodesize = trf_nodesize, stepFactor = trf_SFactor, improve = trf_improve, trace = T, plot = F)
tune_mtry <- as.numeric(names(which(tunerf[,"OOBError"] == min(tunerf[,"OOBError"]))))

cat("\n", "- Running Random Forest:", "\n",
    "  Parameters:", "\n",
    "   mtry = ", tune_mtry, "\n",
    "   ntree = ", rf_ntree, "\n",
    "   nodesize = ", rf_nodesize, "\n\n", sep = "")

rfformula <- as.formula(paste(RanForFactor, paste("~ .", paste(RanForVars, collapse = " + "), sep = " + ")))
#rfformula <- RanForFactor %>% paste("~ .") %>% paste(paste(RanForVars, collapse = " + "), sep = " + ") %>% as.formula()
set.seed(10)
rfmodel <- randomForest(rfformula, data = ranfordat, mtry = tune_mtry, ntree = rf_ntree, nodesize = rf_nodesize, importance = T, plot = T)
OOB_err <- tail(rfmodel$mse, 1)

## B) Prepare output file.
rfout <- capture.output(cat("Random Forest Output: ", subfile, "\n",
                            "---------------------\n\n",
                            "Call:\t", "randomForest(formula = ", deparse(rfformula), ", data = ranfordat, mtry = ", tune_mtry, ", ntree = ", rf_ntree, 
                            ", nodesize = ", rf_nodesize, ", importance = T)\n",
                            "Type of random forest:\t", rfmodel$type, "\n",
                            "Number of trees:\t", rf_ntree, "\n",
                            "No. of variables tried at each split:\t", tune_mtry, "\n\n",
                            "Mean of squared residuals:\t", OOB_err, "\n",
                            "% Var explained:\t", tail(rfmodel$rsq, 1)*100, "\n", sep = ""))

# SAVE FILES -----------------------------------------------------------------------------------------------------------------------------------
cat("\n", "--- Saving OOB table...", "\n")
write.table(rfout, paste0(str_replace(subfile, ".txt", "_out"), ".txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)

cat("\n", "--- PREVALENCE FILTER ENDED ---------------------------------------------", "\n")
