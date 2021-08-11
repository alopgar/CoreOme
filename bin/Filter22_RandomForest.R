options(warn=-1)

suppressMessages(require(tidyverse)); suppressMessages(require(randomForest)); suppressMessages(require(doParallel))

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

rfformula <- as.formula(paste(RanForFactor, paste("~ .", paste(RanForVars, collapse = " + "), sep = " + ")))

cl <- makeCluster(detectCores()[1] - 2)
registerDoParallel(cl)
rf_out <- foreach(i=1:rf_iters) %dopar% {
	require(randomForest)
	tunerf <- tuneRF(x = ranfordat[, !names(ranfordat) %in% pheno_good] , y = ranfordat[,RanForFactor], ntreeTry = trf_ntree, mtryStart = pre_mtry, 
					 nodesize = trf_nodesize, stepFactor = trf_SFactor, improve = trf_improve, trace = T, plot = F)
	tune_mtry <- as.numeric(names(which(tunerf[,"OOBError"] == min(tunerf[,"OOBError"]))))
	rfmodel <- randomForest(rfformula, data = ranfordat, mtry = tune_mtry, ntree = rf_ntree, nodesize = rf_nodesize, importance = T, plot = F)
	OOB_err <- tail(rfmodel$mse, 1)
	list("Model" = rfmodel,
		 "Out_table" = data.frame("mtry" = tune_mtry, "ntree" = rf_ntree, "nodesize" = rf_nodesize, "OOBe" = OOB_err))
}
for(i in 1:rf_iters){names(rf_out)[i] <- paste0("Iteration #", i)}
stopCluster(cl)

rf_table <- do.call(rbind.data.frame, lapply(rf_out, '[[', 2))
rf_table["Mean",] <- apply(rf_table, 2, mean)
rf_table

## B) Prepare output file.
rfout <- capture.output(cat("Random Forest Output: ", subfile, "\n",
							"---------------------\n\n",
							"Call:\t", "randomForest(formula = ", deparse(rfformula), ", data = ranfordat, mtry = c(", 
							unique(rf_table[-nrow(rf_table),"mtry"]), "), ntree = ", rf_table["Mean","ntree"], ", nodesize = ", 
							rf_table["Mean","nodesize"], ", importance = T)\n",
							"Type of random forest:\t", rf_out[[1]]$Model$type, "\n",
							"Number of trees:\t", rf_ntree, "\n",
							"No. of variables tried at each split:\t", unique(rf_table[-nrow(rf_table),"mtry"]), "\n\n",
							"Mean of squared residuals:\t", rf_table["Mean","OOBe"], "\n",
							"% Var explained:\t", mean(unlist(lapply(rf_out, function(x){tail(x$Model$rsq, 1)*100}))), "\n", sep = ""))

# SAVE FILES -----------------------------------------------------------------------------------------------------------------------------------
cat("\n", "--- Saving OOB table...", "\n")
write.table(rfout, paste0(str_replace(subfile, ".txt", "_out"), ".txt"), quote=FALSE, row.names = FALSE, col.names = FALSE)

cat("\n", "--- PREVALENCE FILTER ENDED ---------------------------------------------", "\n")
