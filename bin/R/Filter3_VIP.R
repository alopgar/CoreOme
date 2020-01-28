options(warn=-1)

suppressMessages(require(tidyverse)); suppressMessages(require(randomForest))

## INIT VARIABLES -------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(T)
source(args[1])
best <- args[3]
# best <- "RanFor_thr05.txt"

write.report <- function(str, file = paste0(WD, "Filtered/R_report.txt")){write(str, file = file, append=T)}
tbl.report <- function(tbl, file = paste0(WD, "Filtered/R_report.txt")){capture.output(tbl, file = file, append=T)}

# DATA LOADING ---------------------------------------------------------------------------------------------------------------------------------
if (datatype == "taxonomy"){
  sgt_mdata <- read.delim(paste0(WD, "Filtered/", "F2.tax_reads.tsv"), stringsAsFactors = F)
} else if (datatype == "keggs"){
  sgt_mdata <- read.delim(paste0(WD, "Filtered/", "F1.mreads.tsv"), stringsAsFactors = F)
}
subset_best <- read.table(paste0(WD, "Filtered/RFSubsets/", best), sep = " ", dec = ".", stringsAsFactors = F, header = T) %>% column_to_rownames("ID")
sg_phenot <- read.csv(paste0(WD, "Filtered/", "F1.metadata.csv"), sep = ";", dec = ",", stringsAsFactors = F, row.names = 1)
sumstats <- read.delim(paste0(WD, "Filtered/", "summary_stats.tsv"), stringsAsFactors = F, row.names = 1)
OOB_table <- read.delim(paste0(WD, "Filtered/", "OOB_table.tsv"), stringsAsFactors = F, row.names = 1)

# METADATA MODIFICATION ------------------------------------------------------------------------------------------------------------------------
source(args[2])

# FILTER 3: PREVALENCE FILTERING: L-P RANDOM FOREST --------------------------------------------------------------------------------------------
write.report(c("\n", "--- IMPORTANCE OF LOW-PREVALENCE FEATURES: L-P RANDOM FOREST", ""))
# Lowest OOB error subset:
write.report(c("- Lowest OOB error subset:", ""))
OOB_min <- OOB_table[OOB_table$OOB_error == min(OOB_table[-1,"OOB_error"]),]
tbl.report(OOB_min)

## DATA SUBSET PREPARATION
# Edit counts table rownames to match easily with tax table.
# Also store apart taxonomy/kegg columns from counts table.
write.report(c("", "- Preparing data subsets for Random Forest..."))
if (datatype == "taxonomy"){
  rownames(sgt_mdata) <- paste(substr(dataselect, 1, 2), rownames(sgt_mdata), sep = "_")
  sgtp_mdata <- as.data.frame(t(sgt_mdata[,-which(names(sgt_mdata) %in% taxdenom[taxrng])]))
  sgtp_taxa <- sgt_mdata[,which(names(sgt_mdata) %in% taxdenom[taxrng])]
} else if (datatype == "keggs"){
  sgt_mdata <- sgt_mdata %>% column_to_rownames("KEGG")
  sgtp_mdata <- as.data.frame(t(sgt_mdata[,-which(names(sgt_mdata) %in% keggvars)]))
  sgtp_func <- sgt_mdata[,which(names(sgt_mdata) %in% keggvars)]
}
# Load variables: phenotipic data & non-zero counts:
pheno_good <- c(RanForFactor, RanForVars)
non_zero_counts <- apply(sgtp_mdata, 2, function(c)sum(c!=0))

rfdat_bad <- sgtp_mdata[,setdiff(colnames(sgtp_mdata), colnames(subset_best))]
rfmdat_bad <- merge(sg_phenot[,(colnames(sg_phenot) %in% pheno_good)], rfdat_bad, by = "row.names", all = T) %>% rename(ID=Row.names)
ranfog_order <- c(RanForFactor, "ID", RanForVars, colnames(rfdat_bad))
rfmdat_bad <- rfmdat_bad[,ranfog_order] %>% column_to_rownames("ID")

tbl.report(rfmdat_bad[1:8,1:10])
write.report(paste("Low prevalence features:", dim(rfmdat_bad)[1], "samples and", dim(rfmdat_bad)[2]-length(pheno_good), "features"))

## IMPORTANCE OF LOW-PREVALENCE FEATURES
# Random Forest 2: using only features with prevalence < min_preval (with high zero number).
# Calculate increment in OOB error and keep features with highest importance.
# Use dim[2]/3 to calculate tuning mtry.
pre_mtry <- round(dim(rfmdat_bad)[2]/3)
write.report(c("", "- Tuning Random Forest: mtry optimum.",
               "  Parameters:",
               paste0("   mtry = ", pre_mtry),
               paste0("   ntree = ", trf_ntree),
               paste0("   nodesize = ", trf_nodesize),
               paste0("   StepFactor = ", trf_SFactor),
               paste0("   improve = ", trf_improve)))
set.seed(12)
trf_bad <- tuneRF(x = rfmdat_bad[, !names(rfmdat_bad) %in% RanForFactor], y = rfmdat_bad[,RanForFactor], ntree = 100, 
                  mtry = pre_mtry, nodesize = 3, StepFactor = 100, improve = 0.05, trace = T, plot = F)
tune_mtry <- as.numeric(names(which(trf_bad[,"OOBError"] == min(trf_bad[,"OOBError"]))))

write.report(c("", "- Applying RandomForest...",
               "  Parameters:",
               paste0("   mtry = ", tune_mtry),
               paste0("   ntree = ", rf_ntree),
               paste0("   nodesize = ", rf_nodesize)))
rfformula <- as.formula(paste(RanForFactor, paste("~ .", paste(RanForVars, collapse = " + "), sep = " + ")))
set.seed(14)
rfa_bad <- randomForest(rfformula, data = rfmdat_bad, mtry = tune_mtry, ntree = 200, nodesize = 3, importance = T)

write.report(c("", "- Top VIP values:"))
vImp_bad <- rfa_bad$importance
vImp_bad <- cbind(vImp_bad, rfa_bad$importanceSD)
colnames(vImp_bad)[3] <- "%IncMSE_SD"
vImp_bad[,"%IncMSE"] <- vImp_bad[,"%IncMSE"]/vImp_bad[,"%IncMSE_SD"]
vImp_bad[is.nan(vImp_bad)] <- 0
tbl.report(head(vImp_bad))

# FILTER 3: PREVALENCE FILTERING: SELECTING VARIABLES BY IMPORTANCE ----------------------------------------------------------------------------
write.report(c("\n", "--- IMPORTANCE OF LOW-PREVALENCE FEATURES: SUBSETTING", ""))
## Subset L-P features that increase %MSE at least in 1% (higher importance)
vImp_bad_top <- sort(vImp_bad[vImp_bad[,"%IncMSE"] >= 1, "%IncMSE"], decreasing = T)
vImp_bad_top <- vImp_bad_top[!(names(vImp_bad_top) %in% pheno_good)]
write.report(c("- Subset L-P features that increase %MSE at least in 1% (higher importance):",
               paste("  Number of features:", length(vImp_bad_top))))

## Prevalence, nr of single/doubletons and average RA, from the vImp_bad_top genera.
# Make linear model: CH4 ~ NLACTA + DIASLE + RA (per sample). Add p-value to the results table.
write.report(c("- General stats and GLM for high-importance L-P variables:",
               "  Prevalence",
               "  Single/doubletons",
               "  Average RA",
               paste("  P-value for gaussian GLM:", 
               paste(RanForFactor, paste(paste(RanForVars, collapse = " + "), "Feature", sep = " + "), sep = " ~ ")), 
               ""))
ra_all <- as.data.frame(apply(sgtp_mdata, 1, function(x){x/sum(x)}))
all_meanRA <- apply(ra_all, 1, mean)
all_prev <- non_zero_counts
all_single <- apply(sgtp_mdata, 2, function(c)sum(c==1))
all_double <- apply(sgtp_mdata, 2, function(c)sum(c==2))

## Check if rownames are equivalent for posterior cbind:
#all.equal(rownames(ra_all), names(all_meanRA))
#all.equal(rownames(ra_all), names(all_prev))
#all.equal(rownames(ra_all), names(all_single))
#all.equal(rownames(ra_all), names(all_double))

# Cbind RA, mean RA, prevalence and single/doubletons of full dataset.
ra_all <- as.data.frame(cbind("Mean_RA" = all_meanRA, "Prevalence" = all_prev, "Singletons" = all_single, "Doubletons" = all_double, ra_all))
# Subset ra_all table taking vImp_bad_top genera.
ra_1perc <- ra_all[rownames(ra_all) %in% names(vImp_bad_top),]

tbl.report(ra_all[1:8,1:6])
write.report(paste("Low prevalence features:", dim(ra_all)[1], "features within", dim(ra_all)[2]-4, "samples"))

## GLM: Table with samples at rows and RA of vImp_bad_top genera at columns. Add phenotype info.
write.report(c("", paste("- GLM:", paste(RanForFactor, paste(paste(RanForVars, collapse = " + "), "Feature", sep = " + "), sep = " ~ ")), ""))
glm_ra <- t(ra_1perc[,-c(1:4)])
glm_ra <- merge(sg_phenot[,(colnames(sg_phenot) %in% pheno_good)], glm_ra, by = "row.names") %>% column_to_rownames("Row.names")
# GLM: one glm per feature (CH4 ~ NLACTA + DIASLE + feat[i])
forinit <- length(pheno_good)
imp_glm <- list(); pval <- NULL
for(i in (forinit+1):ncol(glm_ra)){
  glm_formula <- as.formula(paste(RanForFactor, paste(paste(RanForVars, collapse = " + "), colnames(glm_ra)[i], sep = " + "), sep = " ~ "))
  imp_glm[[i-forinit]] <- glm(glm_formula, family = gaussian, data = glm_ra)
  names(imp_glm)[i-forinit] <- colnames(glm_ra)[i]
  # Genus p-value from coef(summary(glm))
  pval[i-forinit] <- coef(summary(imp_glm[[i-forinit]]))[colnames(glm_ra)[i],"Pr(>|t|)"]
  names(pval)[i-forinit] <- colnames(glm_ra)[i]
}

glm_ra_pv <- merge(data.frame(pval), ra_1perc, by = "row.names") %>% column_to_rownames("Row.names")
if (datatype == "taxonomy"){
  glm_ra_pv_tx <- merge(sgtp_taxa, glm_ra_pv, by = "row.names") %>% column_to_rownames("Row.names")
} else if (datatype == "keggs"){
  glm_ra_pv_fn <- merge(glm_ra_pv, sgtp_func, by = "row.names") %>% column_to_rownames("Row.names")
}

# SAVE FILES -----------------------------------------------------------------------------------------------------------------------------------
write.report(c("\n", "--- Saving files:"))
## Preparing files:
# Join core table with taxonomy table:
if (datatype == "taxonomy"){
  write.report(c("", "- Join H-P filtered counts table with taxonomy table"))
  pre_good <- merge(sgtp_taxa, t(subset_best[,!(colnames(subset_best) %in% pheno_good)]), by = "row.names") %>% column_to_rownames("Row.names")
} else if (datatype == "keggs"){
  write.report(c("", "- Join H-P filtered counts table with functionality table"))
  pre_good <- merge(t(subset_best[,!(colnames(subset_best) %in% pheno_good)]), sgtp_func, by = "row.names") %>% column_to_rownames("Row.names")
}

# Create non-core table using table with GLM+mean_RA+prev+single/doubletons, non-core counts and taxonomy tables:
write.report("- Create table with all info from L-P features with RF importance > 1%")
list_ra_glm <- glm_ra_pv[,1:5] %>% rownames_to_column("Feat")
list_rfdat_bad <- as.data.frame(t(rfdat_bad)) %>% rownames_to_column("Feat")
if (datatype == "taxonomy"){
  list_sgtp_taxa <- sgtp_taxa %>% rownames_to_column("Feat")
  imp_list <- list(list_sgtp_taxa, list_ra_glm, list_rfdat_bad)
} else if (datatype == "keggs"){
  list_sgtp_func <- sgtp_func %>% rownames_to_column("Feat")
  imp_list <- list(list_ra_glm, list_rfdat_bad, list_sgtp_func)
}
sgtp_imp <- reduce(imp_list, inner_join, by = "Feat")

# L-P importance table filtering: Subset features with prevalence > prevfilter value (p-value < 0.05 optional) 
write.report(paste("- L-P importance table filtering: Subset features with imp > 1% and prevalence >", prevfilter))
if(isTRUE(glmfilter)){
  sgtp_bad_prev <- glm_ra_pv[glm_ra_pv$Prevalence > prevfilter,]
  subset.exp <- expression(sgtp_imp$Prevalence > prevfilter & sgtp_imp$pval <= 0.05)
} else if(isFALSE(glmfilter)){
  subset.exp <- expression(sgtp_imp$Prevalence > prevfilter)
}
sgtp_bad_good <- sgtp_imp[eval(subset.exp),] %>% remove_rownames() %>% column_to_rownames("Feat")
sgtp_bad_bad <- sgtp_imp[!(sgtp_imp$Feat %in% row.names(sgtp_bad_good)),] %>% remove_rownames() %>% column_to_rownames("Feat")

# Add those feats to core table:
write.report(c("- Adding filtered L-P features to final H-P filtered counts table"))
add_gen <- select(sgtp_bad_good, -one_of(colnames(list_ra_glm)))
add_gen <- add_gen[,colnames(pre_good)]
#colnames(pre_good) == colnames(add_gen)
final_good <- rbind(pre_good, add_gen)

# L-P feature counts by filters:
write.report(c("", "  L-P feature counts with consecutive filter steps:"))
if(isTRUE(glmfilter)){
  lpstats <- c(dim(ra_all)[1], dim(glm_ra_pv)[1], dim(sgtp_bad_prev)[1], dim(sgtp_bad_good)[1])
  names(lpstats) <- c(paste0("< ", OOB_min$Threshold, "% Prevalence"), "Importance > 1%", paste("Prevalence >", prevfilter), "LM p-value < 0.05")
} else if(isFALSE(glmfilter)){
  lpstats <- c(dim(ra_all)[1], dim(glm_ra_pv)[1], dim(sgtp_bad_good)[1])
  names(lpstats) <- c(paste0("< ", OOB_min$Threshold, "% Prevalence"), "Importance > 1%", paste("Prevalence >", prevfilter))
}
tbl.report(data.frame("Features" = lpstats))

## Filtered summary stats:
if (datatype == "taxonomy"){
  # Step-by-step datasets without taxonomy:
  notax_good <- pre_good[,-taxrng]
  notax_end <- final_good[,-taxrng]
  notaxname <- paste0("Reads_thr", OOB_min$Threshold)
  # Step-by-step taxonomies without reads:
  noreads_good <- pre_good[,taxrng]
  noreads_end <- final_good[,taxrng]
  noreadsname <- paste0("Taxa_thr", OOB_min$Threshold)
  # Merge:
  summary_list <- list(list(notax_good, noreads_good),
                       list("Reads_Final" = notax_end, "Taxa_Final" = noreads_end))
  names(summary_list[[1]]) <- c(notaxname, noreadsname)
  # Final summary:
  for( i in seq(summary_list)){
  srow <- c("nsamples"=ncol(summary_list[[i]][[1]]), 
            "ntaxa"=nrow(summary_list[[i]][[1]]), 
            "nreads"=sum(summary_list[[i]][[1]]), 
            apply(summary_list[[i]][[2]], 2, function(x){length(unique(x))}))
  sumstats <- rbind(sumstats, srow)
  row.names(sumstats)[nrow(sumstats)] <- names(sapply(summary_list, "[", 1))[i]
	}
} else if (datatype == "keggs"){
  # Step-by-step datasets without keggs:
  nokegg_good <- select(pre_good, -one_of(c("Total",keggvars)))
  nokegg_end <- select(final_good, -one_of(c("Total",keggvars)))
  nokeggname <- paste0("Reads_thr", OOB_min$Threshold)
  # Step-by-step keggs without reads:
  noreads_good <- select(pre_good, one_of(c("Total",keggvars)))
  noreads_end <- select(final_good, one_of(c("Total",keggvars)))
  noreadsname <- paste0("KEGG_thr", OOB_min$Threshold)
  # Merge:
  summary_list <- list(list(nokegg_good, noreads_good),
                       list("Reads_Final" = nokegg_end, "KEGG_Final" = noreads_end))
  names(summary_list[[1]]) <- c(nokeggname, noreadsname)
  # Final summary:
  for( i in seq(summary_list)){
    srow <- c("nsamples"=ncol(summary_list[[i]][[1]]), 
              "nkeggs"=nrow(summary_list[[i]][[1]]), 
              "nreads"=sum(summary_list[[i]][[1]]), 
              apply(summary_list[[i]][[2]], 2, function(x){length(unique(x))}))
    sumstats <- rbind(sumstats, srow)
    row.names(sumstats)[nrow(sumstats)] <- names(sapply(summary_list, "[", 1))[i]
	}
}

## Save files:
if (datatype == "taxonomy"){
  write.csv(final_good %>% rownames_to_column("Feat"), paste0(WD, "Filtered/", "SQMr_filtered_genera.csv"), quote = F, row.names = F)
  write.csv(sgtp_imp, paste0(WD, "Filtered/", "SQMr_important_Lprev_genera.csv"), quote = F, row.names = F)
} else if (datatype == "keggs"){
  write.csv(final_good %>% rownames_to_column("Feat"), paste0(WD, "Filtered/", "SQMr_filtered_kegg.csv"), quote = F, row.names = F)
  write.csv(sgtp_imp, paste0(WD, "Filtered/", "SQMr_important_Lprev_kegg.csv"), quote = F, row.names = F)
}
write.table(sumstats, paste0(WD, "Filtered/", "summary_stats.tsv"), quote=F, sep='\t', col.names = NA)

write.report(c("\n", "--- L-P IMPORTANCE FILTER ENDED ---------------------------------------------"))

write.report(c("\n", "--- Filtering summary: Final file:", paste0("Saving path: ", WD, "Filtered/")))
tbl.report(sumstats)
