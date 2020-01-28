options(warn=-1)

suppressMessages(require(tidyverse)); suppressMessages(require(doParallel))

## INIT VARIABLES -------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(T)
source(args[1])

write.report <- function(str, file = paste0(WD, "Filtered/R_report.txt")){write(str, file = file, append=T)}
tbl.report <- function(tbl, file = paste0(WD, "Filtered/R_report.txt")){capture.output(tbl, file = file, append=T)}

# DATA LOADING ---------------------------------------------------------------------------------------------------------------------------------
if (datatype == "taxonomy"){
  sgt_mdata <- read.delim(paste0(WD, "Filtered/", "F2.tax_reads.tsv"), stringsAsFactors = F)
} else if (datatype == "keggs"){
  sgt_mdata <- read.delim(paste0(WD, "Filtered/", "F1.mreads.tsv"), stringsAsFactors = F)
}
sg_phenot <- read.csv(paste0(WD, "Filtered/", "F1.metadata.csv"), sep = ";", dec = ",", stringsAsFactors = F, row.names = 1)
sumstats <- read.delim(paste0(WD, "Filtered/", "summary_stats.tsv"), stringsAsFactors = F, row.names = 1)

# METADATA MODIFICATION ------------------------------------------------------------------------------------------------------------------------
source(args[2])

# FILTER 3: PREVALENCE FILTERING ---------------------------------------------------------------------------------------------------------------
write.report(c("\n", "--- Applying FILTER 3: Subset preparation for RandomForest", ""))
## A) Edit counts table rownames (if datatype = taxonomy, do it to match easily with tax table).
if (datatype == "taxonomy"){
  rownames(sgt_mdata) <- paste(substr(dataselect, 1, 2), rownames(sgt_mdata), sep = "_")
  sgtp_mdata <- sgt_mdata[,-which(names(sgt_mdata) %in% taxdenom[1:which(taxdenom == dataselect)])]
  # Store only taxonomy columns from counts table.
  sgtp_taxa <- sgt_mdata[,which(names(sgt_mdata) %in% taxdenom[1:which(taxdenom == dataselect)])]
  # Select randomforest dataset:
  rfdataset <- sgtp_mdata
} else if (datatype == "keggs"){
  sgt_mdata <- remove_rownames(sgt_mdata) %>% column_to_rownames("KEGG")
  sgp_kdata <- select(sgt_mdata, -c(Function, Class))
  sgp_kpaths <- select(sgt_mdata, Function, Class)
  # Select randomforest dataset:
  rfdataset <- sgp_kdata
}

## B) RandomForest table (transpose & add metadata):
rf_data <- t(rfdataset)
# Select only phenotypic data to use
pheno_good <- c(RanForFactor, RanForVars)
# Merge phenotype and data:
rfm_data <- merge(select(sg_phenot, one_of(pheno_good)), rf_data, by = "row.names", all = T) %>% column_to_rownames("Row.names")

## C) Prepare prevalence subsets:
write.report("- Counting non-zeros per feature:")
non_zero_counts <- apply(rf_data, 2, function(c)sum(c!=0))

write.report(c("", "- Preparing subsets: 5% increasing prevalence threshold"))
cl <- makeCluster(detectCores()[1] - 2)
registerDoParallel(cl)

rf_msubs <- foreach(i=0:19) %dopar% {
  require(tidyverse)
  min_preval <- 0.05*i
  bad_cols <- which(non_zero_counts < (nrow(rf_data) * min_preval))
  if(purrr::is_empty(bad_cols)){rf_subset <- rf_data}
  else {rf_subset <- rf_data[,-(bad_cols)]}
  rfm_subset <- merge(select(sg_phenot, one_of(pheno_good)), rf_subset, by = "row.names", all = T) %>% rename(ID = Row.names)
  ranfog_order <- c(RanForFactor, "ID", RanForVars, colnames(rf_subset))
  rfm_subset <- rfm_subset[,ranfog_order]
  list("Subset"=rf_subset, "Subset_metadata"=rfm_subset)
}
for(i in 0:19){names(rf_msubs)[i+1] <- paste0("PrevThr_", (0.05*i)*100, "%")}
stopCluster(cl)

## D) Build subset stats table (future OOB_table):
write.report(c("", "- Building subset stats table (future OOB_table)"))
OOB_table <- data.frame(as.numeric(gsub("^.+_|%", "", names(rf_msubs))))
colnames(OOB_table)[1] <- "Threshold"
OOB_table <- mutate(OOB_table, Min_samples = floor(nrow(rf_data)*(Threshold/100)))
row.names(OOB_table) <- names(rf_msubs)
OOB_table$N_features <- as.numeric(lapply(rf_msubs, function(x){dim(x[[1]])[2]}))
OOB_table$N_reads <- as.numeric(lapply(rf_msubs, function(x){sum(x[[1]])}))
tbl.report(OOB_table)

# SAVE FILES -----------------------------------------------------------------------------------------------------------------------------------
write.report(c("\n", "--- Saving files..."))
#write.table(sg_phenot, paste0(WD, "Filtered/", "F3.metadata.csv"), sep=';', dec='.', quote = F, col.names=NA)
write.table(OOB_table, paste0(WD, "Filtered/", "OOB_table.tmp"), quote=FALSE, sep='\t', col.names=NA)
dir.create(file.path(WD, "Filtered/", "RFSubsets"), showWarnings = FALSE)
for(i in 0:19){
  min_preval <- 0.05*i
  write.table(rf_msubs[[i+1]][[2]], paste0(WD, "Filtered/RFSubsets/RanFor_thr", stringr::str_pad(min_preval*100, 2, pad="0"), ".txt"), sep = ' ', 
              quote = F, row.names = F)
}

write.report(c("\n", "--- PREVALENCE FILTER ENDED ---------------------------------------------"))
