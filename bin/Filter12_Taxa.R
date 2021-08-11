options(warn=-1)

suppressMessages(library(tidyverse))

## INIT VARIABLES -------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(T)
source(args[1])

write.report <- function(str, file = paste0(WD, "Filtered/R_report.txt")){write(str, file = file, append=T)}
tbl.report <- function(tbl, file = paste0(WD, "Filtered/R_report.txt")){capture.output(tbl, file = file, append=T)}

# Regex pattern present in taxonomic ref. db, for removing it. Possible values: "\\[\\w+\\] ", "\\w_\\d__", "\\w_"
taxpattern <- "\\w_" 

# DATA LOADING ---------------------------------------------------------------------------------------------------------------------------------
sg_mdata <- read.delim(paste0(WD, "Filtered/", "F1.mreads.tsv"), stringsAsFactors = F)
sumstats <- read.delim(paste0(WD, "Filtered/", "summary_stats.tsv"), stringsAsFactors = F, row.names = 1)

# FILTER 2: TAXONOMY FILTERING -----------------------------------------------------------------------------------------------------------------
write.report(c("\n", "--- Applying FILTER 2: Taxonomy", ""))
dataselect_col <- which(names(sg_mdata) %in% dataselect)
sgt_mdata <- sg_mdata

# 0) Remove prefix (pattern) from taxonomy:
sgt_mdata[1:dataselect_col] <- apply(sgt_mdata[1:dataselect_col], 2, function(x){gsub(taxpattern, "", x)})

# 1) Remove non-mapped reads (sk = "Unclassified") and Viruses:
write.report("TAX CONTROL I: Total superkingdoms present:")
tbl.report(unique(sgt_mdata$superkingdom))
write.report(c("", "TAX CONTROL I: Filtering Virus/Unclassified..."))
sgt_mdata <- sgt_mdata[sgt_mdata$superkingdom != "Viruses" & sgt_mdata$superkingdom != "Unclassified",]

# 2) Identify animal and plant phyla and remove them:
# First, identify animal/plant phyla.
write.report(c("", "TAX CONTROL II: Total Eukaryotic phyla present:"))
tbl.report(unique(sgt_mdata[sgt_mdata$superkingdom == "Eukaryota" & str_detect(sgt_mdata$phylum, "(no phylum in NCBI)") == FALSE,"phylum"]))
nope <- c("Acanthocephala", "Annelida", "Arthropoda", "Brachiopoda", "Bryozoa", "Chaetognatha", "Chlorophyta", "Chordata", "Cnidaria", "Ctenophora", 
		  "Dicyemida", "Echinodermata", "Entoprocta", "Gnathostomulida", "Hemichordata", "Mollusca", "Nematoda", "Nematomorpha", "Nemertea", "Onychophora", 
		  "Orthonectida", "Placozoa", "Platyhelminthes", "Porifera", "Rhodophyta", "Rotifera", "Streptophyta", "Tardigrada", "Xenacoelomorpha")
# Second, identify animal/plant phyla with no phylum NCBI record.
write.report(c("", "TAX CONTROL II: Animal/plant phyla with no NCBI record:"))
nope_ncbi <- unique(sgt_mdata[sgt_mdata$superkingdom == "Eukaryota" & str_detect(sgt_mdata$phylum, "(no phylum in NCBI)") == TRUE,"phylum", drop = F])
tbl.report(nope_ncbi)
nope_ncbi <- paste(c("Dicyemida", "Priapulimorpha", "Rhopaluridae"), "(no phylum in NCBI)")
# Finally, filter animal/plant genera:
write.report(c("", "TAX CONTROL II: Filtering Animal/Plant..."))
sgt_mdata <- sgt_mdata[!sgt_mdata$phylum %in% c(nope, nope_ncbi),]

# 3) Remove Unclassified Family taxa (and consequently, all superior level unclassifieds):
write.report(c("", "TAX CONTROL III: Filtering Unclassified families..."))
sgt_mdata <- sgt_mdata[str_detect(sgt_mdata$family, "Unclassified") == F,]

# 4) Remove no-NCBI Family + Genus taxa (both conditions are required):
write.report(c("", "TAX CONTROL IV: Filtering no-NCBI family + genus..."))
sgt_mdata <- subset(sgt_mdata, (str_detect(sgt_mdata$family, "(no family in NCBI)") & 
                                  str_detect(sgt_mdata$genus, "(no genus in NCBI)")) == F)

# 5) Replace patterns: final "bacterium/archaeon (no genus in NCBI)" by initial "Unclassified".
# Also, remaining "No genus" will be replaced by "Unclassified family".
write.report(c("", "TAX CONTROL V: Replacing patterns to 'unclassified'..."))
for(t in seq(sgt_mdata$genus)){
  if(str_detect(sgt_mdata[t,"genus"], " bacterium \\(no genus in NCBI\\)")){
    sgt_mdata[t,"genus"] <- gsub(" bacterium \\(no genus in NCBI\\)", "", paste("Unclassified", sgt_mdata[t,"genus"]))
  }
  else if(str_detect(sgt_mdata[t,"genus"], " archaeon \\(no genus in NCBI\\)")){
    sgt_mdata[t,"genus"] <- gsub(" archaeon \\(no genus in NCBI\\)", "", paste("Unclassified", sgt_mdata[t,"genus"]))
  }
  else if(str_detect(sgt_mdata[t,"genus"], " \\(no genus in NCBI\\)")){
    sgt_mdata[t,"genus"] <- gsub(".+ \\(no genus in NCBI\\)", paste("Unclassified", sgt_mdata[t,"family"]), sgt_mdata[t,"genus"])
  }
}

# 6) Group duplicated genera in one row (ones that remain duplicated are named as unique genera):
write.report(c("", "TAX CONTROL VI: Fixing duplicated rows..."))
aggregate.form <- as.formula(paste(". ~", paste(taxdenom[1:which(taxdenom == dataselect)], collapse = " + ")))
sgt_mdata <- aggregate(aggregate.form, sgt_mdata, sum)
sgt_mdata[,dataselect] <- with(sgt_mdata, ave(eval(parse(text = dataselect)), FUN = make.unique))

# FINAL STATS ----------------------------------------------------------------------------------------------------------------------------------
## Filtered stats:
# Step-by-step datasets without taxonomy:
notax_sgt <- sgt_mdata[,-c(1:which(taxdenom %in% dataselect))]
# Step-by-step taxonomies without reads:
noreads_sgt <- sgt_mdata[,c(1:which(taxdenom %in% dataselect))]
# Merge:
summary_list <- list(list("Reads_filt2" = notax_sgt, "Taxa_filt2" = noreads_sgt))

# Taxonomy summary:
for( i in seq(summary_list)){
  srow <- c("nsamples"=ncol(summary_list[[i]][[1]]), 
            "ntaxa"=nrow(summary_list[[i]][[1]]), 
			"nreads"=sum(summary_list[[i]][[1]]), 
            apply(summary_list[[i]][[2]], 2, function(x){length(unique(x))}))
  sumstats <- rbind(sumstats, srow)
  row.names(sumstats)[nrow(sumstats)] <- names(sapply(summary_list, "[", 1))[i]
}

write.report(c("\n", "--- Filtering summary: Part 2 (with taxonomy):"))
tbl.report(sumstats)

# SAVE FILTERED TABLES -------------------------------------------------------------------------------------------------------------------------
write.report(c("\n", "--- Saving files..."))
write.table(sgt_mdata, paste0(WD, "Filtered/", "F2.tax_reads.tsv"), quote = F, sep='\t', row.names=FALSE)
write.table(sumstats, paste0(WD, "Filtered/", "summary_stats.tsv"), quote=F, sep='\t', col.names=NA)

write.report(c("\n", "--- TAXA FILTER ENDED ---------------------------------------------------"))
