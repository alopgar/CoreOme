options(warn=-1)

suppressMessages(library(tidyverse))

## INIT VARIABLES -------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(T)
source(args[1])

write.report <- function(str, file = paste0(WD, "Filtered/R_report.txt")){write(str, file = file, append=T)}
tbl.report <- function(tbl, file = paste0(WD, "Filtered/R_report.txt")){capture.output(tbl, file = file, append=T)}

write(c("SQM reads preprocessing: triple filter pipeline.",
		"-------------------------------------------------------------------------------"),
	  file = paste0(WD, "Filtered/R_report.txt"))

write.report(c("Parameters:",
			   paste0("WD: ", WD),
			   paste0("metadata: ", metadata),
			   paste0("genotypes: ", genotypes),
			   paste0("\n", "datatype: ", datatype),
			   if (datatype == "taxonomy"){paste0("taxrank: ", taxrank)},
			   paste0("\n", "rfmethod: ", rfmethod),
			   paste0("Factor for RF: ", RanForFactor),
			   paste0("Covariables for RF: ", paste(RanForVars, collapse=", "))))

# DATA LOADING ---------------------------------------------------------------------------------------------------------------------------------
write.report(c("\n", "--- Loading data...", ""))
write.report("- Metagenome reads:")
if (datatype == "taxonomy"){
	tabnms <- list.files(paste0(WD, "SQMreads/"), pattern = "combined")
	datlist <- lapply(tabnms, function(x){read.delim(paste0(WD, "SQMreads/", x), stringsAsFactors = F, na.strings = c("", "NA"))})
	names(datlist) <- gsub("combined\\.|\\.allfilter.+|\\.prokfilter.+|\\.abund.+", "", tabnms)
	for(i in seq(datlist)){
		colnames(datlist[[i]]) <- gsub("\\.", "_", colnames(datlist[[i]]))
	}

	# We can rely on the lowest taxonomy table if we want, as it englobes everything above with "Unclassified" tag.
	write.report(paste("Loading dataset:", dataselect, "taxonomic level"))
	a_mdata <- datlist[[dataselect]] %>% separate(X, into = taxdenom, sep = ";")
	oldcols <- colnames(a_mdata)[str_detect(pattern="_old", names(a_mdata))]
	a_mdata <- a_mdata[,!str_detect(pattern="_old", names(a_mdata))]
	a_mdata <- a_mdata[ ,!(names(a_mdata) %in% taxdenom[-(1:which(taxdenom == dataselect))])]
} else if (datatype == "keggs") {
	a_mdata <- read.delim(paste0(WD, "SQMreads/jointkeggtable.txt"), stringsAsFactors = F, na.strings = c("", "NA"))
	a_mdata$KEGG <- str_replace_all(a_mdata$KEGG, ";", "_")
} else {
	write.report("ERROR: datatype must be 'taxonomy' or 'keggs'")
	stop()
}

tbl.report(trunc_mat(a_mdata, n = 5, n_extra = 10))

## PHENOTYPIC DATA: IDs + phenotypes.
write.report(c("", "- Metadata:"))
a_phenot <- read.csv(metadata, sep = ";", dec = ",", stringsAsFactors = F)
names(a_phenot)[1] <- "ID"
a_phenot$ID <- gsub(x = a_phenot$ID, pattern = "-", replacement = "_")
row.names(a_phenot) <- a_phenot$ID

tbl.report(trunc_mat(a_phenot, n = 5, n_extra = 10))

## GENOTYPIC DATA: from +1k animals.
write.report(c("", "- Genotypes:"))
if(isFALSE(genotypes)){
	write.report("Genotypes not provided.")
} else {
	a_genot <- read.table(genotypes, sep = " ", dec = ".", header = F, stringsAsFactors = F)
	tbl.report(trunc_mat(a_genot, n = 5, n_extra = 10))
}

# FILTER 1: SAMPLE FILTERING -------------------------------------------------------------------------------------------------------------------
write.report(c("\n", "--- Applying FILTER 1: Samples filtering", ""))
## a) Remove mdata samples not present in phenot and viceversa:
s_phenot <- a_phenot[which(rownames(a_phenot) %in% colnames(a_mdata)),]
if (datatype == "taxonomy"){
	s_mdata <- a_mdata[,which(colnames(a_mdata) %in% c(taxdenom[taxrng], a_phenot$ID))]
} else if (datatype == "keggs") {
	s_mdata <- a_mdata[,which(colnames(a_mdata) %in% c(keggvars, a_phenot$ID))]
}

## b) Remove samples not present in genotype file:
if(isFALSE(genotypes)){
	sg_phenot <- s_phenot
} else {
	sg_phenot <- s_phenot[which(s_phenot$numero %in% a_genot[,1]),]
	sg_genot <- a_genot[a_genot[,1] %in% sg_phenot$numero,] %>% distinct(V1, .keep_all = TRUE)
}

if (datatype == "taxonomy"){
	sg_mdata <- s_mdata[,which(colnames(s_mdata) %in% c(taxdenom[taxrng], sg_phenot$ID))]
} else if (datatype == "keggs") {
	sg_mdata <- s_mdata[,which(colnames(s_mdata) %in% c(keggvars, sg_phenot$ID))]
}

## c) See unmatched samples:
if (datatype == "taxonomy"){
	write.report("- Old duplicated samples (removed):")
	tbl.report(oldcols)
}

write.report(c("", "- Samples with metagenome data but no phenotipic data:"))
if (datatype == "taxonomy"){
	setdiff(colnames(a_mdata), taxdenom) %>% setdiff(a_phenot$ID) %>% tbl.report()
} else if (datatype == "keggs") {
	setdiff(colnames(a_mdata), c("Total",keggvars)) %>% setdiff(a_phenot$ID) %>% tbl.report()
}

write.report(c("", "- Samples with phenotipic data but no metagenome:"))
tbl.report(setdiff(a_phenot$ID, colnames(a_mdata)))

if(!isFALSE(genotypes)){
	write.report(c("", "- Samples with phenotipic data but no genotype:"))
	filter(s_phenot, numero %in% setdiff(numero, a_genot[,1])) %>% pull(ID) %>% tbl.report()
	
	write.report(c("", "- Samples with metagenome data but no genotype:"))
	if (datatype == "taxonomy"){
		setdiff(colnames(s_mdata), taxdenom) %>% setdiff(sg_phenot$ID) %>% tbl.report()
	} else if (datatype == "keggs") {
		setdiff(colnames(s_mdata), c("Total",keggvars)) %>% setdiff(sg_phenot$ID) %>% tbl.report()
	}
}

## d) Filtered stats:
if (datatype == "taxonomy"){
	# Step-by-step datasets without taxonomy:
	notax_a <- select(a_mdata, -one_of(taxdenom))
	notax_sg <- select(sg_mdata, -one_of(taxdenom))
	# Step-by-step taxonomies without reads:
	noreads_a <- select(a_mdata, one_of(taxdenom))
	noreads_sg <- select(sg_mdata, one_of(taxdenom))
	# Merge:
	summary_list <- list(list("Reads_brute" = notax_a, "Taxa_brute" = noreads_a), 
						list("Reads_filt1" = notax_sg, "Taxa_filt1" = noreads_sg))
	# Taxonomy summary:
	summary_result <- NULL
	for( i in seq(summary_list)){
	srow <- c("nsamples"=ncol(summary_list[[i]][[1]]), 
				"ntaxa"=nrow(summary_list[[i]][[1]]), 
				"nreads"=sum(summary_list[[i]][[1]]), 
				apply(summary_list[[i]][[2]], 2, function(x){length(unique(x))}))
	summary_result <- rbind(summary_result, srow)
	row.names(summary_result)[i] <- names(sapply(summary_list, "[", 1))[i]
	}
} else if (datatype == "keggs"){
	# Step-by-step datasets without keggs:
	nokegg_a <- select(a_mdata, -one_of(c("Total",keggvars)))
	nokegg_sg <- select(sg_mdata, -one_of(c("Total",keggvars)))
	# Step-by-step keggs without reads:
	noreads_a <- select(a_mdata, one_of(c("Total",keggvars)))
	noreads_sg <- select(sg_mdata, -one_of(c("Total",keggvars)))
	# Merge:
	summary_list <- list(list("Reads_brute" = nokegg_a, "Keggs_brute" = noreads_a), 
						list("Reads_filt1" = nokegg_sg, "Keggs_filt1" = noreads_sg))
	# Kegg summary:
	summary_result <- NULL
	for( i in seq(summary_list)){
	srow <- c("nsamples"=ncol(summary_list[[i]][[1]]), 
				"nkeggs"=nrow(summary_list[[i]][[1]]), 
				"nreads"=sum(summary_list[[i]][[1]]), 
				apply(summary_list[[i]][[2]], 2, function(x){length(unique(x))}))
	summary_result <- rbind(summary_result, srow)
	row.names(summary_result)[i] <- names(sapply(summary_list, "[", 1))[i]
	}
}

write.report(c("\n", "--- Filtering summary: Part 1:"))
tbl.report(summary_result)

## e) Collapse taxonomy (for next filter):
if (datatype == "taxonomy"){
	sg_mdata_collapse <- unite(sg_mdata, "Taxonomy", taxrng, sep=";")
	sg_mdata_collapse$Taxonomy <- gsub(pattern = " \\(no \\w.+ in NCBI\\)|k_", replacement = "", sg_mdata_collapse$Taxonomy)
	sg_mdata_collapse$Taxonomy <- gsub(pattern = ";\\w_", replacement = ".+", sg_mdata_collapse$Taxonomy)
	sg_mdata_collapse$Taxonomy <- gsub(pattern = "Unclassified", replacement = "unclassified", sg_mdata_collapse$Taxonomy)
	
	for(i in seq(nrow(sg_mdata_collapse))){
		comps <- unlist(strsplit(sg_mdata_collapse[i,1], split="\\.\\+"))
		newrow <- paste(unique(comps), collapse = '.+')
		sg_mdata_collapse[i,1] <- newrow
		rm(comps, newrow)
	}
}

# SAVE FILTERED TABLES -------------------------------------------------------------------------------------------------------------------------
write.report(c("\n", "--- Saving files..."))
write.table(sg_mdata, paste0(WD, "Filtered/", "F1.mreads.tsv"), quote = F, sep='\t', row.names=FALSE)
write.table(sg_phenot, paste0(WD, "Filtered/", "F1.metadata.csv"), sep=';', dec=',', quote = T, row.names=FALSE)
if(!isFALSE(genotypes)){
	write.table(sg_genot, paste0(WD, "Filtered/", "F1.genotypes.txt"), quote = F, sep=' ')
}
write.table(summary_result, paste0(WD, "Filtered/", "summary_stats.tsv"), quote=FALSE, sep='\t', col.names=NA)
if (datatype == "taxonomy"){
	write.table(sg_mdata_collapse, paste0(WD, "Filtered/", "F1.tax_reads_collapsed.tsv"), quote = F, sep='\t', row.names=FALSE)
}

write.report(c("\n", "--- SAMPLES FILTER ENDED ------------------------------------------------"))
