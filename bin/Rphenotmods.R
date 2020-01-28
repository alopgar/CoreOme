## R PHENOTYPIC DATA MODIFICATION
# You can add modifications for your phenotype data frame columns, or new columns.
sg_phenot$NLACTA <- as.factor(sg_phenot$NLACTA)
sg_phenot$DIASLE <- cut(sg_phenot$DIASLE, breaks = c(0, 70, 150, 1000), labels = c("<70", "70-150", ">150"), right = FALSE)
sg_phenot$CH4_Rank <- cut(sg_phenot$ch4.spkpm_adjusted, 
                          breaks = c(quantile(sg_phenot$ch4.spkpm_adjusted, prob = c(0, 0.25, 0.50, 0.75, 1))),
                          labels=c("LOW","L_MID","H_MID","HIGH"), include.lowest=TRUE)
