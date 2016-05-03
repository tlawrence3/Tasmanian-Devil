source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("yeast2.db")
#biocLite("affyPLM")
library("GEOquery")
library("yeast2.db")
library("stringr")
library("affyPLM")

#Retrieve Rintala's data
Gasch <- getGEO("GDS21")
eset <- GDS2eSet(Gasch)
gene_names <- eset@featureData@data$Platform_ORF

#Check to see if it is already normalized (it is not)
eset_Gasch <- normalize.ExpressionSet.quantiles(eset)
e_Gasch <- exprs(eset_Gasch)
boxplot(log(e_Gasch))


#Retrieve the sample for glucose and ethanol
e_Gasch_matrix <- as.matrix(e_Gasch)
Gasch_glucose <- e_Gasch_matrix[,"GSM990"]
Gasch_ethanol <- e_Gasch_matrix[,"GSM1001"]

#Rename the named vector entries with the gene names and only keep genes with values
names(Gasch_glucose) <- gene_names
names(Gasch_ethanol) <- gene_names
Gasch_glucose_reduced <- Gasch_glucose[!is.na(Gasch_glucose)]
Gasch_ethanol_reduced <- Gasch_ethanol[!is.na(Gasch_ethanol)]
Gasch_glucose_reduced_orfs <- Gasch_glucose_reduced[! names(Gasch_glucose_reduced) %in% c("")]
Gasch_ethanol_reduced_orfs <- Gasch_ethanol_reduced[! names(Gasch_ethanol_reduced) %in% c("")]

#Average the genes that have multiple probes mapped to them
Gasch_glucose_reduced_orfs_average <- tapply(Gasch_glucose_reduced_orfs, names(Gasch_glucose_reduced_orfs), mean)
Gasch_ethanol_reduced_orfs_average <- tapply(Gasch_ethanol_reduced_orfs, names(Gasch_ethanol_reduced_orfs), mean)

write.table(Gasch_glucose_reduced_orfs_average, file="151012_Gasch_glucose.txt", sep = "\t")
write.table(Gasch_ethanol_reduced_orfs_average, file="151012_Gasch_ethanol.txt", sep = "\t")