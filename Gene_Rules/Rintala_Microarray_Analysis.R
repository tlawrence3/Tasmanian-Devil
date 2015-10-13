source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("yeast2.db")
#biocLite("affyPLM")
library("GEOquery")
library("yeast2.db")
library("stringr")
library("affyPLM")

#Retrieve Rintala's data
Rintala <- getGEO("GSE12442")
Rintala <- Rintala[[1]]

#Check to see if it is already normalized (it is)
e_Rintala <- exprs(Rintala)
#eset_Rintala <- normalize.ExpressionSet.quantiles(Rintala)
#e_Rintala <- exprs(eset_Rintala)
boxplot(log(e_Rintala))

#Average the samples for the different conditions
e_Rintala_matrix <- as.matrix(e_Rintala)
Rintala_0 <- rowMeans(e_Rintala_matrix[,c("GSM312529", "GSM312530", "GSM312540", "GSM312541")], na.rm = TRUE)
Rintala_0.5 <- rowMeans(e_Rintala_matrix[,c("GSM312535", "GSM312536", "GSM312537", "GSM312538")], na.rm = TRUE) 
Rintala_1 <- rowMeans(e_Rintala_matrix[,c("GSM312539", "GSM312542", "GSM312543", "GSM312544", "GSM312545", "GSM312546")], na.rm = TRUE)
Rintala_2.8 <- rowMeans(e_Rintala_matrix[,c("GSM312527", "GSM312528", "GSM312531", "GSM312532")], na.rm = TRUE)
Rintala_20.9 <- rowMeans(e_Rintala_matrix[,c("GSM312533", "GSM312534", "GSM312547", "GSM312548")], na.rm = TRUE)

#This is the code to see the gene mappings from yeast2.db
#x <- yeast2ENSEMBL
#mapped_probes <- mappedkeys(x)
#xx <- as.list(x[mapped_probes])
#if(length(xx) > 0){
#  xx[1:5]
#  xx[[1]]
#}

#Map Rintala's probes from the Affymetrix Yeast Genome 2.0 Array
Rintala_genes=unlist(mget(colnames(Rintala_0), yeast2ENSEMBL))
#Remove probes that have NA's coming from S. pombe
Rintala_genes <- Rintala_genes[!is.na(Rintala_genes)]

#Convert the probes into genes
Rintala_probe_to_gene_0 <- c()
Rintala_probe_to_gene_0.5 <- c()
Rintala_probe_to_gene_1 <- c()
Rintala_probe_to_gene_2.8 <- c()
Rintala_probe_to_gene_20.9 <- c()
gene_names <- c()
for(i in names(Rintala_genes)){
  #If there are multiple genes mapped to 1 probe, append the same probe value for each different gene to a vector for the condition
  if (str_sub(i, start= -2)!="at"){
     Rintala_probe_to_gene_0 <- c(Rintala_probe_to_gene_0, Rintala_0[substr(i, 1, nchar(i)-1)])
     Rintala_probe_to_gene_0.5 <- c(Rintala_probe_to_gene_0.5, Rintala_0.5[substr(i, 1, nchar(i)-1)])
     Rintala_probe_to_gene_1 <- c(Rintala_probe_to_gene_1, Rintala_1[substr(i, 1, nchar(i)-1)])
     Rintala_probe_to_gene_2.8 <- c(Rintala_probe_to_gene_2.8, Rintala_2.8[substr(i, 1, nchar(i)-1)])
     Rintala_probe_to_gene_20.9 <- c(Rintala_probe_to_gene_20.9, Rintala_20.9[substr(i, 1, nchar(i)-1)])
     gene_names <- c(gene_names, Rintala_genes[i])
  }
  else{
    Rintala_probe_to_gene_0 <- c(Rintala_probe_to_gene_0, Rintala_0[i])
    Rintala_probe_to_gene_0.5 <- c(Rintala_probe_to_gene_0.5, Rintala_0.5[i])
    Rintala_probe_to_gene_1 <- c(Rintala_probe_to_gene_1, Rintala_1[i])
    Rintala_probe_to_gene_2.8 <- c(Rintala_probe_to_gene_2.8, Rintala_2.8[i])
    Rintala_probe_to_gene_20.9 <- c(Rintala_probe_to_gene_20.9, Rintala_20.9[i])    
    gene_names <- c(gene_names, Rintala_genes[i])
  }
}

#Rename the named vector entries with the gene names
names(Rintala_probe_to_gene_0) <- gene_names
names(Rintala_probe_to_gene_0.5) <- gene_names
names(Rintala_probe_to_gene_1) <- gene_names
names(Rintala_probe_to_gene_2.8) <- gene_names
names(Rintala_probe_to_gene_20.9) <- gene_names

#Average the genes that have multiple probes mapped to them
adjusted_Rintala_probe_to_gene_0 <- tapply(Rintala_probe_to_gene_0, names(Rintala_probe_to_gene_0), mean)
adjusted_Rintala_probe_to_gene_0.5 <- tapply(Rintala_probe_to_gene_0.5, names(Rintala_probe_to_gene_0.5), mean)
adjusted_Rintala_probe_to_gene_1 <- tapply(Rintala_probe_to_gene_1, names(Rintala_probe_to_gene_1), mean)
adjusted_Rintala_probe_to_gene_2.8 <- tapply(Rintala_probe_to_gene_2.8, names(Rintala_probe_to_gene_2.8), mean)
adjusted_Rintala_probe_to_gene_20.9 <- tapply(Rintala_probe_to_gene_20.9, names(Rintala_probe_to_gene_20.9), mean)

write.table(adjusted_Rintala_probe_to_gene_0, file="151006_Rintala_0.txt", sep = "\t")
write.table(adjusted_Rintala_probe_to_gene_0.5, file="151006_Rintala_0.5.txt", sep = "\t")
write.table(adjusted_Rintala_probe_to_gene_1, file="151006_Rintala_1.txt", sep = "\t")
write.table(adjusted_Rintala_probe_to_gene_2.8, file="151006_Rintala_2.8.txt", sep = "\t")
write.table(adjusted_Rintala_probe_to_gene_20.9, file="151006_Rintala_20.9.txt", sep = "\t")