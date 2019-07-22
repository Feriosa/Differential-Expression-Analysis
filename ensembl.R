source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

genes <- read.csv("onlyrosmap-coexpressed-CDK5R2.txt", sep = "\t", row.names=NULL)

genes <- genes[,1]

ensembl_genes<- genes

res= as.data.frame(getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "hgnc_symbol", "description"),
  values= ensembl_genes,
  mart= mart))

write.table(res, file = "onlyrosmap-coexpressed-CDK5R2-genenames.txt", sep = "\t", 
            row.names = F, col.names = T)

