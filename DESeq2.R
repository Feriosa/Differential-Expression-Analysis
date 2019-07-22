source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)
library(EnhancedVolcano)

counts <- read.delim("slr1340_count.txt", header=TRUE, row.names=1)

counts <- as.matrix(counts)

samples <- data.frame(condition = rep(c("NT0", "slr1340"), each=3))

(condition <- factor(c(rep("NT0", 3), rep("slr1340", 3))))

ds <- DESeqDataSetFromMatrix(countData=round(counts), colData=samples, design=~condition)

ds <- DESeq(ds)

res <- results(ds)
table(res$padj<0.05)
res <- res[order(res$padj),]

write.csv(as.data.frame(res), 
          file="slr1340_results.csv")


res<-read.csv('slr1340_results.csv',sep=',',header=TRUE)

hist(res$pvalue, breaks=20, col="grey",main="DESeq2 slr1340 p-value distribution", xlab="Nominal P-value", ylab="Number of genes")

png("slr1340-volcanoplot2.png", 1200, 1000, pointsize=20)
EnhancedVolcano(res,
                
                lab = res$X,
                
                x = "log2FoldChange",
                
                y = "pvalue",
                transcriptPointSize = 3.0,
                transcriptLabSize= 9.0,
                title="slr1340",
                titleLabSize = 32)
dev.off()


