if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

counts <- read.delim("h3k4me3.TPM_normalized_coverages.tsv", header=TRUE, row.names=1)
counts <- as.matrix(counts)
coldata <- data.frame(samples=c('CR025_wo_S16.sam.filtered.sorted.bam','CR028_wo_S3.sam.filtered.sorted.bam','CR031_S6.sam.filtered.sorted.bam', 'CR031_S6.sam.filtered.sorted.bam'), spec = c('human', 'human','chimp','chimp'))
(spec=factor(c("human","human", "chimp", "chimp")))
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=coldata, design =~spec)
dds <- DESeq(dds)
res2 <- results(dds)

res2 <- res2[order(res2$log2FoldChange),]
write.csv(as.data.frame(res2), 
          file="denemee.csv")

res2<-read.csv('denemee.csv',sep=',',header=TRUE)

png("h3k4me3-volcanoplot3.png", 1200, 1000, pointsize=20)
EnhancedVolcano(res2,
                
                lab = res2$X,
                
                x = "log2FoldChange",
                y = "pvalue",
                transcriptPointSize = 3.0,
                transcriptLabSize= 9.0,
                title="h3k4me3",
                titleLabSize = 32)
dev.off()

