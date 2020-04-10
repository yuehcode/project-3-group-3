

library(DESeq2)


# load counts
cnts <- read.csv('AhR_corn_oil_deseq.csv',row.names=1)

# load information sheet
info <- read.csv('group_3_rna_info.csv')
info <-info[c(1:3,10:12),]

# filter out rows that have any zeros for funzies
cnts <- subset(cnts,rowSums(cnts==0)==0)

# sample information
#info <- read.csv('group_3_AhR_rna_info.csv')

# create the DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = cnts,
  colData = info,
  design= ~ mode_of_action
)

# relevel mode_of_action as factor
dds$mode_of_action <- relevel(dds$mode_of_action, ref='Control')

# run DESeq
dds <- DESeq(dds)
res <- results(dds, contrast=c('mode_of_action','AhR','Control'))
res <- lfcShrink(dds, coef=2)

# write out DE results
write.csv(res,'LEFLUNOMIDE_AhR_deseq_results.csv')



#Report the number of genes significant at p-adjust < 0.05
sum(res$padj < 0.05, na.rm=TRUE ) #1379

# A table of the top 10 DE genes 
res<- read.csv(header = T,'LEFLUNOMIDE_AhR_deseq_results.csv')
gene_sortbypvalue <- (res[order(res$pvalue), ])
names(gene_sortbypvalue)[1]<-"Geneid"
write.csv(gene_sortbypvalue[1:10,] , 'top_ten_LEFLUNOMIDE.csv')

# histograms of fold change values from the significant DE genes
resSig <- res[ which(res$padj < 0.05 ), ]
hist( resSig$log2FoldChange, breaks=100, col="grey", main = "Foldchange of Leflunomide group from the significant DE genes", xlab = "log2(foldchange)" )


### make scatter plot
library(ggplot2)
library(ggrepel)

#read data to do the plot
res <- read.csv(header = T,'LEFLUNOMIDE_AhR_deseq_results.csv' )
names(res)[1] <- "Geneid"

# remove na value
res <- na.omit(res)
# make new column for significant values.
res$Significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, 
                          ifelse(res$log2FoldChange > 1, "Up", "Down"), "Stable")

volcano_LEF<- ggplot(
  # data&color
  res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = c("dodgerblue","grey", "red")) +
  # gene_label
  geom_text_repel(
    data = subset(res,  pvalue < 0.05 & abs(res$log2FoldChange) > 1),
    aes(label = Geneid),
    size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")) +
  # dot line
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  # label axis 
  labs(x="Effect size: log2(fold-change)",
       y="-log10 (p-value)") +
  # graph position
  theme(legend.position = "bottom")

volcano_LEF+ggtitle("Volcano plot for Leflunomide group")

