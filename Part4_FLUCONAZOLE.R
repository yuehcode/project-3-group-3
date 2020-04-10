#R

library(DESeq2)
# sample information
info <- read.csv('group_3_rna_info.csv')
info <-info[c(4:6,10:12),]
# load counts
cnts <- read.csv('CAR_corn_oil_deseq.csv',row.names=1)

# filter out rows that have any zeros for funzies
cnts <- subset(cnts,rowSums(cnts==0)==0)

# sample information
#info <- read.csv('groups/group_EX_rna_info.csv')

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
res <- results(dds, contrast=c('mode_of_action','CAR/PXR','Control'))
res <- lfcShrink(dds, coef=2)

#Report the number of genes significant at p-adjust < 0.05
sum(res$padj < 0.05, na.rm=TRUE ) #3506


# write out DE results
write.csv(res,'FLUCONAZOLE_deseq_results.csv')

# A table of the top 10 DE genes 
res<- read.csv(header = T,'FLUCONAZOLE_deseq_results.csv')
gene_sortbypvalue <- (res[order(res$pvalue), ])
names(gene_sortbypvalue)[1]<-"Geneid"
write.csv(gene_sortbypvalue[1:10,] , 'top_ten_FLUCONAZOLE.csv')

# histograms of fold change values from the significant DE genes
resSig <- res[ which(res$padj < 0.05 ), ]
hist( resSig$log2FoldChange, breaks=100, col="grey",main = "Foldchange of Fluconacole group from the significant DE genes", xlab = "log2(foldchange)"  )


####make scatter plot
library(ggplot2)
library(ggrepel)
#read data to do the plot
res <- read.csv(header = T,'FLUCONAZOLE_deseq_results.csv' )

names(res)[1] <- "Geneid" 
# remove na value
res <- na.omit(res)
res$Significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, 
                          ifelse(res$log2FoldChange > 1, "Up", "Down"), "Stable")

volcano_Flu<- ggplot(
  # data&color
  res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = Significant), size=2) +
  scale_color_manual(values = c("dodgerblue","grey", "red")) +
  # gene_label
  geom_text_repel(
    data = subset(res, padj < 0.05 & abs(res$log2FoldChange) > 1),
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

volcano_Flu+ggtitle("Volcano plot for Fluconacole group")



