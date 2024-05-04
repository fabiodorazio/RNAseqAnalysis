#for pieces of code
## for notes
### titles

### simple function to plot custom pca plot


## pca plot
.plot.pca <- function(x){
  coldata <- data.frame(stage = stage, condition = type, treatment = treatment)
  row.names(coldata) <- colnames(x)
  
  dds <- DESeqDataSetFromMatrix(countData = x,
                                colData = coldata,
                                design= ~ stage + condition)
  dds <- DESeq(dds)
  vsd <- vst(dds, blind=TRUE)
  
  pcaData <- plotPCA(vsd, intgroup=c("condition", "stage", "treatment"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  ggplot(pcaData, aes(PC1, PC2, color=group, shape=treatment)) +
    geom_point(size=3) + 
    theme_bw() + geom_point(aes(size = 4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_hline(yintercept=0, linetype="longdash", colour="grey", size=0.8) +
  geom_vline(xintercept=0, linetype="longdash", colour="grey", size=0.8) +
    scale_y_reverse(position = 'top') + scale_x_reverse()
  
}
