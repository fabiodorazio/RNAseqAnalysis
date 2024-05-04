library(ggplot2)

## starts with res object from DEseq2
#res <- results(dds)

#volcano for ATAC
custom.volcano.plot <- function(x, xlim=c(-10,10), ylim=40){
  x$significance <- "NS"
  #Tdrd7.frame$significance[(abs(Tdrd7.frame$Fold) > 3)] <- "FC"
  #x$significance[x$log2FoldChange > 1] <- "FC"
  #x$significance[x$log2FoldChange < -1] <- "FC1"
  
  x$significance[(x$padj < 0.05) & (x$log2FoldChange > 1)] <- "FDR"
  x$significance[(x$padj < 0.05) & (x$log2FoldChange < -1)] <- "FC_FDR"
  x <- data.frame(x)
  g <- ggplot(x, aes(x=x$log2FoldChange, y=-log10(x$padj))) +
    geom_point(aes(color=x$significance, alpha=x$significance, size=x$significance)) +
    scale_color_manual(values=c(NS = "darkgrey", FDR = "purple", FC_FDR ="darkcyan")) +
    scale_size_manual(values = c(2, 2, 1.5)) + scale_alpha_manual(values = c(1, 0.2, 1)) +
    xlim(xlim) + ylim(0,ylim)
  #g
  #return(head(x,10))
    
  g1= g + theme_bw(base_size=2) +theme(legend.background=element_rect(),
            plot.title=element_text(angle=0, size=12, face="bold", vjust=1),
            panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
            axis.text.x=element_text(angle=0, size=12, vjust=1),
            axis.text.y=element_text(angle=0, size=12, vjust=1),
            axis.title=element_text(size=12),legend.position="top",
            legend.key=element_blank(),legend.key.size=unit(0.5, "cm"),
            legend.text=element_text(size=8),title=element_text(size=8),
            legend.title=element_blank()) + xlab(bquote(~Log[2]~ "fold change")) +
    ylab(bquote(~-Log[10]~italic(Pvalue))) +
    geom_vline(xintercept=c(-1, 1), linetype="longdash", colour="black", size=0.4) +
    #yintercept is the -log10(pvalue): if pad < 0.05 -> -log10(0.05) = 1.3
    geom_hline(yintercept=1.3, linetype="longdash", colour="black", size=0.4) 
  g1+ geom_text(aes(label = ifelse((x$padj<0.05) & (abs(x$log2FoldChange)>20), as.character(x$start), "")), 
                hjust = 0, vjust = -0.25)
  
  #return(x) prints new res.all
  g1
}
#custom.volcano.plot(res)

