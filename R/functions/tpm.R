## tpm ##

library(biomaRt)
library(DESeq2)
library(GenomicFeatures)
library(reshape2)

txdb <- makeTxDbFromEnsembl(release = NA,
                            organism = "Danio rerio")
exonsByGene <- exonsBy(txdb1, by="gene")
TxLengths <- sapply(exonsByGene, 
                    function( gene ){ sum( width( reduce(unlist(gene)) ) ) }
)

dds <- estimateSizeFactors(dds)

tpm <- counts_to_tpm(counts(dds), TxLengths)

counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  return(tpm)
}

sampleCorHeatmap <- ggplot(data = sampleCor.m) + 
   geom_tile( aes( x = Var1, y = Var2, fill = value ),
                  colour = "grey60" ) + 
     scale_fill_gradientn( colours = c("blue", "yellow", "red"),
                           guide = guide_colorbar(title = "Correlation\nCoefficient\n(Pearson)") ) + 
     theme_void() + theme( legend.position="right",
                           legend.title = element_text(colour="black" ) )

plotList <- list()
add_to_plot_list <- function( plotList, plot, filename ){
  index <- length(plotList) + 1
  plotList[[index]] <- list( plot = plot, file = filename )
  return( plotList )
}
plotList <- add_to_plot_list( plotList, 
                             sampleCorHeatmap, file.path(rootPath, 'plots', 'Figure.1d.pdf') )
