## PLOT PCA components

plotPCA.comp <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object)) ## calculates row variance
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))] ## return row numbers from higher variance to lower
  pca <- prcomp(t(assay(object)[select, ])) ## principal component analysis
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC3", color = "group")) + geom_point(size = 3)+
    xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
    ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + 
    coord_fixed() 
}

