library(biomaRt)
library(DESeq2)
library(GenomicFeatures)

setwd("~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/perGeneCounts_Stages/")
files <- list.files(pattern = '\\.tab')

# creating a function which reads a csv file and make a matrix
readTagsPerGene <- function(x){
  out <- read.csv(x, skip = 4, header = FALSE, sep = "\t", row.names = 1,
                  col.names = c("","totalCount", "forwardCount", "reverseCount"))
  out
}


allReadsCounts <- lapply(files, readTagsPerGene)
head(allReadsCounts[[1]])

totalReadsCount <- sapply(allReadsCounts, function(x) x$totalCount) # makes a matrix with the values from the totalCount column
rownames(totalReadsCount) <- row.names(allReadsCounts[[1]])
colnames(totalReadsCount) <- c('s256PGC1', 's256Soma1', 's256PGC2', 's256Soma2', 'sHighPGC1', 'sHighSoma1', 'sHighPGC2', 'sHighSoma2',  'sDomePGC1', 'sDomeSoma1', 'sDomePGC2', 'sDomeSoma2', 's10somitesPGC1', 's10somitesSoma1', 's10somitesPGC2', 's10somitesSoma2', 'sPrim5PGC1', 'sPrim5Soma1', 'sPrim5PGC2', 'sPrim5Soma2')
head(totalReadsCount)

## calculate gene length
human <- useMart("ensembl", dataset="drerio_gene_ensembl")
start_pos = getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id" , mart = human, attributesL = "start_position", martL = human, uniqueRows=T)
end_pos = getLDS(attributes = "ensembl_gene_id", filters = "ensembl_gene_id" , mart = human, attributesL = "end_position", martL = human, uniqueRows=T)
gene_L <- merge(start_pos, end_pos, by.x="Gene.stable.ID", by.y="Gene.stable.ID")
gene_L$Length <- gene_L$Gene.end..bp. - gene_L$Gene.start..bp.

gene_lenght <- data.frame(row.names = gene_L$Gene.stable.ID, Length = gene_L$Length)
Lengths_counts <- merge(totalReadsCount, gene_lenght, by = 0)


 #####rpkm#####
#calulate rpm by dividing each read by the 'per million scaling factor'
transcript_level_expr <- sapply(Lenghts_counts, function(x) x/sum(x/1e6))
transcript_level_expr_frame <- data.frame(transcript_level_expr, Lenghts_counts$TxLengths, row.names = rownames(Lenghts_counts))
rpkm <- sapply(transcript_level_expr_frame, function(y) y/(transcript_level_expr_frame$Lenghts_counts.TxLengths/1e3))
rpkm <- data.frame(rpkm, row.names = rownames(Lenghts_counts))

##merge replicates
rpkm$MeanPGC <- rowMeans(rpkm[c('s256PGC1', 's256PGC2')])
rpkm$MeanSOma <- rowMeans(rpkm[c('s256Soma1', 's256Soma2')])
rpkm$MeanPGChigh <- rowMeans(rpkm[c('sHighPGC1', 'sHighPGC2')])
rpkm$MeanSomahigh <- rowMeans(rpkm[c('sHighSoma1', 'sHighSoma2')])
rpkm$MeanPGCprim <- rowMeans(rpkm[c('sPrim5PGC1', 'sPrim5PGC2')])
rpkm$MeanSomaprim <- rowMeans(rpkm[c('sPrim5Soma1', 'sPrim5Soma2')])
#filter rpkm > 1
x <- data.frame(rpkm[,c(23:26)])
y <- subset(x, c(x$MeanPGC > 1 & x$MeanSOma >1, x$MeanSomahigh > 1, x$MeanPGChigh > 1))
y <- na.omit(y)
#transform rpkm in concentrations

PGC.genes <- read.table('Upregulated_MBT_PGC_ony.xls', sep = '\t')
pgc.genes.ID <- rownames(PGC.genes)
rpkm.pgc.genes <- x[pgc.genes.ID,]

 #####tpm#####
#divide each column by the gene length in kb
#transcript_level_expr <- sapply(Lenghts_counts, function(x) x/(Lenghts_counts$TxLenghts/1000))
#rownames(transcript_level_expr) <- rownames(Lenghts_counts)
#transcript_level_expr1 <- data.frame(transcript_level_expr)

transcript_level_expr <- Lenghts_counts/(Lenghts_counts$TxLengths/1e3)
tpm1 <- apply(transcript_level_expr, 2,  function(x) {x/(sum(x, na.rm = T)/1e6)})
tpm1 <- as.data.frame(tpm1)

#average columns
tpm_white$s128mean <- rowMeans(tpm_white[c('s128', 's128.1', 's128.2', 's128.3', 's128.4')])

tpmOrder <- data.frame(tpm1[, c(1,3,2,4,5,7,6,8,9,11,10,12,13,15,14,16,17,19,18,20)])

pheatmap(log(tpmOrder + 1), cluster_rows= TRUE,
         cluster_cols=FALSE, fontsize_col = 5.8, show_rownames = F)

# filter for tpm = 0 in all the stages
tpmByGeneList <- split.data.frame(tpmOrder, 
                                  factor(rownames(tpmOrder), 
                                         levels = rownames(tpmOrder)))
# divide data up by Gene and then Stage

rnaseqSampleInfo <- data.frame(stageName = colnames(tpmOrder))
tpmByGeneByStageList <- lapply(tpmByGeneList,
                               function(x){ 
                                 split.data.frame(t(x), rnaseqSampleInfo$stageName) } )

# filter to genes with non-zero TPM in all replicates at at least one stage
tpmThreshold <- 0
genesAboveThreshold <- sapply( tpmByGeneByStageList,
                               function(geneByStageList){
                                 any( sapply(geneByStageList, 
                                             function(x){ all(x > tpmThreshold) } ) )
                               }
)
tpmFilt <- tpmOrder[genesAboveThreshold, ]

# plot sample correlation
library(ggplot2)
library(reshape2)
sampleCor <- cor(tpmFilt)
sampleCor.m <- melt(sampleCor)

sampleCorHeatmap <- ggplot(data = sampleCor.m) + 
  geom_tile( aes( x = Var1, y = Var2, fill = value ),
             colour = "grey60" ) + 
  scale_fill_gradientn( colours = c("blue", "yellow", "red"),
                        guide = guide_colorbar(title = "Correlation\nCoefficient\n(Pearson)") ) + 
  theme_void() + theme( legend.position="right",
                        legend.title = element_text(colour="black" ) )
plotList <- add_to_plot_list(plotList, sampleCorHeatmap, file.path(rootPath, 'plots', 'Figure.1d.pdf'))

print(sampleCorHeatmap)
