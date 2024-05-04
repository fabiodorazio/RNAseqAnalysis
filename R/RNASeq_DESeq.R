#for pieces of code
## for notes
### titles

### run DESeq on gene counts

library('DESeq2')
library('ggplot2')

setwd("~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/perGeneCounts_Stages/")
files <- list.files(pattern = '\\.tab')

## creating a function which reads a csv file and make a matrix
readTagsPerGene <- function(x){
  out <- read.csv(x, skip = 4, header = FALSE, sep = "\t", row.names = 1,
                  col.names = c("","totalCount", "forwardCount", "reverseCount"))
  out
}


allReadsCounts <- lapply(files, readTagsPerGene)
head(allReadsCounts[[1]])

totalReadsCount <- sapply(allReadsCounts, function(x) x$totalCount) # makes a matrix with the values from the totalCount column
rownames(totalReadsCount) <- row.names(allReadsCounts[[1]])

spl_col <- strsplit(files, "\\Read")
spl_col <- as.character(lapply(spl_col, function(x) x[1]))
colnames(totalReadsCount) <- spl_col
## create lists for coldata columns
stage <- c(rep(c('s256', 'high', 'dome', 'somites10', 'prim5'), each = 4),rep('prim5', times =8))
type <- c(rep(c('PGC', 'Soma'), times = 14))
treatment <- c(rep('wt', times = 20), rep('MO', times = 4), rep('5mm', times = 4))

## deseq requires a reference dataframe named coldata
run.deseq.genes <- function(x){
  coldata <- data.frame(stage = stage, condition = type, treatment = treatment)
  row.names(coldata) <- colnames(x)

  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ stage + condition)
  head(dds)
  dds <- DESeq(dds)
  #normalizedData <- counts(dds, normalized = TRUE) # divides the dds by the size factor. holds the counts in a matrix. Normalized = T devides by the norm factor
  #head(normalizedData)
  #dds$condition <- factor(dds$condition, levels = c("s128","s1k")) # tells deseq what is the control. Most of times not necessery
  
  res <- results(dds) 
  ## for testing logFC above or belove a threshold:
  #res <- results(dds, lfcThreshold = .5, altHypothesis = 'greaterAbs')
  
  #plotMA(res)
  #drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
  #plotMA(res, ylim = c(-10,10)); drawLines()
  
  return(dds)
}

#run the function to obtain differential genes and MA plot
dds <- run.deseq.genes(totalReadsCount)
res <- results(dds)
## plotMA logFC vs normalised reads
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(DEgenes, ylim = c(-10,10)); drawLines()

mcols(res, use.names=TRUE) # check conditions of the experiment

## normalisation for gene-gene comparison
mat <- counts(dds, normalized=TRUE) #normalised read counts (raw/size factor)

ntd <- normTransform(dds) # log2 +1
head(assay(ntd))


#calculate variance for PCA plot
vsd <- vst(dds, blind=TRUE)
plotPCA(vsd, intgroup=c("condition", "stage", "treatment"))

library(ggplot2)
pcaData <- plotPCA(vsd, intgroup=c("condition", "stage", "treatment"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=group, shape=treatment)) +
  geom_point(size=3) + 
  scale_color_manual(values = c("green3", "lightgreen", "deepskyblue4",
                                "cyan", "blue","orchid3", "lightpink2", "mediumpurple4",
                                "mediumorchid4", "violet",  "mediumpurple", "forestgreen",
                                "mediumpurple1", "steelblue3", "yellow")) +
  theme_bw() + geom_point(aes(size = 4)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept=0, linetype="longdash", colour="grey", size=0.8) +
  geom_vline(xintercept=0, linetype="longdash", colour="grey", size=0.8) +
  scale_y_reverse(position = 'top') + scale_x_reverse()

## Rescue DEG
rescue_mm <- totalReadsCount[,c(25,27,29,30)]
rescue_MO <- totalReadsCount[,c(21,23,29,30)]
condition_rescue <- c('MO', 'MO', 'Rescue', 'Rescue')
run.deseq.genes2 <- function(x, condition= condition){
  coldata <- data.frame(condition = condition)
  row.names(coldata) <- colnames(x)
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  res.sub <- subset(res, res$padj < 0.1)
  return(res.sub)
}
res_mm <- run.deseq.genes2(rescue_mm, condition_rescue)
res_MO <- run.deseq.genes2(rescue_MO, condition_rescue)

# extracting values

resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.1)
resUp <- subset(resSig, resSig$log2FoldChange < 0)
resDown <- subset(resSig, resSig$log2FoldChange > 0)

write.table(resOrdered, '24hpf_PGC&Soma_resOrdered.txt', sep='\t') # saveRDS(resOrdered, '24hpf_PGC&Soma_resOrdered.rds')
write.table(resSig, '24hpf_PGC&Soma_resSig.txt', sep='\t') # saveRDS(resSig, '24hpf_PGC&Soma_resSig.rds')

### Heat map
library("pheatmap")
library(RColorBrewer)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:30]
df <- as.data.frame(colData(dds))
rownames(df) <- colnames(ntd)

breaksList = seq(0, 10, by = 1)

pheatmap(assay(ntd)[polycomb2,], color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList, cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=df)

## select high pvalues genes for heatmap
s <- row.names(resSig)
s <- head(s, 30)
assay(dds[s])
s <- log2(assay(dds[s]))
pheatmap(s, cluster_rows=FALSE,
         cluster_cols=FALSE, annotation_col=df)

# heatmap of k-mean clustered tpms

km256 <- kmeans(kmean256, 5)
km256 <- cbind(kmean256, km256$cluster)
km256 <- km256[order(km256$`km256$cluster`),]
km256$`km256$cluster` <- NULL
pheatmap(log(km256 + 1), cluster_cols = FALSE, cluster_rows = FALSE)

library(biomaRt)
library(pheatmap)
## heatmap for TFs
GermDevGenes <- assay(ntd)[c("ENSDARG00000053569", "ENSDARG00000051926", "ENSDARG00000031777", "ENSDARG00000005915",
       "ENSDARG00000071497", "ENSDARG00000070913", "ENSDARG00000104307", "ENSDARG00000011235",
       "ENSDARG00000077467", "ENSDARG00000015536", "ENSDARG00000069342", "ENSDARG00000003293",
       "ENSDARG00000069726"),]

GermDevGenes <- assay(ntd)[c('ENSDARG00000014373', 'ENSDARG00000022813',
                             'ENSDARG00000013453', 'ENSDARG00000071450', 'ENSDARG00000040732',
                             'ENSDARG00000077523', 'ENSDARG00000036214',
                             'ENSDARG00000008454', 'ENSDARG00000075113',
                             'ENSDARG00000044774', 'ENSDARG00000070913', 
                             'ENSDARG00000079922',
                             'ENSDARG00000091001', 'ENSDARG00000016999','ENSDARG00000103379',
                             'ENSDARG00000068567', 'ENSDARG00000043457', 'ENSDARG00000055158',
                             'ENSDARG00000028071', 'ENSDARG00000014321', 'ENSDARG00000038990',
                             'ENSDARG00000055554', 'ENSDARG00000014420', 'ENSDARG00000021184', 'ENSDARG00000079305',
                             'ENSDARG00000058917', 'ENSDARG00000037507',
                             'ENSDARG00000102731', 'ENSDARG00000041959', 'ENSDARG00000053569',
                             'ENSDARG00000015906', 'ENSDARG00000043802', 'ENSDARG00000078754',
                             'ENSDARG00000002601',
                             'ENSDARG00000010571', 'ENSDARG00000019995',
                             'ENSDARG00000045139', 'ENSDARG00000039310', 'ENSDARG00000027957', 
                             'ENSDARG00000018404', 'ENSDARG00000041952'),]
housekeeping <- assay(ntd)[c('ENSDARG00000007682', 'ENSDARG00000101061', 'ENSDARG00000004754', 
                             'ENSDARG00000044267', 'ENSDARG00000032175', 'ENSDARG00000044301', 'ENSDARG00000038835'),]
#assign gene names from biomart
ensembl = useEnsembl(biomart="ensembl", dataset='drerio_gene_ensembl', host = "www.ensembl.org")
mart <- useMart("drerio_gene_ensembl", useMart("ensembl"),verbose = FALSE)

s4 <- rownames(GermDevGenes)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=s4,mart= mart)
row.names(G_list) <- G_list$ensembl_gene_id
G_list$ensembl_gene_id <- NULL
GermDevGenes <- merge(GermDevGenes,G_list, by = 'row.names')
row.names(GermDevGenes) <- GermDevGenes$external_gene_name
GermDevGenes$Row.names <- NULL
GermDevGenes$external_gene_name <- NULL
GermDevGenesOrder <- GermDevGenes[order(GermDevGenes$PGC1),]
pheatmap(GermDevGenesOrder, cluster_rows = FALSE, cluster_cols = FALSE)


#biomart
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset='drerio_gene_ensembl')
mart <- useDataset("drerio_gene_ensembl", useMart("ensembl"), verbose = FALSE)
#obtain gene lengths
attr <- listAttributes(mart)
attr[grep('length', attr$name),]
#obtain gene names
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=s4,mart= mart)


library(EDASeq)
library(plyr)
s256PGC <- read.csv('q256Cells_PGCsBucRep1_filt_Zv10ChrERCC92M1000ReadsPerGene.out.tab', skip = 4, header = F, sep = '\t')
gene_list <- data.frame(s256PGC$V1)
gene_lengths <- getGeneLengthAndGCContent(gene_list$s256PGC.V1, 'dre') # try a couple of time if doesn't work
join(gen, s256PGC, by = 'V1', type = 'left') -> IDs_lentghs
IDs_lentghs$V3 <- NULL
IDs_lentghs$V4 <- NULL
IDs_lentghs$gc <- NULL
IDs <- na.omit(IDs_lentghs)


#plot go terms
library(org.Dr.eg.db)
library(clusterProfiler)
ego <- enrichGO(rownames(resUp), OrgDb = org.Dr.eg.db, keyType = 'ENSEMBL', ont = 'BP', pvalueCutoff = 0.1, universe = rownames(totalReadsCount))
dotplot(ego, showCategory = 25, font.size = 8)
ego[2]$geneID #visualize gene ids for one category

# Get the stacked barplot
frame.count <- data.frame(c(311,1618), row.names = c('zygotic', 'maternal and zygotic'))
frame.count <- data.matrix(frame.count)
barplot(frame.count, col=c('blue', 'orange') , border=NULL, space=0.04, font.axis=2, legend=c('zygotic', 'maternal and zygotic'))

## plot correlation of fold change
## the fold change of B in respect to A is B/A

fold.change.tpm <- data.frame(foldchangePGC = (heyn.tpm$sHighPGC1+1)/(heyn.tpm$s256PGC1+1),
                              foldchangeSoma = (heyn.tpm$sHighSoma1+1)/(heyn.tpm$s256Soma1+1))
plot(log(fold.change.tpm))


## correlation plots
normalizedData.oocyte.pgc <- normalizedData[,c(1,17)]
normalizedData.oocyte.latePGC <- data.frame(EarlyPGC = rowMeans(normalizedData[c(1,3)],))
normalizedData.oocyte.latePGC$LatePGC <- rowMeans(normalizedData[c(17,19)],)
## subset genes
normalizedData.oocyte.latePGC.sub <- normalizedData.oocyte.latePGC[rownames(pgc.threshold.bigger.than.5),]
x <- merge(pgc.threshold.bigger.than.5, cage.PGC.merged, by.x = 0, by.y = 'ENSEMBL')

plot(normalizedData.oocyte.latePGC, cex = 0.4, 
     pch = 16, col = rgb(red = 0, green = 0.2, blue = 0.8, alpha = 0.3))
points(normalizedData.oocyte.latePGC[unique(x$Row.names),], col = 'red', cex = 0.6, pch = 16)


## zygotic vs maternal transcripts
tpm <- read.csv('~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables/tpmStagesFiltered.txt', sep = '\t')
zygotic_transcripts <- read.csv('~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables/ZygoticTranscripts_tpm_threshold2.txt', sep = '\t')
up_high_dome <- read.csv('~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables/Diff_Expressed_Genes/UPREGinPGCfromHighToDome.txt', sep = '\t')
up_256_high <- read.csv('~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables/MBT/UpregGenesInPGConly256toHighpadj01.txt', sep = '\t')
# make a list and subset the genes 
# arrange genes based on zygotic or maternal
my_transcript_list <- list(up_256_high, up_high_dome)
names(my_transcript_list) <- c('Up_256_high', 'Up_high_dome')

subset_transcripts <- function(x){
  up_tpm <- na.omit(x[rownames(tpm),])
  up_tpm$Expression <- ifelse(rownames(up_tpm) %in% rownames(zygotic_transcripts),'Maternal', 'Zygotic')
  return(up_tpm)
}

my_transcript_list <- lapply(my_transcript_list, subset_transcripts)
# save
write.table(my_transcript_list[[1]], '~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables/Upregulated_Transcript_List_256_high_PGConly_Maternal_vs_zygotic. txt', sep = '\t')
write.table(my_transcript_list[[2]], '~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables/Upregulated_Transcript_List_high_dome_PGConly_Maternal_vs_zygotic. txt', sep = '\t')
