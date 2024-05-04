
library(DESeq2)
setwd('~/Desktop/Postdoc/Data_Analysis/RNAseq/Dnd/')
files_2 <- list.files(pattern = '\\.tab')
## creating a function which reads a csv file and make a matrix
readTagsPerGene <- function(x){
  out <- read.csv(x, skip = 4, header = FALSE, sep = "\t", row.names = 1,
                  col.names = c("","totalCount", "forwardCount", "reverseCount"))
  out
}

allReadsCounts <- lapply(files_2, readTagsPerGene)
# keeps only the total count column
totalReadsCount_2 <- sapply(allReadsCounts, function(x) x$totalCount)
rownames(totalReadsCount_2) <- row.names(allReadsCounts[[1]])

#spl_col <- strsplit(files, "\\Read")
#x <- lapply(spl_col, function(x) x[1])

## get column names
spl_names <- strsplit(files_2, "\\_Z")
col_names <- as.character(lapply(spl_names, function(u) u[[1]][1]))
colnames(totalReadsCount_2) <- col_names

totalReadsCount_2 <- totalReadsCount_2[,-c(6,8,12,14)]
colnames(totalReadsCount_2)

condition <- c(rep('PGC', times = 6), rep(rep(c('Soma', 'PGC'), each = 2), times = 2), 'PGC', rep(c('PGC', 'Soma'), each = 2))
treatment <- c(rep('Dndmm', times = 2), rep('DndMO', times = 2), rep('Tdrd7mm', times = 4), 
               rep('Tdrd7MO', times = 4), rep('Rescue', times = 3), rep('wt', times = 4))


run.deseq.genes <- function(x, condition= condition, treatment= NULL){
  if(is.null(treatment)){
    coldata <- data.frame(condition = condition)
    row.names(coldata) <- colnames(x)
    dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ condition)
  }
  else{
    coldata <- data.frame(treatment = treatment,condition = condition)
    row.names(coldata) <- colnames(x)
    dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ condition + treatment)
  }
  dds <- DESeq(dds)
  
  return(dds)
}
good_rep <- run.deseq.genes(totalReadsCount_2, condition, treatment)
##correlation heatmap
pheatmap(cor(totalReadsCount_2)[-c(1:2,13:15,16:19),-c(1:2,13:15,16:19)], fontsize = 5)

### diff expressed genes compared to mm
dnd_deg <- totalReadsCount_2[,c(3:6)]
condition_dnd <- c('MO', 'MO', 'mm', 'mm')
dnd_dds <- run.deseq.genes(dnd_deg, condition_dnd)
dnd_res <- results(dnd_dds)
## subset deg 
up_dnd <- subset(dnd_res, dnd_res$log2FoldChange > 0 & dnd_res$padj < 0.1)
down_dnd <- subset(dnd_res, dnd_res$log2FoldChange < 0 & dnd_res$padj < 0.1)

## compare to Tdrd7 deg
tdrd7_deg <- totalReadsCount_2[,c(9,10, 5,6)]
condition_tdrd7 <- c('MO', 'MO', 'mm', 'mm')
tdrd7_dds <- run.deseq.genes(tdrd7_deg, condition_tdrd7)
tdrd7_res <- results(tdrd7_dds)
## subset deg 
up_tdrd7 <- subset(tdrd7_res, tdrd7_res$log2FoldChange > 0 & tdrd7_res$padj < 0.1)
down_tdrd7 <- subset(tdrd7_res, tdrd7_res$log2FoldChange < 0 & tdrd7_res$padj < 0.1)

## overlap genes
x <- rownames(down_tdrd7)
dim(na.omit(as.data.frame(down_dnd)[x,]))
write.table(down_dnd, "Tables/DownregDndMOvsTdrd7mm.txt", sep = '\t')
write.table(down_tdrd7, "Tables/DownregTdrd7MOvsTdrd7mm.txt", sep = '\t')

y <- rownames(up_tdrd7)
dim(na.omit(as.data.frame(up_dnd)[y,]))
write.table(up_dnd, "Tables/UpregDndMOvsTdrd7mm.txt", sep = '\t')
write.table(up_tdrd7, "Tables/UpregTdrd7MOvsTdrd7mm.txt", sep = '\t')

## GO:0048729,GO:0060322
dev_genes <- read.csv("~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/Dev_genesGO0048729_GO0060322.txt")
dev_genes <- unique(dev_genes$Gene.stable.ID)

