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
colnames(totalReadsCount) <- c('s256PGC1', 's256Soma1', 's256PGC2', 's256Soma2', 'sHighPGC1', 'sHighSoma1', 'sHighPGC2', 'sHighSoma2',  'sDomePGC1', 'sDomeSoma1', 'sDomePGC2', 'sDomeSoma2', 's10somitesPGC1', 's10somitesSoma1', 's10somitesPGC2', 's10somitesSoma2', 'sPrim5PGC1', 'sPrim5Soma1', 'sPrim5PGC2', 'sPrim5Soma2', 'Tdrd7pgcMO1', 'Tdrd7somaMO1', 'Tdrd7pgcMO2', 'Tdrd7somaMO2', 'Tdrd7pgc5mm1', 'Tdrd7soma5mm1', 'Tdrd7pgc5mm2', 'Tdrd7soma5mm2')
head(totalReadsCount)
## create lists for coldata columns
condition <- rep(c('PGC', 'Soma'), times = 2)


## load tpm
tpm <- read.csv("~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables/tpmStagesFiltered.txt", sep = '\t')

## deseq requires a reference dataframe named coldata
run.deseq.genes <- function(input1, input2 = FALSE, condition = condition){
  if(input2==FALSE){
    x.input3 <- grep(input1, colnames(totalReadsCount))
  }

  else{
    x.input2 <- grep(input2, colnames(totalReadsCount))
    x.input1 <- grep(input1, colnames(totalReadsCount))
    x.input3 <- c(x.input1, x.input2)
  }
  total.r <- totalReadsCount[,x.input3]
    
  ## make coldata
  coldata <- data.frame(condition = condition)
  row.names(coldata) <- colnames(total.r)

  dds <- DESeqDataSetFromMatrix(countData = total.r, colData = coldata, design = ~ condition)
  dds <- DESeq(dds)

  res <- results(dds) 
  resSig <- subset(res, padj < 0.1)

  return(resSig)
}

#run the function to obtain differential genes and MA plot
condition <- rep(c('PGC', 'Soma'), times = 2)
tpm.256 <- run.deseq.genes('256', condition = condition)
tpm.high <- run.deseq.genes('High', condition = condition)
tpm.dome <- run.deseq.genes('Dome', condition = condition)
tpm.10somites <- run.deseq.genes('10somites', condition = condition)
tpm.prim5 <- run.deseq.genes('Prim5', condition = condition)

## PGC
condition <- rep(c('early', 'late'), each = 2)
tpm.256_PGC <- run.deseq.genes('256PGC', 'HighPGC', condition = condition)
tpm.HIGH_PGC <- run.deseq.genes('HighPGC', 'DomePGC', condition = condition)

## get tpm
get_tpm <- function(input1, input2, input3, sample_names = list()){
  x.input1 <- grep(sample_names[1], colnames(totalReadsCount))
  y.input2 <- grep(sample_names[2], colnames(totalReadsCount))
  if(length(sample_names) == 3){
    z.input3 <- grep(sample_names[3], colnames(totalReadsCount))
    list_of_samples <- c(x.input1, y.input2, z.input3)
    names.ensembl <- c(rownames(input1), rownames(input2), rownames(input3))
    
  }
  else{
    list_of_samples <- c(x.input1, y.input2)
    names.ensembl <- c(rownames(input1), rownames(input2))
    
  }

  tpm <- tpm[,list_of_samples]
  tpm.sub <- tpm[names.ensembl,]
  tpm.sub <- na.omit(tpm.sub)
}

heatmap_late <- get_tpm(tpm.dome, tpm.10somites, tpm.prim5, sample_names = c('Dome','10somites', 'Prim5'))
heatmap_early <- get_tpm(tpm.256,tpm.high, tpm.dome, sample_names = c('256','High', 'Dome'))

## for PGC only the function above does not work as the columns are reshuffled
tpm.pgc <- tpm[,c(1,2,5,6,9,10)]
list_names_pgc <- c(rownames(tpm.256_PGC), rownames(tpm.HIGH_PGC))
tpm.sub.pgc <- tpm.pgc[list_names_pgc,]
tpm.sub.pgc <- na.omit(tpm.sub.pgc)

## PGC only early
tpm_early <- tpm[,c(1,2,5,6,9,10)]
tpm_early_Mean <- data.frame(row.names = rownames(tpm_early), s256 = rowMeans(tpm_early[,c(1,2)]),
                               High = rowMeans(tpm_early[,c(3,4)]), Dome = rowMeans(tpm_early[,c(5,6)]))

pheatmap(na.omit(tpm_early_Mean), color = hmcol, scale = 'row', show_rownames = F, cluster_cols = F)

## plot
library(RColorBrewer)
library(pheatmap)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

pheatmap(na.omit(log(heatmap_late + 1)), color = hmcol, show_rownames = F, cluster_cols = F)
pheatmap(na.omit(log(heatmap_early + 1)), color = hmcol, show_rownames = F, cluster_cols = F)

pheatmap(na.omit(heatmap_early), color = hmcol, scale = 'row', show_rownames = F, cluster_cols = F)


## plot only mean of replicates
tpm.256_PGC_sub <- subset(tpm.256_PGC, tpm.256_PGC$padj < 0.1)
tail_tpm <- tail(tpm.256_PGC_sub[order(tpm.256_PGC_sub$log2FoldChange),],20)
diff.expressed <- rownames(tail_tpm)
                    
heatmap_late_Mean <- data.frame(row.names = rownames(heatmap_late), 
                                sDomePGC = rowMeans(heatmap_late[,c('sDomePGC1', 'sDomePGC2')]),
                                sDomeSoma = rowMeans(heatmap_late[,c('sDomeSoma1', 'sDomeSoma2')]),
                                s10somitesPGC = rowMeans(heatmap_late[,c('s10somitesPGC1', 's10somitesPGC2')]),
                                s10somitesSoma = rowMeans(heatmap_late[,c('s10somitesSoma1', 's10somitesSoma2')]),
                                sPrim5PGC = rowMeans(heatmap_late[,c('sPrim5PGC1', 'sPrim5PGC2')]),
                                sPrim5Soma = rowMeans(heatmap_late[,c('sPrim5Soma1', 'sPrim5Soma2')])
                                )

pheatmap(heatmap_late_Mean, color = hmcol, scale = 'row', show_rownames = F, cluster_cols = F)



## cherry pick genes
GermDevGenes <- heatmap_early[diff.expressed,]

GermDevGenes <- heatmap_early[c('ENSDARG00000075113', 'ENSDARG00000078754', 'ENSDARG00000044774', 'ENSDARG00000036214',
'ENSDARG00000022813', 'ENSDARG00000014373', 'ENSDARG00000045139'),]
                             
#assign gene names from biomart
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset='drerio_gene_ensembl')
mart <- useDataset("drerio_gene_ensembl", useMart("ensembl"), verbose = FALSE)

s4 <- rownames(GermDevGenes)
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=s4,mart= mart)
row.names(G_list) <- G_list$ensembl_gene_id
G_list$ensembl_gene_id <- NULL
GermDevGenes <- merge(GermDevGenes,G_list, by = 'row.names')
row.names(GermDevGenes) <- GermDevGenes$external_gene_name
GermDevGenes$Row.names <- NULL
GermDevGenes$external_gene_name <- NULL

pheatmap(log(GermDevGenes[,c(1,5,9)]+1), color = hmcol, cluster_cols = F)
