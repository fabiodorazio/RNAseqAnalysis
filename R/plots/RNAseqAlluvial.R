#for pieces of code
## for notes
### titles

### run DESeq on gene counts

library('DESeq2')
library('ggplot2')
library(dplyr)
library(alluvial)

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
  
  ## for alluvial plot
  res$stage <- paste0('s',input1)
  return(res)
}

#run the function to obtain differential genes and MA plot
condition <- rep(c('PGC', 'Soma'), times = 2)

tpm.dome <- run.deseq.genes('Dome', condition = condition)
tpm.10somites <- run.deseq.genes('10somites', condition = condition)
tpm.prim5 <- run.deseq.genes('Prim5', condition = condition)

list.stages <- list(tpm.dome, tpm.10somites, tpm.prim5)

## add column
list.stages <- lapply(list.stages, function(x) {data.frame(x)})
list.stages <- lapply(list.stages, function(x) {x$Class = 'NS'; return(x)})
list.stages <- lapply(list.stages, function(x) {x$Class[x$log2FoldChange < 0 & x$padj < 0.1] <- 'Hyper'; return(x)})
list.stages <- lapply(list.stages, function(x) {x$Class[x$log2FoldChange > 0 & x$padj < 0.1] <- 'Hypo'; return(x)})
list.stages <- lapply(list.stages, function(x) {colnames(x)[8] <- x$stage[1]; return(x)})

stages.final <- merge(list.stages[[1]], list.stages[[2]], by = 0)
rownames(stages.final) <- stages.final$Row.names
stages.final$Row.names <- NULL
stages.final <- merge(stages.final, list.stages[[3]], by = 0)

stages.final <- stages.final[,c(1,8,9,17,25)]

alluv_tb <- as_tibble(stages.final) %>%
  count(sDome, s10somites, sPrim5)
  

alluvial(alluv_tb %>% select(-n),
         freq=alluv_tb$n, border=NA, alpha = 0.5,
         col=case_when(alluv_tb$sPrim5 == "Hyper" ~ "red",
                       alluv_tb$sPrim5 == "Hypo" ~ "blue",
                       TRUE ~ "grey"),
         cex=0.75,
         axis_labels = c("dome", "10somites", 'prim5'))


########
ggplot(as.data.frame(UCBAdmissions),
       aes(y = Freq, axis1 = Gender, axis2 = Dept)) +
  geom_alluvium(aes(fill = Admit), width = 1/40) +
  geom_stratum(width = 1/40, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label =  Dept)) +
  scale_x_discrete(limits = c("Gender", "Dept"), expand = c(.03, .03)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("UC Berkeley admissions and rejections, by sex and department")




