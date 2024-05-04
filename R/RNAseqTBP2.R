## TBP2
rm(list = ls())

library(DESeq2)
setwd('~/Desktop/Postdoc/Data_Analysis/RNAseq/TBP2_counts/')
files_2 <- list.files(pattern = '\\.tab')
## creating a function which reads a csv file and make a matrix
readTagsPerGene <- function(x){
  out <- read.csv(x, skip = 4, header = FALSE, sep = "\t", row.names = 1,
                  col.names = c("","totalCount", "forwardCount", "reverseCount"))
  out
}

allReadsCounts <- lapply(files_2, readTagsPerGene)
# keeps only the total count column
totalReadsTBP2_all <- sapply(allReadsCounts, function(x) x$totalCount)
rownames(totalReadsTBP2_all) <- row.names(allReadsCounts[[1]])

## get column names
spl_names <- strsplit(files_2, "\\_Z")
col_names <- as.character(lapply(spl_names, function(u) u[[1]][1]))
colnames(totalReadsTBP2_all) <- col_names

## use old mm replicates
totalReadsTBP2_tdrd7 <- totalReadsTBP2_all[,c(5,6,3,4)]
## use new mm replicates
totalReadsTBP2 <- totalReadsTBP2_all[,c(1:4)]

condition_TBP2 <- rep(c('5mm', 'MO'), each=2)

run.deseq.genes2 <- function(x, condition= condition){
  coldata <- data.frame(condition = condition)
  row.names(coldata) <- colnames(x)
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  #res.sub <- subset(res, res$padj < 0.1)
  return(dds)
}

dds_TBP2 <- run.deseq.genes2(totalReadsTBP2, condition = condition_TBP2)
dds_TBP2_tdrd7 <- run.deseq.genes2(totalReadsTBP2_tdrd7, condition = condition_TBP2)
## correlation heatmap
library(pheatmap)
pheatmap(cor(totalReadsTBP2))

## filter low and high count genes
dds.original <- dds_TBP2
dds <- dds_TBP2[rowSums(counts(dds_TBP2)) > 1, ]
dds <- dds_TBP2[rowSums(counts(dds_TBP2)) < 40000, ]
nrow(dds)

plot(density(as.data.frame(counts(dds))$TBP2_5mm_PGC_rep1), xlim = c(0, 1500),
     ylim = c(0, .01), main = "Filtered Densities", col = "Blue")
lines(density(as.data.frame(counts(dds))$TBP2_5mm_PGC_rep2), col = "Red")
lines(density(as.data.frame(counts(dds))$TBP2_MO_PGC_rep1), col = "Orange")
lines(density(as.data.frame(counts(dds))$TBP2_MO_PGC_rep2), col = "Green")

plot(density(as.data.frame(counts(dds.original))$TBP2_5mm_PGC_rep1), xlim = c(0, 2500), ylim = c(0, .01), main = "Filtered Densities", col = "Blue")
lines(density(as.data.frame(counts(dds.original))$TBP2_5mm_PGC_rep2), col = "Red")
lines(density(as.data.frame(counts(dds.original))$TBP2_MO_PGC_rep1), col = "Orange")
lines(density(as.data.frame(counts(dds.original))$TBP2_MO_PGC_rep2), col = "Green")


## calculate fold change
res_TBP2 <- results(dds)
res_TBP2_tdrd7 <- results(dds_TBP2_tdrd7)
norm_data <- counts(dds,normalized=TRUE)
# save
write.table(res_TBP2, "~/Desktop/Postdoc/Data_Analysis/Tables/PGC_Analysis/Res_TBP2_newMM.txt", sep = '\t')
## volcano plot
source('~/Desktop/Postdoc/R_scripts/R graphics scripts/volcanoPlot.R')
custom.volcano.plot(res_TBP2, col1 = 'dodgerblue3', col2 = 'firebrick1')

## effect of variance linear transformation
ntd <- normTransform(dds_TBP2) # f(count(dds,normalized=TRUE) + pseudocount)
vsd <- vst(dds_TBP2, blind=FALSE)
rld <- rlog(dds_TBP2, blind=FALSE)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

## correlation
library(pheatmap)
pheatmap(cor(totalReadsTBP2))
## PCA plot
vsd <- readRDS('~/Desktop/Postdoc/Data_Analysis/RNAseq/TBP2_counts/vsdAllSamplesTBP2-Tdrd7.rds')
vsd <- vsd[,c(9:14, 17:19,21, 23:25)]
pcaData <- plotPCA(vsd, intgroup=c("condition", "sample", "treatment"), returnData = TRUE)
## plot PC3
source('~/Desktop/Postdoc/R_scripts/RNAseq scripts/PCA_componentsPC3.R')
plotPCA.comp(vsd, intgroup=c("condition", "sample", "treatment"))


## subset genes
subset_res <- function(x){
  x$degenes <- 'NS'
  x$degenes[x$log2FoldChange < -2 & x$padj < 0.05] <- 'TBP2_up'
  x$degenes[x$log2FoldChange > 2 & x$padj < 0.05] <- 'TBP2_down'
  return(x)
}
res_TBP2 <- subset_res(res_TBP2)
res_TBP2_tdrd7_2 <- subset_res(res_TBP2_tdrd7)


write.table(res_TBP2, '../../Tables/ResTBP2_5mmVsMO.txt', sep = '\t')

## gene ontology
library(clusterProfiler)
library(org.Dr.eg.db)

ego_down <- enrichGO(rownames(subset(res_TBP2, res_TBP2$degenes == 'TBP2_down')), 
                OrgDb = org.Dr.eg.db, keyType = 'ENSEMBL', ont = 'BP', pvalueCutoff = 0.1, 
                universe = rownames(totalReadsTBP2))
barplot(ego_down, showCategory = 25, font.size = 8)
ego_up <- enrichGO(rownames(subset(res_TBP2, res_TBP2$degenes == 'TBP2_up')), 
                     OrgDb = org.Dr.eg.db, keyType = 'ENSEMBL', ont = 'BP', pvalueCutoff = 0.1, 
                     universe = rownames(totalReadsTBP2))
barplot(ego_up, showCategory = 25, font.size = 8)


######
## correlation with upset
get_go_ID <- function(x, subset = subset){
  ego <- enrichGO(rownames(subset(x, x$degenes == subset)), 
                       OrgDb = org.Dr.eg.db, keyType = 'ENSEMBL', ont = 'BP', pvalueCutoff = 0.1, 
                       universe = rownames(totalReadsTBP2))
  f <- head(ego@result$ID, 150)
  return(f)
}
ego_down_TBP2 <- get_go_ID(res_TBP2, 'TBP2_down')
ego_up_TBP2 <- get_go_ID(res_TBP2, 'TBP2_up')
ego_down_tdrd7 <- get_go_ID(res_TBP2_tdrd7, 'TBP2_down')
ego_up_tdrd7 <- get_go_ID(res_TBP2_tdrd7, 'TBP2_up')

get_ENSEMBL_ID <- function(x, subset = subset){
  s <- rownames(subset(x, x$degenes == subset))
  return(s)
}

ens_down_TBP2 <- get_ENSEMBL_ID(res_TBP2, 'TBP2_down')
ens_up_TBP2 <- get_ENSEMBL_ID(res_TBP2, 'TBP2_up')
ens_down_tdrd7 <- get_ENSEMBL_ID(res_TBP2_tdrd7, 'TBP2_down')
ens_up_tdrd7 <- get_ENSEMBL_ID(res_TBP2_tdrd7, 'TBP2_up')

## plot
library(UpSetR)
library(wesanderson)
library(RColorBrewer)

pal= brewer.pal(11, "Set3")
pal2=wes_palette("GrandBudapest2")

# or if you have numbers just run
input_go <- c(TBP2_up=length(ego_up_TBP2), TBP2_down=length(ego_down_TBP2), Tdrd7_up=length(ego_up_tdrd7),
          Tdrd7_down=length(ego_down_tdrd7), 
          "TBP2_up&Tdrd7_up"= length(intersect(ego_up_TBP2, ego_up_tdrd7)),
          "TBP2_down&Tdrd7_down"=length(intersect(ego_down_TBP2, ego_down_tdrd7)),
          "TBP2_up&Tdrd7_down"=length(intersect(ego_up_TBP2, ego_down_tdrd7)),
          "TBP2_down&Tdrd7_up"=length(intersect(ego_down_TBP2, ego_up_tdrd7)))

input<- c(TBP2_up=length(ens_up_TBP2), TBP2_down=length(ens_down_TBP2), Tdrd7_up=length(ens_up_tdrd7),
          Tdrd7_down=length(ens_down_tdrd7), 
          "TBP2_up&Tdrd7_up"= length(intersect(ens_up_TBP2, ens_up_tdrd7)),
          "TBP2_down&Tdrd7_down"=length(intersect(ens_down_TBP2, ens_down_tdrd7)),
          "TBP2_up&Tdrd7_down"=length(intersect(ens_up_TBP2, ens_down_tdrd7)),
          "TBP2_down&Tdrd7_up"=length(intersect(ens_down_TBP2, ens_up_tdrd7)))
par(family="serif")
##!!check numbers in set bar!!
p<- upset(fromExpression(input_go), 
          nintersects = 14, 
          nsets = 4, 
          order.by = "freq", 
          decreasing = T, 
          mb.ratio = c(0.6, 0.4),
          number.angles = 0, 
          text.scale = 2,
          point.size = 2.8, 
          line.size = 2,
          sets.bar.color= pal2,
          matrix.color = "black")
p

#### cross correlation logfold change
log_fold_cor <- merge(as.data.frame(res_TBP2), as.data.frame(res_TBP2_tdrd7), by = 0)
ggplot(log_fold_cor,aes(x=log_fold_cor$log2FoldChange.x, y=log_fold_cor$log2FoldChange.y)) +
  geom_point(aes(alpha=0.05, color = log_fold_cor$genes)) + theme_bw() +
  geom_vline(xintercept=0, linetype="longdash", colour="red", size=0.4) +
  #yintercept is the -log10(pvalue): if pad < 0.05 -> -log10(0.05) = 1.3
  geom_hline(yintercept=0, linetype="longdash", colour="red", size=0.4)

##only significant
log_fold_cor_sig <- log_fold_cor[!log_fold_cor$padj.x < 0.05 | !log_fold_cor$padj.y < 0.05,]
ggplot(log_fold_cor_sig,aes(x=log_fold_cor_sig$log2FoldChange.x, y=log_fold_cor_sig$log2FoldChange.y)) +
  geom_point(aes(alpha=0.05, color = log_fold_cor_sig$genes)) + theme_bw() +
  geom_vline(xintercept=0, linetype="longdash", colour="red", size=0.4) +
  #yintercept is the -log10(pvalue): if pad < 0.05 -> -log10(0.05) = 1.3
  geom_hline(yintercept=0, linetype="longdash", colour="red", size=0.4)


###
## GENES EXPRESSED IN OOCYTE AND LATE PGCS
setwd('~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables/')
tpm <- read.csv('tpmStagesFiltered.txt', sep = '\t')
upreg.insoma.late <- read.csv('UpregulatedInSomavsPGCatPrim5.txt', sep = '\t', row.names = 1)
upreg.inpgcs.late <- read.csv('UpregulatedGenesInPGCvsSomaPrim5.txt', sep = '\t', row.names = 1)

tpm$EarlyMean <- rowMeans(tpm[c('s256PGC1', 's256PGC2', 's256Soma1', 's256Soma2')],)
#tpm_oocyte <- subset(tpm, tpm$EarlyMean > 20)
## no maternal
tpm_oocyte <- subset(tpm, tpm$EarlyMean < 20)

pgc.threshold.bigger.than.5 <- subset(upreg.inpgcs.late, upreg.inpgcs.late$log2FoldChange < -5)
soma.threshold.bigger.than.5 <- subset(upreg.insoma.late, upreg.insoma.late$log2FoldChange > 5)

oocyte_and_pgc <- merge(tpm_oocyte, pgc.threshold.bigger.than.5, by = 0)

oocyte_and_pgc_names <- oocyte_and_pgc$Row.names
## get oocyte and pgc genes
d <- na.omit(as.data.frame(res_TBP2)[oocyte_and_pgc_names,])
p <- na.omit(as.data.frame(res_TBP2)[rownames(pgc.threshold.bigger.than.5),])
s <- na.omit(as.data.frame(res_TBP2)[rownames(soma.threshold.bigger.than.5),])
  
overlap.histograms <- function(x,y){
  data.frame1 <- data.frame(values = x, hist = 'hist1')
  data.frame2 <- data.frame(values = y, hist = 'hist2')
  
  data.frame.final <- rbind(data.frame1, data.frame2)
  return(data.frame.final)
}

data.frame.final <- overlap.histograms(res_TBP2, s)
data.frame.final.pgc <- overlap.histograms(res_TBP2, p)

ggplot(data.frame.final, aes(values.log2FoldChange, fill = hist)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') +
  geom_density(aes(color = hist), alpha = .1) +
  scale_fill_manual(values = c('blue', 'red')) +
  scale_color_manual(values = c('blue', 'red')) +
  theme_bw()

ggplot(data.frame.final.pgc, aes(values.log2FoldChange, fill = hist)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') +
  geom_density(aes(color = hist), alpha = .1) +
  scale_fill_manual(values = c('blue', 'lightgreen')) +
  scale_color_manual(values = c('blue', 'lightgreen')) +
  theme_bw()
## get CAGE peaks
library(CAGEr)
library(GenomicRanges)
library(ChIPpeakAnno)
library(ChIPseeker)
library(AnnotationDbi)
## load cage set
myCAGEset <- readRDS("~/Desktop/Postdoc/Data_Analysis/CAGEsets/CAGEset_PGC_soma_Early_Late.rds")
## get tag cluster for PGC prim5
sample.PGC <- unname(sampleLabels(myCAGEset)[3])
pgc_tag_cluster <- tagClusters(myCAGEset, sample = sample.PGC, 
                               returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
## annotate peaks
source('~/Desktop/Postdoc/R_scripts/gcCAGE/R/Functions.R')
txdb <- loadDb('~/Desktop/Postdoc/Data_Analysis/annotation/txdb_DanRer7.sqlite')
pgc_tag_cluster_anno <- annotate.cage.peaks(pgc_tag_cluster)

pgc_tag_cluster_anno_merged <- merge(pgc_tag_cluster_anno, data.frame(res_TBP2), by.x = 'geneId', by.y = 0)
pgc_tag_cluster_anno_merged$geneId <- NULL
## resize promoter window for heatmap and order by interquartile width
resize_prom <- function(x){
  x <- data.frame(x)
  y <- GRanges(seqnames = x$seqnames, IRanges(start = x$dominant_ctss, end = x$dominant_ctss), strand = x$strand)
  values(y) <- x[,6:ncol(x)]
  
  resize.tc <- resize(y, 150, fix = 'center')
  resize.tc <- resize.tc[order(resize.tc$degenes),]
  return(resize.tc)
}

prom_for_heatmap <- resize_prom(pgc_tag_cluster_anno_merged)

## plot cage signal
library(heatmaps)

## plot nucleotide frequencies
library(BSgenome.Drerio.UCSC.danRer7)
.centerTSS <- function(x, subset = 'NS'){
  new.data.frame <- data.frame(x)
  new.data.frame <- subset(new.data.frame, new.data.frame$degenes == subset)
  zebrafishPromotersTSS<-GRanges(seqnames = new.data.frame$seqnames,
                                 ranges=IRanges(start = new.data.frame$dominant_ctss, 
                                                end = new.data.frame$dominant_ctss),
                                 strand = new.data.frame$strand,
                                 IQ_width = new.data.frame$interquantile_width,
                                 DEgene = new.data.frame$degenes,
                                 seqlengths = seqlengths(Drerio))

  zebrafishPromotersTSSflank <- promoters(zebrafishPromotersTSS, upstream = 150,
                                        downstream = 150)
  zebrafishPromotersTSSflankSeq <- getSeq(Drerio, zebrafishPromotersTSSflank)
  return(zebrafishPromotersTSSflankSeq)
}
zebrafishPromotersTBP2up <- .centerTSS(prom_for_heatmap, subset = 'TBP2_up')
plotPatternDensityMap(zebrafishPromotersTBP2up,
                                   "WW",
                                   flankUp = 150, flankDown = 150)
zebrafishPromotersTBP2down <- .centerTSS(prom_for_heatmap, subset = 'TBP2_down')
plotPatternDensityMap(zebrafishPromotersTBP2down,
                      "WW",
                      flankUp = 150, flankDown = 150)#, bandWidth = 1)


up <- zebrafishPromotersTSSflank$DEgene == 'TBP2_up'
down <- zebrafishPromotersTSSflank$DEgene == 'TBP2_down'

plotPatternOccurrenceAverage(regionsSeq = zebrafishPromotersTBP2up,
                            
                               patterns = c("WW", "SS"), flankUp = 150, flankDown = 150,
                             smoothingWindow = 3, color = c("red3", "blue3"), cex.axis = 0.9)

### Chromatin
library(rtracklayer)
setwd('~/Desktop/Postdoc/Data_Analysis/ATACseq/')

atac.fold.change.pgc <- readRDS('ATAC_FoldChange_Cov120.rds')
  #readRDS('ATAC_FoldChange_PGC_5mm_rep2_resize20.rds')
atac.fold.change.soma <- readRDS('ATAC_FoldChange_Soma_MO_rep1_resize20.rds')
## summit to center the peaks
file_summit <- list.files()[grep('summit', list.files())]

.import.summit <- function(f_summit){
  summit <- read.csv(f_summit, sep = '\t', header = F)[,c(1:3)]
  colnames(summit) <- c('seqnames', 'start', 'end')
  s_range <- GRanges(seqnames = summit$seqnames, IRanges(start = summit$start, end = summit$start))
  return(s_range)
}

pgc.summit <- .import.summit(file_summit[1])
soma.summit <- .import.summit (file_summit[2])

## get the tss
## extend 5kb up and downstream
.centerTSS2 <- function(x, width = 1000, class = class){
  new.data.frame <- data.frame(x)
  new.data.frame <- subset(new.data.frame, new.data.frame$degene == class)
  new.data.frame <- na.omit(new.data.frame)
  
  if(class == 'TBP2_up'){
    new.data.frame <- head(new.data.frame[order(new.data.frame$log2FoldChange),], 200)
  }
  else{
    new.data.frame <- tail(new.data.frame[order(new.data.frame$log2FoldChange),], 200)
  }
  zebrafishPromotersTSS<-GRanges(seqnames = new.data.frame$seqnames,
                                 ranges=IRanges(start = new.data.frame$dominant_ctss, 
                                                end = new.data.frame$dominant_ctss),
                                 strand = new.data.frame$strand,
                                 IQ_width = new.data.frame$interquantile_width,
                                 DEgene = new.data.frame$degenes,
                                 seqlengths = seqlengths(Drerio))
  TSS <- resize(zebrafishPromotersTSS, width = width, fix = 'end')
  
  return(TSS)
}
TSS_5kb_up <- .centerTSS2(pgc_tag_cluster_anno_merged, width = 10000, class = 'TBP2_up')
TSS_500_up <- .centerTSS2(pgc_tag_cluster_anno_merged, width = 1000, class = 'TBP2_up')

TSS_5kb_down <- .centerTSS2(pgc_tag_cluster_anno_merged, width = 10000, class = 'TBP2_down')
TSS_500_down <- .centerTSS2(pgc_tag_cluster_anno_merged, width = 1000, class = 'TBP2_down')

## get atac peaks  and make grange with these as start
.clean.peaks <- function(tss5k, tss500){
  centered_atac_peaks <- subsetByOverlaps(pgc.summit, tss5k)
  ## remove promoters
  non_promoters_atac_peaks <- subsetByOverlaps(centered_atac_peaks, tss500, invert = TRUE)
  
  #non_promoters_atac_peaks <- subset(non_promoters_atac_peaks, non_promoters_atac_peaks$score > 0)
  #non_promoters_atac_peaks <- reduce(non_promoters_atac_peaks, min.gapwidth = 100)
  non_promoters_atac_peaks <- resize(non_promoters_atac_peaks, width = 1500, fix = 'center')
}
cleaned_peak_up <- .clean.peaks(TSS_5kb_up, TSS_500_up)
cleaned_peak_down <- .clean.peaks(TSS_5kb_down, TSS_500_down)

## cluster peaks

## plot heatmaps
library(genomation)
m_TBP2_up <- ScoreMatrix(target = atac.fold.change.pgc, windows = cleaned_peak_up, weight.col = 'score', strand.aware = T)
m_TBP2_down <- ScoreMatrix(target = atac.fold.change.pgc, windows = cleaned_peak_down, weight.col = 'score', strand.aware = T)

heatMatrix(m_TBP2_up, winsorize = c(0,95), clustfun = function(x) kmeans(x, centers=2)$cluster, 
           col = brewer.pal(n = 5, "Blues"))
## plot meta
df <- data.frame(cbind(colMeans(m_TBP2_up), colMeans(m_TBP2_down)))
colnames(df) <- c("PGCup", "PGCdown")

df.gg <- df %>% gather(Expression, Tn5CutSites, PGCup:PGCdown)
df.gg$Index <- rep(c(-750:-1, 1:750), 2)

ggplot(df.gg, aes(x=Index, y=Tn5CutSites, Group=factor(Expression))) +
  geom_line(aes(colour=factor(Expression)), size = 1.5) +
  scale_colour_manual(values = c("seagreen4", "slateblue")) + ###7CAE00 ###F8766D +
  geom_vline(xintercept=0, colour="grey", linetype="longdash") +
  theme_classic() + theme(legend.position='none')
## random assign strand
n <- dim(as.data.frame(non_promoters_atac_peaks))[1]
s_strand <- sample(c('-','+'), replace=TRUE, size=n)
strand(non_promoters_atac_peaks) <- s_strand



### exclude genes with low count
## work only on upreg genes


#######
## analysis of pentamers
library(seqPattern)
library(purrr)
library(dplyr)
library(GenomicFeatures)
library(ChIPpeakAnno)
library(AnnotationDbi)

## load cage data
myCAGEset <- readRDS("~/Desktop/Postdoc/Data_Analysis/CAGEsets/CAGEset_PGC_soma_Early_Late.rds")

## load diff expressed genes and get the promoter coordinates
TBP2_dataframe <- data.frame(res_TBP2)
TBP2_dataframe$degenes <- 'NS'
TBP2_dataframe$degenes[TBP2_dataframe$log2FoldChange < -1 & TBP2_dataframe$padj < 0.05] <- 'TBP2_up'
TBP2_dataframe$degenes[TBP2_dataframe$log2FoldChange > 1 & TBP2_dataframe$padj < 0.05] <- 'TBP2_down'
## ENSEMBL ID in column
TBP2_dataframe$ENSEMBL <- rownames(TBP2_dataframe)

## use high and soma prim promoters
sample.names <- unname(sampleLabels(myCAGEset)[3])

txdb <- toGRanges(loadDb('~/Desktop/Postdoc/Data_Analysis/annotation/txdb_DanRer7.sqlite'))

## merge deseq table with gene coordinates
## matches DGE class and logFoldChange with gene coordinates
deseq.coordinates <- merge(as.data.frame(txdb), TBP2_dataframe, by = 0) %>%
  mutate(cc_id = paste(seqnames,start,end,width,strand,sep = "_"))
deseq.res.gr <- GRanges(seqnames = deseq.coordinates$seqnames,
                        ranges = IRanges(start = deseq.coordinates$start,
                                         end = deseq.coordinates$end),
                        strand = deseq.coordinates$strand)
values(deseq.res.gr) <- deseq.coordinates[,7:ncol(deseq.coordinates)]
## match tag clusters with DEG
tcs <- lapply(sample.names, function(x){
  
  a <- tagClusters(myCAGEset, sample = x, returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
  a <- a[a$chr %in% paste("chr",1:25, sep = ""),]
  a$sampleID <- x
  a <- mutate(a, promoterID = paste(chr,start,end,strand, sep = "_"))
  gr <- GRanges(seqnames = a$chr,
                ranges = IRanges(start = a$start,
                                 end = a$end),
                strand = a$strand, iq = a$interquantile_width)
  values(gr) <- a[,c(6:12,ncol(a))]
  return(gr)
})
names(tcs) <- sample.names

### merge promoter tss with DEG class and logFoldChange
tcs_res <- purrr::map(tcs, function(x,y){
  # determine all overlaps and have granges map to TC coordinates
  hits <- findOverlaps(x,y)
  x_sub <- x[queryHits(hits)]
  y_sub <- y[subjectHits(hits)]
  tc_deseq <- x_sub
  values(tc_deseq) <- cbind(values(tc_deseq), values(y_sub))
  # keep only highest:
  df <- tbl_df(data.frame(promoterID = tc_deseq$promoterID,
                          tpm = tc_deseq$tpm,
                          cc_id = tc_deseq$cc_id)) %>%
    group_by(cc_id) %>%
    dplyr::slice(which.max(tpm))
  # select only these promoterIDs
  sel <- which(tc_deseq$promoterID %in% df$promoterID & tc_deseq$cc_id %in% df$cc_id)	
  return(tc_deseq[sel])
}, y = deseq.res.gr)

saveRDS(tcs_res, '../../Tables/TSSofUp_DownregulatedGenesInTBP2.rds')
## center to sample TCs to dominant TSS and extend 120 bp up and 50 downstream  
# function
center_on_tss <- function(x, up, down){
  y <- x
  ranges(y) <- IRanges(start = x$dominant_ctss, end = x$dominant_ctss)
  y <- trim(promoters(y, upstream = range_extend[1], downstream = range_extend[2]))
  return(y)}
# execution
range_extend <- c(120,50)
tc_tss <- map(tcs_res, center_on_tss, up = range_extend[1], down = range_extend[2])
names(tc_tss) <- sample.names

## generate random pentamers with A and Ts
library(tcR)
pents <- generate.kmers(.k = 5, .alphabet = c("A","T"))

## Get the numbers of matches for each pentamer
# prerequisites
library(BSgenome.Drerio.UCSC.danRer7)
bs = BSgenome.Drerio.UCSC.danRer7
# list of long dataframes
pent_list <- list()
for(i in seq_along(tc_tss)){
  sample_name <- names(tc_tss)[i]
  print(paste0("Analysing sample: ",sample_name, "..."))
  tc_gr <- tc_tss[[i]]
  # get seqs
  tc.seq <- getSeq(bs,tc_gr)
  wwbox <- lapply(pents, function(x){
    message(paste0("Pattern matching to: ", x, ".."))
    ## calculate 100% match  score and sum for hit 
    pat <- getPatternOccurrenceList(regionsSeq = tc.seq, patterns = x)[[1]]
    # create empty matrix and fill in to select the right range
    mat <- matrix(data = 0,
                  nrow = length(tc.seq),
                  ncol = width(tc.seq)[1])
    
    mat[cbind(pat$sequence, pat$position)] <- pat$value
    colnames(mat) <- paste0("pos_",1:ncol(mat))
    rownames(mat) <- tc_gr$cc_id
    ## count occurences between -40 and -20 
    # window to count (-120+80 = -40)
    sel = 80:100
    motif.df <- tbl_df(mat)[,sel]
    motif.df$cc_id <- rownames(mat)
    # bind the two frames
    meta <- tbl_df(tc_gr) %>% mutate(sampleID = sample_name)
    info_motif <- inner_join(meta[,c("ENSEMBL","degenes","tpm","cc_id", "log2FoldChange")], motif.df, by = "cc_id")
    # find the max per row (exclude first 3 columns)
    one <- sapply(1:nrow(info_motif),function(y){ max(info_motif[y,-c(1:5)])})
    # replace other columns
    df <- cbind(info_motif[,1:5],tibble(motif = x,motif_value = one))
    # return
    return(df)
  })
  pent_list[[i]] <- do.call(rbind, wwbox)
}
# all together
penties <- do.call(rbind, pent_list)

##bar plot
library(tidyr)
library(ggplot2)
## divided by the number of promoters which have that occurrence
bar_tab <- penties %>%
  group_by(degenes,motif) %>%
  summarize(N = n(),
            motif_rel = sum(motif_value != 0)/N)
# only three samples:
bar_gg <- bar_tab[(bar_tab$degenes == "TBP2_up") | (bar_tab$degenes == "TBP2_down"),]
# ordering 
or<- bar_gg %>%
  ungroup() %>%
  dplyr::select(degenes,motif,motif_rel) %>%
  spread(degenes,motif_rel) %>%
  group_by(motif) %>% 
  mutate(diff = max(abs(diff(c(TBP2_up, TBP2_down))))) %>%
  arrange(desc(diff))
bar_gg$motif <- factor(bar_gg$motif, levels = rev(or$motif))
## normalize by size
dim(subset(TBP2_dataframe[TBP2_dataframe$degenes == 'TBP2_down',]))
#1005
dim(subset(TBP2_dataframe[TBP2_dataframe$degenes == 'TBP2_up',]))
#2044
bar_gg$N[bar_gg$degenes == 'TBP2_down'] <- bar_gg$N[bar_gg$degenes == 'TBP2_down']/1005
bar_gg$N[bar_gg$degenes == 'TBP2_up'] <- bar_gg$N[bar_gg$degenes == 'TBP2_up']/2044

##plot
library(viridis)
vallie <- sample(viridis(31))
# orders
bar_gg$degenes <- factor(bar_gg$degenes, levels = c("TBP2_down","TBP2_up"))
p <- ggplot(bar_gg,aes(x = motif, y = motif_rel, fill = degenes)) +
  geom_col(position = 'dodge') +
  theme_classic() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = vallie) +
  coord_flip()
p

###### distinguish canonical vs non canonical
### divide 
penties$TATAtype <- 'non_canonical'
penties$TATAtype[penties$motif == 'TATAA'] <- 'canonical'
penties$TATAtype[penties$motif == 'ATATA'] <- 'canonical'
penties$TATAtype[penties$motif == 'ATAAA'] <- 'canonical'
penties$TATAtype[penties$motif == 'TAAAA'] <- 'canonical'

## only matches
penties <- subset(penties, penties$motif_value > 0)
### merge with row counts
d <- merge(totalReadsTBP2, penties, by.x = 0, by.y = 'ENSEMBL')

## proportion of deg with tata boxes
d_can <- subset(d,d$TATAtype == 'canonical') %>% group_by(degenes) %>% summarise(N = n())

d <- d[order(d$log2FoldChange),]
d_up <- head(d, 300)
d_down <- tail(d, 300)
d1 <- bind_rows(d_up, d_down) 
## plot
ggplot(d1, aes(x=log(TBP2_5mm_PGC_rep1+1), y=log(TBP2_MO_PGC_rep1+1))) +
  geom_point(aes(colour=factor(TATAtype)),size=2) + xlim(c(-2,12)) + ylim(c(-2,9)) +
  scale_color_manual(values = c("lightseagreen", "darkorange1")) +
  theme_bw()

## 
d_sign <- subset(d, d$degenes == 'TBP2_up' | d$degenes == 'TBP2_down')
ggplot(d_sign, aes(x=log(TBP2_5mm_PGC_rep1), y=log(TBP2_MO_PGC_rep1))) +
  geom_point(aes(colour=factor(TATAtype)),size=2) + xlim(c(0,12)) + ylim(c(0,9))

d2 <- d %>% group_by(degenes) %>% group_by(TATAtype) %>% summarise(N = n())

##
## initiators
setwd("~/Desktop/Postdoc/Data_Analysis/CAGEsets/")
danRer7CAGEset <- readRDS("CAGEsetPGCsomaEarlyLateMerged.replicates.rsd")

## generate a list of tag clusters
sample.PGC <- unname(sampleLabels(danRer7CAGEset)[c(3)])
cage.peaks <- lapply(sample.PGC, function(x){ tagClusters(danRer7CAGEset, sample = x, returnInterquantileWidth = TRUE,  
                                                       qLow = 0.1, qUp = 0.9)})
tbp2 <- as.data.frame(subset(res_TBP2, res_TBP2$degenes == 'TBP2_up'))
notbp2 <- as.data.frame(subset(res_TBP2, res_TBP2$degenes == 'TBP2_down'))
ns <- as.data.frame(subset(res_TBP2, res_TBP2$degenes == 'NS'))

tc.list <- list(tbp2, notbp2, ns)
## annotate
source('~/Desktop/Postdoc/R_scripts/gcCAGE/R/Functions.R')
txdb <- loadDb('~/Desktop/Postdoc/Data_Analysis/annotation/txdb_DanRer7.sqlite')
pgc_tag_cluster_anno <- lapply(cage.peaks, annotate.cage.peaks)
pgc_tag_cluster_anno <- data.frame(pgc_tag_cluster_anno)
## merge with RNAseq
func <- function(x,y){merge(x, y, by.x=0, by.y='geneId')}
pgc_tag_cluster_anno_merged <- lapply(tc.list, func, pgc_tag_cluster_anno)

# write function
.getSequences <- function(x, range = c(1,1), SeqLengths, 
                          bs_genome, remove.na = FALSE){
  x <- x[sample(1:nrow(x), nrow(x)),] # rewrite x by shuffling sequences 
  # into granges from dominant
  tc.gr <- GRanges(seqnames = x$seqnames, ranges = IRanges(start = x$dominant_ctss, end = x$dominant_ctss), 
                   strand = x$strand, interquantile_width = x$interquantile_width, seqlengths = SeqLengths)
  # extend grange to the flanking regions
  tc.flank <- trim(promoters(tc.gr, upstream = range[1], downstream = range[2])) ## 1 up and 1 down
  tc.flank <- resize(tc.flank, 3, fix = 'center')
  tc.flank <- resize(tc.flank, 2, fix = 'end')
  selection <- width(tc.flank) == sum(range) ## only those with width = 2bp
  tc.flank <- tc.flank[selection]
  tc.flank <- tc.flank[order(tc.flank@elementMetadata$interquantile_width)] ## order by IQ width
  # get sequence
  tc.seq <- getSeq(bs_genome, tc.flank)
  #if(remove.na == TRUE){
  #	tc.seq <- clean(tc.seq)
  #}
  return(tc.seq)
}
library(BSgenome.Drerio.UCSC.danRer7)
# get seqs for 2bp around dominant peaks: range = c(1,1)
seqs.tc <- lapply(pgc_tag_cluster_anno_merged, .getSequences, range = c(1,1), 
                  SeqLengths = seqlengths(BSgenome.Drerio.UCSC.danRer7), 
                  bs_genome = BSgenome.Drerio.UCSC.danRer7, remove.na = TRUE)
# determine proportions
seqs.table <- lapply(seqs.tc, function(x){as.data.frame(table(x))}) ## table counts each combination
names(seqs.table) <- c('TBP2 +', 'TBP2 -', 'NS')
seqs.prop <- lapply(seqs.table, function(x){x$prop <- (x$Freq/sum(x$Freq))*100;return(x)})
# file for plotting
initiators <- do.call("rbind",mapply(cbind, seqs.prop, "SampleID"=names(seqs.table), SIMPLIFY=F))
rownames(initiators) <- 1:nrow(initiators) 
## remove NN
#grep('N', initiators$x)-> nn
#initiators <- initiators[-nn,]

## Plots
# all in one plot:
col <- c('firebrick1', 'dodgerblue3', 'grey')
p <- ggplot(initiators, aes(x = x, y = prop, fill = SampleID)) +
  scale_fill_manual(values = col) + geom_col(position = "dodge") +  theme_bw()
p


