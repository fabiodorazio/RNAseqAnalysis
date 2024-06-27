#@include run.deseq.genes.R
## RNAseq Tdrd7

## import count files
library(DESeq2)
setwd('../data/Tdrd7_MO/')
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

## remove bad replicates
library(pheatmap)
pheatmap(cor(totalReadsCount_2))
totalReadsCount_3 <- totalReadsCount_2[,-c(3:14,16,18,28,22,24,30:33)]
condition <- c(rep('Dnd', times=2), rep('5mm', times=4), rep('MO', times=4), rep("Rescue", times=2))
sample <- c(rep('PGC',times=2),rep(c(rep('PGC', times=2), rep('Soma', times=2)), times=2), rep('PGC', times=2))

## only new
only_new <- totalReadsCount_2[,-c(13,15,17:19,21,23:27)]
condition <- c(rep('Rescue', times=2), rep('wt', times=10), rep("5mm", times=2), rep("MO", times=2))
sample <- c(rep('PGC', times=2), rep(c('PGC','Soma'), times=5), rep("PGC", times=4))

##only new2
only_new_2 <- totalReadsCount_2[,c(3,4,14,16,19,21,25,27)]
condition <- c(rep('wt', times=2),rep("5mm", times=2), rep("MO", times=2), rep('Rescue', times=2))

#########
## DGE for rescue (MO2018, Rescue1/Rescue2)
rescue_MO <- totalReadsCount_2[,c(21,23,27,29)]
condition_MO <- c(rep('MO', times=2), rep("Rescue", times=2))
#(mm2018, Rescue1/Rescue2)
rescue_mm <- totalReadsCount_2[,c(15,17,27,29)]
condition_mm <- c(rep('5mm', times=2), rep("Rescue", times=2))
#########

## all samples
all_samples <- totalReadsCount_2[,-c(14,16,20,22)]
condition <- c(rep('Rescue', times=2), rep('wt', times=10), rep("5mm", times=4), rep("MO", times=4), rep('Rescue', times=3))
sample <- c(rep('PGC', times=2), rep(rep(c('PGC', 'Soma'), each=2),times =1), "PGC", "Soma", rep(rep(c('PGC', 'Soma'), each=2),times =3), rep('PGC', times=3))

run.deseq.genes <- function(x, condition= condition){
  coldata <- data.frame(sample = sample, treatment = treatment,condition = condition)
  row.names(coldata) <- colnames(x)
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ sample + condition)
  dds <- DESeq(dds)

  return(dds)
}

good_rep <- run.deseq.genes(totalReadsCount_3, condition = condition)
good_rep_res <- results(good_rep)
good_rep_res_sub <- subset(good_rep_res, good_rep_res$padj < 0.1)
#calculate variance for PCA plot 
vsd <- vst(good_rep, blind=TRUE)
plotPCA(vsd, intgroup=c("condition", "sample"))


## run deseq for sample pairs
run.deseq.genes2 <- function(x, condition= condition){
  coldata <- data.frame(condition = condition)
  row.names(coldata) <- colnames(x)
  dds <- DESeqDataSetFromMatrix(countData = x, colData = coldata, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds)
  #res.sub <- subset(res, res$padj < 0.1)
  return(dds)
  
}

rescue_mm_dds <- run.deseq.genes2(rescue_mm, condition = condition_mm)
rescue_MO_dds <- run.deseq.genes2(rescue_MO, condition = condition_MO)

## number of DEG
res_mm_res <- results(rescue_mm_dds)
res_mm_sub <- subset(res_mm_res, res_mm_res$padj < 0.1)

res_MO_res <- results(rescue_MO_dds)
res_MO_sub <- subset(res_MO_res, res_MO_res$padj < 0.1)

dim(res_MO_sub[res_MO_sub$log2FoldChange < 1])
## count DEG
library(dplyr)
res_mm_sub$DEG <- ifelse(res_mm_sub$log2FoldChange > 0, 'up', 'down')
count_mm <- data.frame(res_mm_sub) %>% group_by(DEG) %>% summarise(N= n())
res_MO_sub$DEG <- ifelse(res_MO_sub$log2FoldChange > 0, 'up', 'down')
count_MO <- data.frame(res_MO_sub) %>% group_by(DEG) %>% summarise(N= n())
bar_table <- data.frame(row.names = c('down', 'up'), mm = count_mm$N, MO = count_MO$N)
barplot(as.matrix(bar_table), col = c('darkblue', 'light blue'), ylim=c(0,3000))
## logFold correlation
log_fold_cor <- merge(as.data.frame(res_MO_res), as.data.frame(res_mm_res), by = 0)

smoothScatter(log_fold_cor$log2FoldChange.x, log_fold_cor$log2FoldChange.y, xlim = c(-10,10))
ggplot(log_fold_cor,aes(x=log_fold_cor$log2FoldChange.x, y=log_fold_cor$log2FoldChange.y)) +
  geom_point(aes(alpha=0.05, color = log_fold_cor$genes)) + theme_bw()
## import gos from GO:0019953, GO:0022414, GO:0048609, GO:0022412 (reproduction)
setwd('~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/')
germ_cell_go <- read.csv('GermcellDevGO.txt', sep = '\t')
germ_cell_go <- as.character(unique(germ_cell_go[,1]))
## import tpm matrix
tpm <- read.csv('../tpmStagesFiltered.txt', sep = '\t')
tpm_germ <- tpm[germ_cell_go,]
tpm_germ_sub <- subset(tpm_germ, tpm_germ$s256PGC1 < 1)
## add column to logfold corr
log_fold_cor$genes <- 'NS'
rownames(log_fold_cor) <- log_fold_cor$Row.names
log_fold_cor[rownames(tpm_germ_sub),]$genes <- 'germ cells'
log_fold_cor[rownames(up_MO),]$genes <- 'upMO'
log_fold_cor[rownames(down_MO),]$genes <- 'downMO'

library(ggplot2)
p1 <- log_fold_cor[log_fold_cor$genes == 'ge',]
p <- ggplot(log_fold_cor,aes(x=log_fold_cor$log2FoldChange.x, y=log_fold_cor$log2FoldChange.y, color = genes))
p + geom_point(alpha=0.7) + facet_wrap(~genes, ncol=4)

## correlation heatmap
library(pheatmap)
mat <- counts(dds, normalized=TRUE) #normalised read counts (raw/size factor)
ntd <- normTransform(dds_new)


## expression of cherry-picked genes
library(biomaRt)
library(pheatmap)
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
                             'ENSDARG00000045139', 'ENSDARG00000039310', 'ENSDARG00000027957', 'ENSDARG00000018404', 'ENSDARG00000041952'),]
housekeeping <- assay(ntd)[c('ENSDARG00000007682', 'ENSDARG00000101061', 'ENSDARG00000004754', 'ENSDARG00000044267', 'ENSDARG00000032175', 'ENSDARG00000044301', 'ENSDARG00000038835'),]
#assign gene names from biomart
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
GermDevGenesOrder <- GermDevGenes[order(GermDevGenes$Tdrd7_5mm_PGC_rep1_2018),]
pheatmap(GermDevGenesOrder, cluster_rows = FALSE, cluster_cols = FALSE)
## get genes upregulated by MO
setwd('../data/')
up_MO <- read.csv('UpregulatedByTdrd7MOinPGCs.txt', sep = '\t', row.names = 1)
down_MO <- read.csv('DownregulatedByTdrd7MOinPGCs.txt', sep = '\t', row.names = 1)

log_fold_cor[rownames(up_MO),]$genes <- 'upMO'
log_fold_cor[rownames(down_MO),]$genes <- 'downMO'

library(ggplot2)
p <- ggplot(log_fold_cor,aes(x=log_fold_cor$log2FoldChange.x, y=log_fold_cor$log2FoldChange.y, color = genes))
p + geom_point(alpha=0.7) + facet_wrap(~genes) + theme_bw() +
  geom_vline(xintercept=0, linetype="longdash", colour="black", size=0.4) +
  #yintercept is the -log10(pvalue): if pad < 0.05 -> -log10(0.05) = 1.3
  geom_hline(yintercept=0, linetype="longdash", colour="black", size=0.4) 

## levels of gene expression
norm_data <- counts(good_rep, normalized=TRUE)
norm_data_upMO <- as.data.frame(norm_data)[rownames(up_MO),]
norm_data_downMO <- as.data.frame(norm_data)[rownames(down_MO),]
norm_data_germ_cell <- as.data.frame(norm_data)[germ_cell_go,]

boxplot(norm_data_upMO$Tdrd7_5mm_PGC_rep1_2018, norm_data_upMO$Tdrd7_Rescue_PGC_rep1,
        norm_data_upMO$Tdrd7_MO_PGC_rep1_2018, outline=F, col=c('darkcyan', 'cyan3', 'darkgreen'),
        ylim=c(0,2000))
boxplot(norm_data_downMO$Tdrd7_5mm_PGC_rep1_2018, norm_data_downMO$Tdrd7_Rescue_PGC_rep1,
        norm_data_downMO$Tdrd7_MO_PGC_rep1_2018, outline=F, col=c('darkcyan', 'cyan3', 'darkgreen'),
        ylim=c(0,1200))
boxplot(norm_data_germ_cell$Tdrd7_5mm_PGC_rep1_2018, norm_data_germ_cell$Tdrd7_Rescue_PGC_rep1,
        norm_data_germ_cell$Tdrd7_MO_PGC_rep1_2018, outline=F, col=c('darkcyan', 'cyan3', 'darkgreen'),
        ylim=c(0,700))

## wilcoxon test
wilcox.test(norm_data_upMO$Tdrd7_5mm_PGC_rep1_2018, norm_data_upMO$Tdrd7_Rescue_PGC_rep1)
wilcox.test(norm_data_upMO$Tdrd7_5mm_PGC_rep1_2018, norm_data_upMO$Tdrd7_MO_PGC_rep1_2018)

wilcox.test(norm_data_downMO$Tdrd7_5mm_PGC_rep1_2018, norm_data_downMO$Tdrd7_Rescue_PGC_rep1)
wilcox.test(norm_data_downMO$Tdrd7_5mm_PGC_rep1_2018, norm_data_downMO$Tdrd7_MO_PGC_rep1_2018)

wilcox.test(norm_data_germ_cell$Tdrd7_5mm_PGC_rep1_2018, norm_data_germ_cell$Tdrd7_Rescue_PGC_rep1)
wilcox.test(norm_data_germ_cell$Tdrd7_5mm_PGC_rep1_2018, norm_data_germ_cell$Tdrd7_MO_PGC_rep1_2018)

### correlation with ATAC seq

setwd('~/Desktop/Postdoc/Data_Analysis/ATACseq/')
atac.fold.change.pgc <- readRDS('ATAC_FoldChange_PGC_5mm_rep2_resize20.rds')
atac.fold.change.soma <- readRDS('ATAC_FoldChange_Soma_MO_rep1_resize20.rds')

## define window
library(GenomicRanges)
library(biomaRt)
res_upRescue_ord <- head(res_MO_sub[order(res_MO_sub$log2FoldChange),], 200)
res_downRescue_ord <- tail(res_MO_sub[order(res_MO_sub$log2FoldChange),], 200)
## get fish promoters
mart <- useMart(host='http://oct2014.archive.ensembl.org',
                biomart='ENSEMBL_MART_ENSEMBL',
                dataset = "drerio_gene_ensembl")

fish_genes <- getBM(attributes = c("ensembl_gene_id",
                                    "ensembl_transcript_id",
                                    "chromosome_name",
                                    "transcript_start",
                                    "transcript_end",
                                    "strand",
                                    "gene_biotype"),
                     mart = mart)
## change strand names
fish_genes$strand <- ifelse(fish_genes$strand == 1, '+', '-')
fish_genes$chromosome_name <- paste0('chr', fish_genes$chromosome_name)
## subset
res_upRescue_ord_coord <- merge(data.frame(res_upRescue_ord), fish_genes, by.x = 0, by.y = 'ensembl_gene_id')
res_downRescue_ord_coord <- merge(data.frame(res_downRescue_ord), fish_genes, by.x= 0, by.y = 'ensembl_gene_id')

res_upRescue_ord_range <- GRanges(seqnames = res_upRescue_ord_coord$chromosome_name,
                                  IRanges(start = res_upRescue_ord_coord$transcript_start, end = res_upRescue_ord_coord$transcript_end),
                                  strand = res_upRescue_ord_coord$strand)

res_downRescue_ord_range <- GRanges(seqnames = res_downRescue_ord_coord$chromosome_name,
                                    IRanges(start = res_downRescue_ord_coord$transcript_start, end = res_downRescue_ord_coord$transcript_end),
                                    strand = res_downRescue_ord_coord$strand)

PromotersUpRescue <- promoters(res_upRescue_ord_range, upstream = 500,
                                        downstream = 500)
PromotersDownRescue <- promoters(res_downRescue_ord_range, upstream = 500,
                                downstream = 500)

## chromatin accessibility of subset of genes
library(genomation)
sm.pgc.up <- ScoreMatrix(target = atac.fold.change.pgc, windows = PromotersUpRescue, weight.col = 'score', strand.aware = T)
sm.pgc.down <- ScoreMatrix(target = atac.fold.change.pgc, windows = PromotersDownRescue, weight.col = 'score', strand.aware = T)

sm.soma.up <- ScoreMatrix(target = atac.fold.change.soma, windows = PromotersUpRescue, weight.col = 'score', strand.aware = T)
sm.soma.down <- ScoreMatrix(target = atac.fold.change.soma, windows = PromotersDownRescue, weight.col = 'score', strand.aware = T)

## metaplot
library(dplyr)
### ggplot meta
df <- data.frame(cbind(colMeans(sm.pgc.up), colMeans(sm.pgc.down)))
colnames(df) <- c("PGCup", "PGCdown")

df.gg <- df %>% gather(Expression, Tn5CutSites, PGCup:PGCdown)
df.gg$Index <- rep(c(-500:-1, 1:500), 2)

ggplot(df.gg, aes(x=Index, y=Tn5CutSites, Group=factor(Expression))) +
  geom_line(aes(colour=factor(Expression)), size = 1.5) +
  scale_colour_manual(values = c("seagreen4", "slateblue")) + ###7CAE00 ###F8766D +
  geom_vline(xintercept=0, colour="grey", linetype="longdash") +
  theme_classic() + theme(legend.position='none')


#######
## enhancers
##requires: resMOsub, pgc_tag_cluster
library(rtracklayer)
setwd('../Data_Analysis/ATACseq/')

atac.fold.change.pgc <- readRDS('ATAC_FoldChange_Coverage120ATAC.PGC_24hpf_Tdrd7_MO5mismatch.rep2.rds')
#readRDS('ATAC_FoldChange_PGC_5mm_rep2_resize20.rds')
atac.fold.change.soma <- readRDS('ATAC_FoldChange_Coverage120ATAC.somatic_24hpf_Tdrd7_MO.rep1.rds')
## summit to center the peaks
setwd('~/Desktop/Postdoc/Data_Analysis/ATACseq/Summit_MACS/')
file_summit <- list.files()[grep('summit', list.files())]

.import.summit <- function(f_summit){
  summit <- read.csv(f_summit, sep = '\t', header = F)[,c(1:3)]
  colnames(summit) <- c('seqnames', 'start', 'end')
  s_range <- GRanges(seqnames = summit$seqnames, IRanges(start = summit$start, end = summit$start))
  return(s_range)
}

pgc.summit <- .import.summit(file_summit[1])
soma.summit <- .import.summit (file_summit[2])
#import k27ac data
k27ac <- read.csv('../ATAC analysis/enhancers/ChIP-seq_Skarmeta_Lab_H3K27ac_0002AS.DCD000652SQ.USERpanosfirbas.R1.filt.deduplicated.nodup_pooled.naive_overlap.narrowPeak.gz', sep = '\t', header = F)
colnames(k27ac) <- c('chromosome', 'start', 'end', 'name', 'V5', 'strand', 'fold', 'pvalue', 'qvalue', 'V10')
k27ac <- toGRanges(k27ac)

## only overlapping k27
pgc.summit <- subsetByOverlaps(pgc.summit, k27ac)
soma.summit <- subsetByOverlaps(soma.summit, k27ac)

## add degene class
res_MO_sub$degenes <- 'NS'
res_MO_sub$degenes[res_MO_sub$log2FoldChange < 0 & res_MO_sub$padj < 0.05] <- 'Tdrd7_up'
res_MO_sub$degenes[res_MO_sub$log2FoldChange > 0 & res_MO_sub$padj < 0.05] <- 'Tdrd7_down'

## load cage set
myCAGEset <- readRDS("~/Desktop/Postdoc/Data_Analysis/CAGEsets/CAGEset_PGC_soma_Early_Late.rds")
## get tag cluster for PGC prim5
sample.PGC <- unname(sampleLabels(myCAGEset)[3])
pgc_tag_cluster <- tagClusters(myCAGEset, sample = sample.PGC, returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
## annotate peaks
source('~/Desktop/Postdoc/R_scripts/gcCAGE/R/Functions.R')
txdb <- loadDb('~/Desktop/Postdoc/Data_Analysis/annotation/txdb_DanRer7.sqlite')
pgc_tag_cluster_anno <- annotate.cage.peaks(pgc_tag_cluster)

pgc_tag_cluster_anno_merged_Tdrd7 <- merge(pgc_tag_cluster_anno, data.frame(res_MO_sub), by.x = 'geneId', by.y = 0)
pgc_tag_cluster_anno_merged_Tdrd7$geneId <- NULL

## get the tss
## extend 5kb up and downstream
.centerTSS2_Tdrd7 <- function(x, width = 1000, class = class){
  new.data.frame <- data.frame(x)
  new.data.frame <- subset(new.data.frame, new.data.frame$degene == class)
  new.data.frame <- na.omit(new.data.frame)
  
  if(class == 'Tdrd7_up'){
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
  TSS <- resize(zebrafishPromotersTSS, width = width, fix = 'center')
  
  return(TSS)
}
TSS_5kb_up_Tdrd7 <- .centerTSS2_Tdrd7(pgc_tag_cluster_anno_merged_Tdrd7, width = 5000, class = 'Tdrd7_up')
TSS_500_up_Tdrd7 <- .centerTSS2_Tdrd7(pgc_tag_cluster_anno_merged_Tdrd7, width = 1000, class = 'Tdrd7_up')

TSS_5kb_down_Tdrd7 <- .centerTSS2_Tdrd7(pgc_tag_cluster_anno_merged_Tdrd7, width = 5000, class = 'Tdrd7_down')
TSS_500_down_Tdrd7<- .centerTSS2_Tdrd7(pgc_tag_cluster_anno_merged_Tdrd7, width = 1000, class = 'Tdrd7_down')

## get atac peaks  and make grange with these as start
.clean.peaks <- function(tss5k, tss500){
  centered_atac_peaks <- subsetByOverlaps(pgc.summit, tss5k)
  ## remove promoters
  non_promoters_atac_peaks <- subsetByOverlaps(centered_atac_peaks, tss500, invert = TRUE)
  
  #non_promoters_atac_peaks <- subset(non_promoters_atac_peaks, non_promoters_atac_peaks$score > 0)
  #non_promoters_atac_peaks <- reduce(non_promoters_atac_peaks, min.gapwidth = 100)
  non_promoters_atac_peaks <- resize(non_promoters_atac_peaks, width = 1500, fix = 'center')
}
cleaned_peak_up_Tdr7 <- .clean.peaks(TSS_5kb_up_Tdrd7, TSS_500_up_Tdrd7)
cleaned_peak_down_Tdrd7 <- .clean.peaks(TSS_5kb_down_Tdrd7, TSS_500_down_Tdrd7)

## cluster peaks

## plot heatmaps
library(genomation)
library(tidyr)
m_Tdrd7_up <- ScoreMatrix(target = atac.fold.change.pgc, windows = cleaned_peak_up_Tdr7, weight.col = 'score', strand.aware = T)
m_Tdrd7_down <- ScoreMatrix(target = atac.fold.change.pgc, windows = cleaned_peak_down_Tdrd7, weight.col = 'score', strand.aware = T)

heatMatrix(m_Tdrd7_up, winsorize = c(0,95), clustfun = function(x) kmeans(x, centers=2)$cluster, 
           col = brewer.pal(n = 5, "Blues"))
## plot meta
df <- data.frame(cbind(colMeans(m_Tdrd7_up), colMeans(m_Tdrd7_down)))
colnames(df) <- c("Up_MO", "Up_Rescue")

df.gg <- df %>% gather(Expression, Tn5CutSites, Up_MO:Up_Rescue)
df.gg$Index <- rep(c(-750:-1, 1:750), 2)

ggplot(df.gg, aes(x=Index, y=Tn5CutSites, Group=factor(Expression))) +
  geom_line(aes(colour=factor(Expression)), size = 1.5) +
  scale_colour_manual(values = c("slateblue", "seagreen4")) + ###7CAE00 ###F8766D +
  geom_vline(xintercept=0, colour="grey", linetype="longdash") +
  theme_classic() + theme(legend.position='bottom')

### differential openess 
atacUp <- read.csv('../ATAC analysis/TablesATAC/All_ATAC_Upreg_PGC_Prim5.txt', sep = '\t')
atacDown <- read.csv('../ATAC analysis/TablesATAC/All_ATAC_Upreg_Soma_Prim5.txt', sep = '\t')

convert.to.grange <- function(x){
  x$chr <- as.character(lapply(strsplit(rownames(x), "\\:"), function(u) u[[1]][1]))
  x$range <- as.character(lapply(strsplit(rownames(x), "\\:"), function(u) u[[2]][1]))
  x$start <- as.character(lapply(strsplit(x$range, "\\-"), function(u) u[[1]][1]))
  x$end <- as.character(lapply(strsplit(x$range, "\\-"), function(u) u[[2]][1]))
  x$range <- NULL
  return(x)
}
atacUpRange <- convert.to.grange(atacUp)
atacDownRange <- convert.to.grange(atacDown)

library(AnnotationDbi)
txdb <- loadDb('../annotation/txdb_DanRer7.sqlite')
ensembl.ids.atac <- function(x){
  grange.obj <- toGRanges(x)
  anno.ensembl.id <- annotatePeakInBatch(grange.obj, AnnotationData=toGRanges(txdb), 
                                         output="nearestBiDirectionalPromoters",
                                         bindingRegion=c(-10000, 10000))
  anno.ensembl.id.order <- anno.ensembl.id[order(abs(anno.ensembl.id$log2FoldChange)),]
  anno.ensembl.id.200 <- tail(anno.ensembl.id.order, 50)
  ensembl.ids <- anno.ensembl.id.200$feature
  return(ensembl.ids)
  
}
## ENSEMBL IDs of genes with Upreg ATAC peaks
ids.PGC <- ensembl.ids.atac(atacUpRange)
ids.Soma <- ensembl.ids.atac(atacDownRange)

## plot log
res_MO_res_Up <- as.data.frame(res_MO_res)[ids.PGC,]
res_MO_res_Down <- as.data.frame(res_MO_res)[ids.Soma,]
res_MO_res_Up[is.na(res_MO_res_Up)] <- 0
res_MO_res_Down[is.na(res_MO_res_Down)] <- 0

boxplot(res_MO_res_Up$log2FoldChange, res_MO_res_Down$log2FoldChange, outline = F)
## plot normalized reads

mat_mm <- counts(rescue_mm_dds, normalized=TRUE) #normalised read counts (raw/size factor)
mat_MO <- counts(rescue_MO_dds, normalized=TRUE) #normalised read counts (raw/size factor)

mat_MO_Up <- na.omit(as.data.frame(mat_MO)[ids.PGC,])
mat_MO_Down <- na.omit(as.data.frame(mat_MO)[ids.Soma,])

boxplot(log(mat_mm_Up$Tdrd7_5mm_PGC_rep1_2018+1), log(mat_mm_Down$Tdrd7_Rescue_PGC_rep1+1))

ntd <- normTransform(dds_new)



