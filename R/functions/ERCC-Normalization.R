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

totalReadsCount <- sapply(allReadsCounts, function(x) x$totalCount)
rownames(totalReadsCount) <- row.names(allReadsCounts[[1]])
totalReadsCount <- totalReadsCount[,c(1:20)]
colnames(totalReadsCount) <- c('s256PGC1', 's256Soma1', 's256PGC2', 's256Soma2', 'sHighPGC1', 'sHighSoma1', 'sHighPGC2', 'sHighSoma2',  'sDomePGC1', 'sDomeSoma1', 'sDomePGC2', 'sDomeSoma2', 's10somitesPGC1', 's10somitesSoma1', 's10somitesPGC2', 's10somitesSoma2', 'sPrim5PGC1', 'sPrim5Soma1', 'sPrim5PGC2', 'sPrim5Soma2')
head(totalReadsCount)
totalReadsCount <- tail(totalReadsCount[,c(1,2,5,6,7,8,17,18,19,20)], 92)

#import gtf file for gene lenghts
setwd("~/Desktop/PhD-March-2019-backup/FabioRNAseq (2)/2108/RNAseqtables/ERCC/")
gtf.with.lenghts <- read.csv('ERCC92/ERCC92.gtf', sep = '\t', header = F)
gtf.with.lenghts <- data.frame(row.names = gtf.with.lenghts$V1, 'lenght' = gtf.with.lenghts$V5)

norm.lenght.function <- function(x){
   x/(gtf.with.lenghts$lenght/1e3)
}
norm.lenght <- apply(totalReadsCount, 2, norm.lenght.function)

#get normalization factors

scaling.factors <- apply(totalReadsCount, 2, function(y) (sum(y))/1e6)
# rpkm
rpkm <- norm.lenght/scaling.factors[col(norm.lenght)]

rpkm.256 <- data.frame(row.names = rownames(rpkm), 'PGC' = rpkm[,1], 'Soma' = rpkm[,2])
rpkm.high <- data.frame(row.names = rownames(rpkm), 'PGC' = rowMeans(as.data.frame(rpkm)[c('sHighPGC1', 'sHighPGC2')],), 'Soma' = rowMeans(as.data.frame(rpkm)[c('sHighSoma1', 'sHighSoma2')],))
rpkm.prim <- data.frame(row.names = rownames(rpkm), 'PGC' = rowMeans(as.data.frame(rpkm)[c('sPrim5PGC1', 'sPrim5PGC2')],), 'Soma' = rowMeans(as.data.frame(rpkm)[c('sPrim5Soma1', 'sPrim5Soma2')],))

rpkm.256.high.PGC.Soma <- data.frame(row.names = rownames(rpkm), 'PGC' = rpkm[,1], 'Soma' = rowMeans(as.data.frame(rpkm)[c('sHighSoma1', 'sHighSoma2')],))
rpkm.256.high.Soma.PGC <- data.frame(row.names = rownames(rpkm), 'PGC' = rowMeans(as.data.frame(rpkm)[c('sHighPGC1', 'sHighPGC2')],), 'Soma' = rpkm[,2])
rpkm.highPGC.primSoma <- data.frame(row.names = rownames(rpkm), 'PGC' = rowMeans(as.data.frame(rpkm)[c('sHighPGC1', 'sHighPGC2')],), 'Soma' = rowMeans(as.data.frame(rpkm)[c('sPrim5Soma1', 'sPrim5Soma2')],))

#set treshold of 1
rpkm.list <- list(rpkm.256, rpkm.high, rpkm.prim, rpkm.256.high.PGC.Soma, rpkm.256.high.Soma.PGC, rpkm.highPGC.primSoma)
df <- function(x){
  treshold <- subset(x, x$PGC & x$Soma > 1)
}
treshold.all <- lapply(rpkm.list, df)

# combine PGC and Soma

Obs.Ratio.256 <- data.frame(row.names = rownames(treshold.all[[1]]), treshold.all[[1]]$PGC/treshold.all[[1]]$Soma)
Obs.Ratio.high <- data.frame(row.names = rownames(treshold.all[[2]]), treshold.all[[2]]$PGC/treshold.all[[2]]$Soma)
Obs.Ratio.prim <- data.frame(row.names = rownames(treshold.all[[3]]), treshold.all[[3]]$PGC/treshold.all[[3]]$Soma)
Obs.Ratio.PGC256.somaHigh <- data.frame(row.names = rownames(treshold.all[[4]]), treshold.all[[4]]$PGC/treshold.all[[4]]$Soma)
Obs.Ratio.PGCHigh.soma256 <- data.frame(row.names = rownames(treshold.all[[5]]), treshold.all[[5]]$PGC/treshold.all[[5]]$Soma)
Obs.Ratio.PGCHigh.somaprim <- data.frame(row.names = rownames(treshold.all[[6]]), treshold.all[[6]]$PGC/treshold.all[[6]]$Soma)

Obs.Ratio.all <- list(Obs.Ratio.256, Obs.Ratio.high, Obs.Ratio.prim, Obs.Ratio.PGC256.somaHigh, Obs.Ratio.PGCHigh.soma256, Obs.Ratio.PGCHigh.somaprim)

#import control file and merge obs vs expected ratio
control <- read.csv('ERCC_Control_Analysis.txt', sep = '\t', header = T)

#assess concentration
#multiply concentration for the dilution factor
#1 ul of 1:10000 dilution in 5 ng = 0.0001
control$NormConcMix1 <- control$concentration.in.Mix.1..attomoles.ul.*0.0001
control$NormConcMix2 <- control$concentration.in.Mix.2..attomoles.ul.*0.0001

Mix1.Mix2 <- function(x){
  Obs.vs.Exp.Ratio <- merge(x, control, by.x = 0, by.y = 'ERCC.ID')
  Obs.vs.Exp.Ratio
}
Mix1.Mix2.control <- lapply(Obs.Ratio.all, Mix1.Mix2)

#import generate data tables
Exp.vs.obs.FC.256 <- data.frame(row.names = rownames(Mix1.Mix2.control[[1]]), 'Observed_FC' = Mix1.Mix2.control[[1]]$treshold.all..1...PGC.treshold.all..1...Soma, 'Expected_FC' = Mix1.Mix2.control[[1]]$expected.fold.change.ratio)
Exp.vs.obs.FC.high <- data.frame(row.names = rownames(Mix1.Mix2.control[[2]]), 'Observed_FC' = Mix1.Mix2.control[[2]]$treshold.all..2...PGC.treshold.all..2...Soma, 'Expected_FC' = Mix1.Mix2.control[[2]]$expected.fold.change.ratio)
Exp.vs.obs.FC.prim <- data.frame(row.names = rownames(Mix1.Mix2.control[[3]]), 'Observed_FC' = Mix1.Mix2.control[[3]]$treshold.all..3...PGC.treshold.all..3...Soma, 'Expected_FC' = Mix1.Mix2.control[[3]]$expected.fold.change.ratio)
Exp.vs.obs.FC.PGC256.somaHigh <- data.frame(row.names = rownames(Mix1.Mix2.control[[4]]), 'Observed_FC' = Mix1.Mix2.control[[4]]$treshold.all..4...PGC.treshold.all..4...Soma, 'Expected_FC' = Mix1.Mix2.control[[4]]$expected.fold.change.ratio)
Exp.vs.obs.FC.PGChigh.soma256 <- data.frame(row.names = rownames(Mix1.Mix2.control[[5]]), 'Observed_FC' = Mix1.Mix2.control[[5]]$treshold.all..5...PGC.treshold.all..5...Soma, 'Expected_FC' = Mix1.Mix2.control[[5]]$expected.fold.change.ratio)
Exp.vs.obs.FC.PGChigh.somaprim <- data.frame(row.names = rownames(Mix1.Mix2.control[[6]]), 'Observed_FC' = Mix1.Mix2.control[[6]]$treshold.all..6...PGC.treshold.all..6...Soma, 'Expected_FC' = Mix1.Mix2.control[[6]]$expected.fold.change.ratio)


list.exp.obs <- list(Exp.vs.obs.FC.256, Exp.vs.obs.FC.high, Exp.vs.obs.FC.prim, Exp.vs.obs.FC.PGC256.somaHigh, Exp.vs.obs.FC.PGChigh.soma256, Exp.vs.obs.FC.PGChigh.somaprim)
list.exp.obs.titles <- c('Exp.vs.obs.FC.256', 'Exp.vs.obs.FC.high', 'Exp.vs.obs.FC.prim')

#plot and correlation
test <- function(x){
  plot(x$Observed_FC, x$Expected_FC, pch = 16)
  abline(lm(x$Expected_FC ~ x$Observed_FC))
  
}
plot.exp.vs.obs <- lapply(list.exp.obs, test)

#print correlations
for(i in list.exp.obs){
  print(cor(i$Observed_FC, i$Expected_FC))
}

## generate calibration curve for concentration and rpkm
#generate data.frames of dose.response
dose.response <- function(rpkm, control = control, y, z){
  control.dose <- merge(rpkm, control, by.x = 0, by.y = 'ERCC.ID')
  control.dose <- subset(control.dose, control.dose[,y] > 0)
  
  plot(log(control.dose[,y]), log(control.dose[,z]), pch = 16, col = 'blue')
  line.dose <- lm(log(control.dose[,z]) ~ log(control.dose[,y]))
  abline(line.dose)
  cor(log(control.dose[,y]), log(control.dose[,z]))
  return(line.dose)
}
line256PGC1 <- dose.response(rpkm.256, control, 'PGC', 'concentration.in.Mix.1..attomoles.ul.')
line256Soma1 <- dose.response(rpkm.256, control, 'Soma', 'concentration.in.Mix.2..attomoles.ul.')
lineHighPGC1 <- dose.response(rpkm.high, control, 'PGC', 'concentration.in.Mix.1..attomoles.ul.')
lineHighSoma1 <- dose.response(rpkm.high, control, 'Soma', 'concentration.in.Mix.2..attomoles.ul.')


dose.response.PGC <- function(x){
  dose.resp.PGC <- merge(x, control, by.x = 0, by.y = 'ERCC.ID')
}

list.dose.response.PGC <- lapply(rpkm.list, dose.response.PGC)



## normalise rpkm based on the equation
## import rpkm genome-wide
setwd('../MBT/')
transcript_level_expr <- read.csv('rpkm_Allstages.txt', sep = '\t')

setwd('../New Analysis MBT/')
up.pgc <- read.csv('UpregInPGCOnlyAtHigh.txt', sep = '\t')
germ_genes <- rownames(up.pgc)
GermDevGenes <- transcript_level_expr[germ_genes,]
GermDevGenes <- na.omit(GermDevGenes)

#equation of a line for concentration
#divide all the concentration values by 0.0005 (dilution factor used in the experiment)
#this gives concentration in attomoles/1ng of total RNA
equation.of.line <- function(x,y){
  intercept <- summary(x)$coefficients[1]
  slope <- summary(x)$coefficients[2]

  norm.conc.256.soma <- GermDevGenes[,y]
  answer <- (norm.conc.256.soma - intercept)/slope
  return(answer)
  
}
conc.256.pgc.germ <- equation.of.line(line256PGC1, 's256PGC1')
conc.256.soma.germ <- equation.of.line(line256Soma1, 's256Soma1')
conc.high.pgc.germ <- equation.of.line(lineHighPGC1, 'sHighPGC1')
conc.high.soma.germ <- equation.of.line(lineHighSoma1, 'sHighSoma1')



#generates lists of values that can be plotted
boxplot(conc.256.pgc.germ, conc.256.soma.germ, conc.high.pgc.germ, conc.high.soma.germ,
        outline = F, col = c('pale green', 'purple', 'spring green', 'purple'), ylim = c(0,500))


x <- data.frame( PGC = answer, Soma = answer.soma)
write.table(x, 'ERCC-norm-rpkm-256-PGC-Soma.txt', sep = '\t')

