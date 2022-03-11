library(dplyr)
library(purrr)

path.to.file <- '../Intron_Retention/Result_Tables/'
files <- list.files()

load.IR.table <- function(x){
  out <- read.csv(paste0(path.to.file, x), sep = '\t', header = T)
  names <- x
  out$sample <- names
  out
}

IR.tables <- lapply(files, load.IR.table)
## filter out low coverage, extracting ENSEMBL IDs and sum all values with the same identifier
subset.by.coverage <- function(x){
  out <- x #subset(x, x$Warnings == '-') ## if filter out low coverage is necessary only
  ## extract string between slashes
  out$ENSEMBL <- strsplit(as.character(out$Name), "/", TRUE) %>% map(`[`(2)) %>% unlist()
  # merge rows with same identifier and sum them
  out <- tapply(out$IRratio, out$ENSEMBL, sum)
  out <- data.frame(out)
  out <- data.frame(IRratio = out$out, ENSEMBL = rownames(out), row.names = NULL)
  return(out)
}
IR.tables.sub <- lapply(IR.tables, subset.by.coverage)

## split list in single data frames
for (i in seq(IR.tables.sub))
  assign(paste0("IR.table", i), IR.tables.sub[[i]])

## merges replicates
test.function <- function(table1, table2){
  frame <- rbind(table1,table2)
  frame <- tapply(frame$IRratio, frame$ENSEMBL, mean)
  frame <- data.frame(frame)
  return(frame)
}


## subset by early genes
path.to.early.genes <- "../New Analysis MBT/"
filter.early.genes <- function(x){
  early.pgc <- read.csv(paste0(path.to.early.genes,"UpregInPGCOnlyAtHigh.txt"), sep = '\t')
  early.pgc.names <- merge(x, early.pgc, by = 0)
  #table <- subset(x, x[early.pgc.names,]
  #table <- na.omit(table)
  return(early.pgc.names)
}



## random sample and permutation
permute.early.genes <- function(x, iter = 10000, size = 152){
  distr <- integer(iter)
  for(i in 1:iter) {
    distr[i] <- mean(x[sample(1:nrow(x), size=size, replace=F),])
  }
  #table <- subset(x, x[early.pgc.names,]
  #table <- na.omit(table)
  return(distr)
}



# how many out of iter have lower score than observed -> empirical P
p_Emp <- length(which(distr > mean(merged.PGC.High$frame))) / iter

