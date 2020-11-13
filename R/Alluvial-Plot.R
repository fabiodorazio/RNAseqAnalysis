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
         freq=fly$n, border=NA, alpha = 0.5,
         col=case_when(alluv_tb$sPrim5 == "Hyper" ~ "red",
                       alluv_tb$sPrim5 == "Hypo" ~ "blue",
                       TRUE ~ "grey"),
         cex=0.75,
         axis_labels = c("dome", "10somites", 'prim5'))


