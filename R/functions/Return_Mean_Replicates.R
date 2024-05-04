## Mean for replicates
# create an empty data frame
return_mean_replicates <- function(x){
  tpm1 <- data.frame(row.names = rownames(x))
  for(n in 1:ncol(x)){
    tryCatch({
      ## n1 is the odd columns
      n1 <- n[n %% 2 != 0]
      tpm1[,n] <- rowMeans(x[,c(n1,n1+1)])
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  colnames(tpm1) <- colnames(x)
  tpm1 <- t(na.omit(t(tpm1)))
  return(data.frame(tpm1))
}
