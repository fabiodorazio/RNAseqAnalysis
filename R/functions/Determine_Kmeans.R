-------------------
  #Determining number of k
  -------------------
  
  #Within group Sum of Squares (WSS) Criterion : k6-~k15
  clusterdata.n <- totalReadsCount
  wss <- per <- rep(0,20)
  wss[1] <- (nrow(clusterdata.n)-1) * sum(apply(clusterdata.n, 2, var)) #== kmeans(clusterdata.n,1)$withinss
  per[1] <- 1 - wss[1]/wss[1]
  for (i in 2:20) {
    wss[i] <- sum(kmeans(clusterdata.n, centers = i, iter.max = 500, nstart=50, algorithm="Lloyd")$withinss)
    per[i] <- round(1 - wss[i]/wss[1], 3)*100
  }
  pdf('./plots/K_estimation_WSS_criterion.pdf')
  plot(1:20, wss, type='b', col='cadetblue', xlab='Number of Clusters', ylab='Within group sum of squares', 
       main='Variance within Clusters\n(Normalized Data)', lwd=2.0, pch=16)
  abline(v=c(5,15), col='firebrick2', lwd=2)
  dev.off()
  pdf('./plots/K_estimation_Variance_WSS_criterion.pdf')
  plot(1:25, per, type='b', col='cadetblue', xlab='Number of Clusters', ylab='Percentage of variance explained', 
       main='Percentage of variance explained', lwd=2.0, pch=16)
  abline(v=c(5,15), col='firebrick2', lwd=2)
  text(c(5,15),68, paste0(as.character(per[c(5,15)]), '%'), col= 'firebrick2', cex=1.5)
  dev.off()
