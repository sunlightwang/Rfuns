library(matrixStats)
#TODO: filter out bimodel genes first
center.scale <- function(data)  { ###TODO: center scale maybe not good for asymmetric situation; alt. 5-95% interval rescale: check PCA..
  (data - rowMeans(data)) / sqrt(rowVars(data)) }

cart2pol <- function(data) {
  do.call(rbind, sapply(1:nrow(data), function(x) {
    r <- sqrt(data[x,1]^2 + data[x,2]^2)
    theta <- atan(data[x,2] / data[x,1])
    if(data[x,2]<0 & data[x,1]<0) {theta=theta-pi}
    if(data[x,2]>0 & data[x,1]<0) {theta=theta+pi}
    return(cbind(r, theta))
  }, simplify=F))
}
  
pol2cart <- function(data) {
  do.call(rbind, sapply(1:nrow(data), function(x) {
    xx <- data[x,1] * cos(data[x,2])
    yy <- data[x,1] * sin(data[x,2])
    return(cbind(xx, yy))
  }, simplify=F))
}

cycG_dect <- function(GP) {
  GP.scaled <- center.scale( GP )
  #plot(t(GP.scaled))
  pca <- prcomp(t(GP.scaled)) 
  # pca$x = t(GP.scaled) %*% pca$rotation
  pca.rotated <- t(GP.scaled) %*% pca$rotation
  pca.rotated.scaled <- t( t(pca.rotated) / pca$sdev )
  #plot(pca.rotated.scaled)
  #abline(v=0,h=0,lty=2,col="grey80")
  GP.pole <- cart2pol(pca.rotated.scaled)
  #plot(GP.pole)
  #hist(GP.pole[,1])
  GP.pole
}

cycG_dect_wrapper <- function(data, topN=100) { # data normalized, row - genes, col - samples
  cmp.no <- choose(nrow(data), 2)
  cmp <- rep(NA, cmp.no)
  pole1 <- matrix(NA, nrow=cmp.no, ncol=ncol(data))
  pole2 <- matrix(NA, nrow=cmp.no, ncol=ncol(data))
  gene.names <- rownames(data)
  n <- 0
  
  for(i in 1:(nrow(data)-1)) { 
    for(j in (i+1):nrow(data)) {
      n <- n + 1
      #print(paste(i,j))
      cycGD.res <- cycG_dect(rbind(data[i,], data[j,]))
      cmp[n] <- paste0(gene.names[i], ".vs.", gene.names[j])
      pole1[n,] <- cycGD.res[,1]
      pole2[n,] <- cycGD.res[,2]
    }
  }
  rownames(pole1) <- cmp
  colnames(pole1) <- colnames(data)
  rownames(pole2) <- cmp
  colnames(pole2) <- colnames(data)
  #write.table(pole1, "H1.pole1.txt", sep="\t", quote=F)
  #write.table(pole2, "H1.pole2.txt", sep="\t", quote=F)
  cv <- sqrt(apply(pole1, 1, var)) / apply(pole1, 1, mean)
  cmp.sel <- names(head(sort(cv), n=topN))
  g1 <- sapply(cmp.sel, function(x) unlist(strsplit(x,'[.]'))[1])
  g2 <- sapply(cmp.sel, function(x) unlist(strsplit(x,'[.]'))[3])
  for(i in 1:topN) { plot(data[g1[i],], data[g2[i],], pch=20, xlab=g1[i], ylab=g2[i]) }
}

cycG_dect_wrapper.p <- function(data, topN=100, p=8) { # data normalized, row - genes, col - samples
  require(doParallel)
  cl <- makeCluster(p)
  registerDoParallel(cl)
  cmp.no <- choose(nrow(data), 2)
  gene.names <- rownames(data)
  cv.vec <- foreach(i = 1:(nrow(data)-1), .combine = "c", 
                    .export=ls(envir=globalenv()), .packages='matrixStats') %dopar% {
    cv <- rep(NA, nrow(data)-i)
    n <- 0
    for(j in (i+1):nrow(data)) {
      n <- n+1
      cycGD.res <- cycG_dect(rbind(data[i,], data[j,]))
      cmp <- paste0(gene.names[i], ".vs.", gene.names[j])
      pole1 <- cycGD.res[,1]
      cv[n] <- sqrt(var(pole1)) / mean(pole1)
    }
    cv
  }
  cmp.sel <- names(head(sort(cv.vec), n=topN))
  g1 <- sapply(cmp.sel, function(x) unlist(strsplit(x,'[.]'))[1])
  g2 <- sapply(cmp.sel, function(x) unlist(strsplit(x,'[.]'))[3])
  for(i in 1:topN) { plot(data[g1[i],], data[g2[i],], pch=20, xlab=g1[i], ylab=g2[i]) }
  cv.vec 
}
               
cycG_simu <- function(seed=1000) {
  set.seed(seed)
  G1 <- rnorm(1000)
  G2 <- rnorm(1000)*20
  GP <- rbind(G1, G2)
  plot(t(GP))
  GP
}

cycG_simu2 <- function() {
  ####
  r <- rnorm(1000,1,0.2)
  th <- runif(1000, -pi, pi)
  circ <- t(pol2cart(cbind(r,th)))
  plot(t(circ))
  circ.scaled <- center.scale( circ )
  plot(t(circ.scaled))
  pca <- prcomp(t(circ.scaled))
  pca.rotated <- t(circ.scaled) %*% pca$rotation
  pca.rotated.scaled <- t( t(pca.rotated) / pca$sdev )
  plot(pca.rotated.scaled)
  abline(v=0,h=0,lty=2,col="grey80")
  GP.pole <- cart2pol(pca.rotated.scaled)
  plot(GP.pole)
  hist(GP.pole[,1])
}
