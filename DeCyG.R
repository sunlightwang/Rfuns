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
  gene.names <- rownames(data)
  cv.vec <- foreach(i = 1:(nrow(data)-1), .combine = "c", 
                    .export=ls(envir=globalenv()), .packages='matrixStats') %dopar% {
    cv <- rep(NA, nrow(data)-i)
    cmp <- rep(NA, nrow(data)-i)
    n <- 0
    for(j in (i+1):nrow(data)) {
      n <- n+1
      cycGD.res <- cycG_dect(rbind(data[i,], data[j,]))
      cmp[n] <- paste0(gene.names[i], ".vs.", gene.names[j])
      pole1 <- cycGD.res[,1]
      cv[n] <- sqrt(var(pole1)) / mean(pole1)
    }
    names(cv) <- cmp
    cv
  }
  cmp.sel <- names(head(sort(cv.vec), n=topN))
  g1 <- sapply(cmp.sel, function(x) unlist(strsplit(x,'[.]'))[1])
  g2 <- sapply(cmp.sel, function(x) unlist(strsplit(x,'[.]'))[3])
  for(i in 1:topN) { plot(data[g1[i],], data[g2[i],], pch=20, xlab=g1[i], ylab=g2[i]) }
  cv.vec 
}
               
cycG_perm <- function(data, nonexpr.filter=F, nonexpr.q=0.1, 
                      p=8, seed=9999) { # data normalized, row - genes, col - samples
  require(doParallel)
  cl <- makeCluster(p)
  registerDoParallel(cl)
  gene.names <- rownames(data)
  results <- foreach(i = 1:(nrow(data)-1), .combine = "rbind", 
                    .export=ls(envir=globalenv()), .packages='matrixStats') %dopar% {
    temp <- matrix(NA, nrow(data)-i, 12)
    cmp <- rep(NA, nrow(data)-i)
    n <- 0
    for(j in (i+1):nrow(data)) {
      n <- n+1
      cmp[n] <- paste0(gene.names[i], ".vs.", gene.names[j])
      if(nonexpr.filter) {
        idx <- rep(TRUE, ncol(data))
        if(quantile(data[i,], 0) == quantile(data[i,], nonexpr.q)) idx <- idx & data[i,] > quantile(data[i,], 0)
        if(quantile(data[j,], 0) == quantile(data[j,], nonexpr.q)) idx <- idx & data[j,] > quantile(data[j,], 0)
        cycGD.res <- cycG_dect(rbind(data[i,idx], data[j,idx]))
        set.seed(seed)
        perm.cycGD.res <- cycG_dect(rbind(data[i,idx], sample(data[j,idx])))
      } else {
        cycGD.res <- cycG_dect(rbind(data[i,], data[j,]))
        set.seed(seed)
        perm.cycGD.res <- cycG_dect(rbind(data[i,], sample(data[j,])))
      }
      pole1 <- cycGD.res[,1] / mean(cycGD.res[,1])
      perm.pole1 <- perm.cycGD.res[,1] / mean(perm.cycGD.res[,1])
      temp[n, 1] <- wilcox.test(pole1, perm.pole1, alternative="greater")$p.val #wilcox.p
      temp[n, 2]  <- ks.test(pole1, perm.pole1, alternative="less")$p.val       #ks.p
      temp[n, 3:7] <- quantile(pole1, c(0.1, 0.25, 0.5, 0.75, 0.9))
      temp[n, 8:12] <- quantile(perm.pole1, c(0.1, 0.25, 0.5, 0.75, 0.9))
    }
    rownames(temp) <- cmp
    temp 
  }
  colnames(results) <- c("wilcox.p", "ks.p", "Q10", "Q25", "Q50", "Q75", "Q90", "P_Q10", "P_Q25", "P_Q50", "P_Q75", "P_Q90")
  results 
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
