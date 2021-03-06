library(matrixStats)
#NOTE: bimodel genes detection 
#TODO: geneset: with correlation filter out genes
center.scale <- function(data, method=c("quantile","mean"), quantile_low=0.05, quantile_high=0.95)  { 
  if(method %in% "mean")  (data - rowMeans(data)) / sqrt(rowVars(data)) 
  if(method %in% "quantile") {
    d1 <- apply(data, 1, function(x) {
      lb <- quantile(x, quantile_low)
      ub <- quantile(x, quantile_high)
      x[x<lb] <- lb
      x[x>ub] <- ub
      (2 * x - lb - ub) / (ub - lb) 
    })
    t(d1)
  }
}

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
  pca.rotated.scaled <- t(center.scale(t(pca.rotated)))
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
  ret.value <- c("wilcox.p", "ks.l.p", "ks.g.p", "Rsq", "Q10", "Q25", "Q50", "Q75", "Q90", 
                 "P_Q10", "P_Q25", "P_Q50", "P_Q75", "P_Q90")
  require(doParallel)
  cl <- makeCluster(p)
  registerDoParallel(cl)
  gene.names <- rownames(data)
  results <- foreach(i = 1:(nrow(data)-1), .combine = "rbind", 
                    .export=ls(envir=globalenv()), .packages='matrixStats') %dopar% {
    temp <- matrix(NA, nrow(data)-i, length(ret.value))
    cmp <- rep(NA, nrow(data)-i)
    n <- 0
    for(j in (i+1):nrow(data)) {
      n <- n+1
      cmp[n] <- paste0(gene.names[i], ".vs.", gene.names[j])
      if(nonexpr.filter) {
        idx <- rep(TRUE, ncol(data))
        if(quantile(data[i,], 0) == quantile(data[i,], nonexpr.q)) idx <- idx & data[i,] > quantile(data[i,], 0)
        if(quantile(data[j,], 0) == quantile(data[j,], nonexpr.q)) idx <- idx & data[j,] > quantile(data[j,], 0)
        if( sum(idx) < 10 ) next()
        cycGD.res <- cycG_dect(rbind(data[i,idx], data[j,idx]))
        set.seed(seed)
        perm.cycGD.res <- cycG_dect(rbind(data[i,idx], sample(data[j,idx])))
        R <- cor(data[i,idx], data[j,idx])
      } else {
        cycGD.res <- cycG_dect(rbind(data[i,], data[j,]))
        set.seed(seed)
        perm.cycGD.res <- cycG_dect(rbind(data[i,], sample(data[j,])))
        R <- cor(data[i,], data[j,])
      }
      pole1 <- cycGD.res[,1] / mean(cycGD.res[,1])
      perm.pole1 <- perm.cycGD.res[,1] / mean(perm.cycGD.res[,1])
      temp[n, 1] <- wilcox.test(pole1, perm.pole1, alternative="greater")$p.val #wilcox.p
      temp[n, 2]  <- ks.test(pole1, perm.pole1, alternative="less")$p.val       #ks.l.p
      temp[n, 3]  <- ks.test(pole1, perm.pole1, alternative="greater")$p.val    #ks.g.p
      temp[n, 4]  <- R ^ 2                                                      #Rsq
      temp[n, 5:9] <- quantile(pole1, c(0.1, 0.25, 0.5, 0.75, 0.9))
      temp[n, 10:14] <- quantile(perm.pole1, c(0.1, 0.25, 0.5, 0.75, 0.9))
    }
    rownames(temp) <- cmp
    temp 
  }
  colnames(results) <- ret.value
  results 
}
    
cand_scatter_plot <- function(cycG.result, expr.data, cols=c(2,3), topN=100) { 
  feat <- apply(log10(cycG.result[, cols, drop=F]), 1, sum) 
  cycG.result.ordered <- cycG.result[order(feat), ]
  gene.pairs <- rownames(cycG.result.ordered)[1:topN]
  g1 <- sapply(gene.pairs, function(x) unlist(strsplit(x, "[.]"))[1])
  g2 <- sapply(gene.pairs, function(x) unlist(strsplit(x, "[.]"))[3])
  for(i in 1:topN) {
    plot(expr.data[g1[i],], expr.data[g2[i],], pch=20, xlab=g1[i], ylab=g2[i])
  }
  cycG.result.ordered
}
               
cand_scatter_plot.I80P <- function(cycG.result, expr.data, q.col=c(5,9), topN=200) { 
  feat <- apply(cycG.result, 1, function(x) diff(x[q.col]) )
  cycG.result.ordered <- cycG.result[order(feat), ]
  gene.pairs <- rownames(cycG.result.ordered)[1:topN]
  g1 <- sapply(gene.pairs, function(x) unlist(strsplit(x, "[.]"))[1])
  g2 <- sapply(gene.pairs, function(x) unlist(strsplit(x, "[.]"))[3])
  for(i in 1:topN) {
    ce <- cor(expr.data[g1[i],], expr.data[g2[i],])
    plot(expr.data[g1[i],], expr.data[g2[i],], pch=20, xlab=g1[i], ylab=g2[i], main=paste("r =",signif(ce, 3)))
  }
  cycG.result.ordered
}

cand_scatter_plot.IQR <- function(cycG.result, expr.data, q.col=c(6,8), topN=200) { 
  feat <- apply(cycG.result, 1, function(x) diff(x[q.col]) )
  cycG.result.ordered <- cycG.result[order(feat), ]
  gene.pairs <- rownames(cycG.result.ordered)[1:topN]
  g1 <- sapply(gene.pairs, function(x) unlist(strsplit(x, "[.]"))[1])
  g2 <- sapply(gene.pairs, function(x) unlist(strsplit(x, "[.]"))[3])
  for(i in 1:topN) {
    plot(expr.data[g1[i],], expr.data[g2[i],], pch=20, xlab=g1[i], ylab=g2[i])
  }
  cycG.result.ordered
}
       
cand_scatter_plot.IQR_perm <- function(cycG.result, expr.data, q.col=c(6,8), perm.q.col=c(11,13), topN=200) { 
  feat <- apply(cycG.result, 1, function(x) diff(x[q.col]) / diff(x[perm.q.col]) )
  cycG.result.ordered <- cycG.result[order(feat), ]
  gene.pairs <- rownames(cycG.result.ordered)[1:topN]
  g1 <- sapply(gene.pairs, function(x) unlist(strsplit(x, "[.]"))[1])
  g2 <- sapply(gene.pairs, function(x) unlist(strsplit(x, "[.]"))[3])
  for(i in 1:topN) {
    plot(expr.data[g1[i],], expr.data[g2[i],], pch=20, xlab=g1[i], ylab=g2[i])
  }
  cycG.result.ordered
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
