### this is an R function library for scRNAseq data analysis
### TODO: make an object throughout the analysis, and using flags to indicate analysis steps; one can use the predifined pipeline and can also customize their own pipelines 

library(ggplot2)
library(DESeq2)
library(matrixStats)
library(statmod)
library(tsne)
library(grid)
library(ggthemes)
library(gridExtra)
library(scales)
library(destiny)
source("https://raw.githubusercontent.com/sunlightwang/Rfuns/master/theme.R")

geneScatterplot <- function( x, y, xlab, ylab, col, xylim=c(-1,4.0)) {
  plot( NULL, xlim=xylim, ylim=xylim,
        xaxt="n", yaxt="n", xaxs="i", yaxs="i", asp=1,
        xlab=xlab, ylab=ylab )
  abline( a=-1, b=1, col = "lightgray", lwd=2 )
  abline( a=0, b=1, col = "lightgray", lwd=2 )
  abline( a=1, b=1, col = "lightgray", lwd=2 )
  abline( h=c(0,2,4,6), v=c(0,2,4,6), col = "lightgray", lwd=2 )
  points(
    ifelse( x > 0, log10(x), -.7 ),
    ifelse( y > 0, log10(y), -.7 ),
    pch=19, cex=.5, col = col )
  axis( 1, c( -.7, 0:6 ),
        c( "0", "1", "10", "100", expression(10^3), expression(10^4),
           expression(10^5), expression(10^6) ) )
  axis( 2, c( -.7, 0:6 ),
        c( "0", "1", "10", "100", expression(10^3), expression(10^4),
           expression(10^5), expression(10^6) ), las=2 )
  axis( 1, -.35, "//", tick=FALSE, line=-.7 )
  axis( 2, -.35, "\\\\", tick=FALSE, line=-.7 )
}

Count.norm <- function(counts, method=c("sizefactor", "mean", "total","none")) {
  method <- match.arg(method, c("sizefactor", "mean", "total", "none"))
  if(method == "sizefactor") {
    counts.sf <- estimateSizeFactorsForMatrix( counts )
  } else if(method == "mean") { 
    mean.counts <- colMeans(counts)
    counts.sf <- mean.counts / mean ( mean.counts )
  } else if(method == "total") { 
    counts.sf <- colSums(counts) 
  } else { 
    counts.sf <- rep(1, ncol(counts))
  }
  t( t(counts) / counts.sf )
}

gene.biovar <- function(ERCC.cnt, Gene.cnt, minBiolDisp=0.5^2, winsorize=T, outlier.rm=F, outlier.coef=5) {
  if( is.null(ERCC.cnt) ) { ERCC.cnt <- Gene.cnt }
  sf.ERCC <- estimateSizeFactorsForMatrix( ERCC.cnt )
  sf.Gene <- estimateSizeFactorsForMatrix( Gene.cnt )
  ERCC.cnt.norm <- t( t(ERCC.cnt) / sf.ERCC )
  Gene.cnt.norm <- t( t(Gene.cnt) / sf.Gene )
  if(winsorize) { 
    winsorization <- function(vec) { sec_max <- sort(vec,decreasing=T)[2]; vec[which.max(vec)] <- sec_max; vec}
    ERCC.cnt.norm <- t(apply(ERCC.cnt.norm, 1, winsorization))
    Gene.cnt.norm <- t(apply(Gene.cnt.norm, 1, winsorization))
  }
  means.ERCC <- rowMeans( ERCC.cnt.norm )
  vars.ERCC <- rowVars( ERCC.cnt.norm )
  cv2.ERCC <- vars.ERCC / means.ERCC^2
  if(outlier.rm) { 
    r <- apply(Gene.cnt.norm, 1, function(gene_cnt) {
      gene_cnt_log2 <- log2(1+gene_cnt)
      Upper <- quantile(gene_cnt_log2, 3/4)
      Lower <- quantile(gene_cnt_log2, 1/4)
      IQR <- Upper - Lower
      Upper.outlier <- Upper + outlier.coef * IQR
      Lower.outlier <- Lower - outlier.coef * IQR
      gene_cnt <- gene_cnt[gene_cnt_log2 > Lower.outlier & gene_cnt_log2 < Upper.outlier]
      c(mean(gene_cnt), var(gene_cnt))
    })
    means.Gene <- r[1,]
    vars.Gene <-  r[2,]
  } else {
    means.Gene <- rowMeans( Gene.cnt.norm )
    vars.Gene <- rowVars( Gene.cnt.norm )
  }
  cv2.Gene <- vars.Gene / means.Gene^2
  # minimum mean value 
  minMeanForFit <- unname( quantile( means.ERCC[ which( cv2.ERCC > .3 ) ], .8) )
  useForFit<- means.ERCC >= minMeanForFit
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means.ERCC[useForFit] ), cv2.ERCC[useForFit] )
  residual <- var( log( fitted.values(fit) ) - log( cv2.ERCC[useForFit] ) )
  total <- var( log( cv2.ERCC[useForFit] ) )
  # # explained variances
  explained <- 1 - residual / total
  print(paste("ERCC curve fit explained:", explained))
  ######
  xi <- mean( 1 / sf.ERCC )
  m <- ncol(Gene.cnt.norm)
  #psia1theta <- mean( 1 / sf.ERCC ) + ( coefficients(fit)["a1tilde"] - xi ) * mean( sf.ERCC / sf.Gene) ## from scLVM
  psia1theta <- mean( 1 / sf.Gene ) + coefficients(fit)["a1tilde"] * mean( sf.ERCC / sf.Gene)
  cv2th <- coefficients(fit)["a0"] + minBiolDisp + coefficients(fit)["a0"] * minBiolDisp
  testDenom <- ( means.Gene * psia1theta + means.Gene^2 * cv2th ) / ( 1 + cv2th/m )
  Gene.bio_var <- vars.Gene * (m-1) / testDenom
  return(Gene.bio_var)
} 


ERCC_noise_model <- function(ERCC.cnt, plot=T, normalization=c("sizefactor", "none", "total"), lowCV4fit=TRUE, 
                             winsorize=T, seq_effi=c(0.7, 0.5, 0.3)) {
  normalization <- match.arg(normalization, c("sizefactor", "none", "total"))
  if(normalization == "sizefactor") {
    sf.ERCC <- estimateSizeFactorsForMatrix( ERCC.cnt )
    ERCC.cnt.norm <- t( t(ERCC.cnt) / sf.ERCC )
  } else if (normalization == "none") { 
    sf.ERCC <- rep(1, ncol(ERCC.cnt))
    ERCC.cnt.norm <- ERCC.cnt
  } else { # total
    # mean.ERCC <- colMeans(ERCC.cnt)
    # sf.ERCC <- mean.ERCC / mean(mean.ERCC)
    sf.ERCC <- colSums(ERCC.cnt)
    ERCC.cnt.norm <- t( t(ERCC.cnt) / sf.ERCC )
  }
  if(winsorize) { 
    winsorization <- function(vec) { sec_max <- sort(vec,decreasing=T)[2]; vec[which.max(vec)] <- sec_max; vec}
    ERCC.cnt.norm <- t(apply(ERCC.cnt.norm, 1, winsorization))
  }
  means.ERCC <- rowMeans( ERCC.cnt.norm )
  vars.ERCC <- rowVars( ERCC.cnt.norm )
  cv2.ERCC <- vars.ERCC / means.ERCC^2
  if(lowCV4fit) {
    # minimum mean value 
    minMeanForFit <- unname( quantile( means.ERCC[ which( cv2.ERCC > 0.3 ) ], 0.8) )
    useForFit<- means.ERCC >= minMeanForFit
  } else {
    useForFit <- rep(TRUE, length(cv2.ERCC))
  }
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means.ERCC[useForFit] ), cv2.ERCC[useForFit] )
  residual <- var( log( fitted.values(fit) ) - log( cv2.ERCC[useForFit] ) )
  total <- var( log( cv2.ERCC[useForFit] ) )
  # # explained variances
  # 1 - residual / total
  if(plot & normalization == "total") {
    plot( means.ERCC, cv2.ERCC, log="xy", col=1+useForFit, main="", xlim = c( 1e-5, 1 ), ylim = c( 1e-3, 100),  xaxt="n", yaxt="n" )
    axis( 1, 10^(-5:1), c(expression(10^-5),expression(10^-4), "0.001", "0.01", "0.1", "1", "10"))
    axis( 2, 10^(-3:2), c("0.001", "0.01", "0.1", "1", "10" ,"100"), las=2 )
    abline( h=10^(-5:5), v=10^(-5:5), col="#D0D0D0", lwd=2 )
    xg <- 10^seq( -5, 5, length.out=101 )
    lines( xg, coefficients(fit)["a0"] + coefficients(fit)["a1tilde"]/xg )
    segments( means.ERCC[useForFit], cv2.ERCC[useForFit],
              means.ERCC[useForFit], fit$fitted.values, col="gray" )
    legend("bottomleft", legend=c(paste0("a0: ", signif(coefficients(fit)["a0"], 3)), 
                                  paste0("a1: ", signif(coefficients(fit)["a1tilde"], 3)),
                                  paste0("explained variances: ", signif(1 - residual / total, 3))), bty="n")
    ### theoretical model
    a <- (1 + seq_effi) * mean(1 / sf.ERCC) 
    null <- sapply(1:length(a), function(i) lines( xg, - a[i] + a[i]/xg, col=2+i, lty=1 ) )
    null <- sapply(1:length(a), function(i) lines( xg, a[i]/xg, col=2+i, lty=2 ) )
    legend("topright", c(paste0("Binomial+SeqEff=",seq_effi), paste0("Poisson+SeqEff=",seq_effi)), lty=rep(1:2, each=length(seq_effi)), 
           col=2+1:length(seq_effi), bty="n")
  }
  if(plot & normalization %in% c("sizefactor", "none")) {
    plot( means.ERCC, cv2.ERCC, log="xy", col=1+useForFit, main="", xlim = c( 1e-3, 10000 ), ylim = c( 1e-2, 100),  xaxt="n", yaxt="n" )
    axis( 1, 10^(-3:4), c("0.001", "0.01", "0.1", "1", "10", "100", "1000", expression(10^4)))
    axis( 2, 10^(-2:2), c("0.01", "0.1", "1", "10" ,"100"), las=2 )
    abline( h=10^(-5:5), v=10^(-5:5), col="#D0D0D0", lwd=2 )
    xg <- 10^seq( -5, 5, length.out=101 )
    lines( xg, coefficients(fit)["a0"] + coefficients(fit)["a1tilde"]/xg )
    segments( means.ERCC[useForFit], cv2.ERCC[useForFit],
              means.ERCC[useForFit], fit$fitted.values, col="gray" )
    legend("bottomleft", legend=c(paste0("a0: ", signif(coefficients(fit)["a0"], 3)), 
                                  paste0("a1: ", signif(coefficients(fit)["a1tilde"], 3)),
                                  paste0("explained variances: ", signif(1 - residual / total, 3))), bty="n")
  }
  return()
} 

HVG.identifier <- function(ERCC.cnt, Gene.cnt, plot=T, normalization=c("sizefactor", "none", "mean"), minBiolDisp=0.5^2, 
                           padjcutoff=0.1, winsorize=T, topN=NULL, HVGnames=F) {
  if( is.null(ERCC.cnt) ) { ERCC.cnt <- Gene.cnt }
  normalization <- match.arg(normalization, c("sizefactor", "none", "mean"))
  if(normalization == "sizefactor") {
    sf.ERCC <- estimateSizeFactorsForMatrix( ERCC.cnt )
    sf.Gene <- estimateSizeFactorsForMatrix( Gene.cnt )
    ERCC.cnt.norm <- t( t(ERCC.cnt) / sf.ERCC )
    Gene.cnt.norm <- t( t(Gene.cnt) / sf.Gene )
  } else if (normalization == "none") { 
    sf.ERCC <- rep(1, ncol(ERCC.cnt))
    sf.Gene <- rep(1, ncol(Gene.cnt))
    ERCC.cnt.norm <- ERCC.cnt
    Gene.cnt.norm <- Gene.cnt
  } else { 
    mean.ERCC <- colMeans(ERCC.cnt)
    mean.Gene <- colMeans(Gene.cnt)
    sf.ERCC <- mean.ERCC / mean(mean.ERCC)
    sf.Gene <- mean.Gene / mean(mean.Gene)
    ERCC.cnt.norm <- t( t(ERCC.cnt) / sf.ERCC )
    Gene.cnt.norm <- t( t(Gene.cnt) / sf.Gene )
  }
  
  if(winsorize) { 
    winsorization <- function(vec) { sec_max <- sort(vec,decreasing=T)[2]; vec[which.max(vec)] <- sec_max; vec}
    ERCC.cnt.norm <- t(apply(ERCC.cnt.norm, 1, winsorization))
    Gene.cnt.norm <- t(apply(Gene.cnt.norm, 1, winsorization))
  }
  means.ERCC <- rowMeans( ERCC.cnt.norm )
  vars.ERCC <- rowVars( ERCC.cnt.norm )
  cv2.ERCC <- vars.ERCC / means.ERCC^2
  means.Gene <- rowMeans( Gene.cnt.norm )
  vars.Gene <- rowVars( Gene.cnt.norm )
  cv2.Gene <- vars.Gene / means.Gene^2
  # minimum mean value 
  minMeanForFit <- unname( quantile( means.ERCC[ which( cv2.ERCC > .3 ) ], .8) )
  useForFit<- means.ERCC >= minMeanForFit
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means.ERCC[useForFit] ), cv2.ERCC[useForFit] )
  residual <- var( log( fitted.values(fit) ) - log( cv2.ERCC[useForFit] ) )
  total <- var( log( cv2.ERCC[useForFit] ) )
  # # explained variances
  # 1 - residual / total
  if(plot) {
    plot( means.ERCC, cv2.ERCC, log="xy", col=1+useForFit, main="", xlim = c( 1e-2, 1e5 ), ylim = c( 0.01, 100) )
    xg <- 10^seq( -2, 5, length.out=100 )
    lines( xg, coefficients(fit)["a0"] + coefficients(fit)["a1tilde"]/xg )
    segments( means.ERCC[useForFit], cv2.ERCC[useForFit],
              means.ERCC[useForFit], fit$fitted.values, col="gray" )
    legend("bottomleft", legend=c(paste0("a0: ", signif(coefficients(fit)["a0"], 3)), 
                                  paste0("a1tilde: ", signif(coefficients(fit)["a1tilde"], 3)),
                                  paste0("explained variances: ", signif(1 - residual / total, 3))), bty="n")
  }
  ######
  xi <- mean( 1 / sf.ERCC )
  m <- ncol(Gene.cnt.norm)
  #psia1theta <- mean( 1 / sf.ERCC ) + ( coefficients(fit)["a1tilde"] - xi ) * mean( sf.ERCC / sf.Gene) ## from scLVM
  psia1theta <- mean( 1 / sf.Gene ) + coefficients(fit)["a1tilde"] * mean( sf.ERCC / sf.Gene)
  cv2th <- coefficients(fit)["a0"] + minBiolDisp + coefficients(fit)["a0"] * minBiolDisp
  testDenom <- ( means.Gene * psia1theta + means.Gene^2 * cv2th ) / ( 1 + cv2th/m )
  if(is.null(topN)) { 
    p <- 1 - pchisq( vars.Gene * (m-1) / testDenom, m-1 )
    padj <- p.adjust( p, "BH" )
    HVG.stat <- table( padj < padjcutoff & !is.nan(padj))
    HVG <- names(padj) [padj < padjcutoff & !is.nan(padj)]
  } else { 
    Gene.bio_var <- vars.Gene * (m-1) / testDenom
    HVG <- names(sort(Gene.bio_var, decreasing=T)[1:topN])
    HVG.stat <- table(names(means.Gene) %in% HVG)
  }
  if(plot) {
    plot( NULL, xaxt="n", yaxt="n",
          log="xy", xlim = c( 1e-2, 1e5 ), ylim = c( 0.01, 100 ),
          xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)" )
    axis( 1, 10^(-2:5), c( "0.01", "0.1", "1", "10", "100", "1000",
                           expression(10^4), expression(10^5) ) )
    axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10" ,"100"), las=2 )
    abline( h=10^(-2:1), v=10^(-2:5), col="#D0D0D0", lwd=2 )
    # Plot the plant genes, use a different color if they are highly variable
    if(is.null(topN)) {
      points( means.Gene, cv2.Gene, pch=20, cex=.2,
            col = ifelse( padj < padjcutoff, "#C0007090", "#70500040" ) )
    } else {
      points( means.Gene, cv2.Gene, pch=20, cex=.2,
            col = ifelse( names(means.Gene) %in% HVG, "#C0007090", "#70500040" ) )
    }
    if(HVGnames) { 
      text(means.Gene[HVG], cv2.Gene[HVG], label=HVG, cex=0.5)
    }
    # Add the technical noise fit, as before
    xg <- 10^seq( -2, 6, length.out=1000 )
    lines( xg, coefficients(fit)["a1tilde"] / xg + coefficients(fit)["a0"], col="#FF000080", lwd=3 )
    # Add a curve showing the expectation for the chosen biological CV^2 thershold
    lines( xg, psia1theta/xg + coefficients(fit)["a0"] + minBiolDisp,
           lty="dashed", col="#C0007090", lwd=3 )
    # Add the normalised ERCC points
    points( means.ERCC, cv2.ERCC, pch=20, cex=1, col="#0060B8A0" )
    legend("bottomleft", legend=paste0("HVG: ", HVG.stat[2]," out of ", sum(HVG.stat)), bty="n")
  }
  return(HVG)
} 

PCA.analysis <- function(Gene.cnt.scaled, plot=T, pca.perm.n=100, pca.padj.cutoff=0.01, plot_ngene=0, plot_nrow=3, permute=TRUE) {
  pca.real <- prcomp(t(Gene.cnt.scaled))
  pca.explained <- pca.real$sdev^2 / sum(pca.real$sdev^2) * 100
  if(plot) {
    # plot
    type <- factor(celltypes(colnames(Gene.cnt.scaled)))
    #pca.rst <- data.frame(PC1 = t(t(pca.real$rotation[, 1]) %*% Gene.cnt.scaled), PC2 = t(t(pca.real$rotation[, 2]) %*% Gene.cnt.scaled), type=type)
    pca.dist <- apply(as.matrix(dist(pca.real$x[,1:2], diag = T, upper = T)),1,function(x) min(x[x!=0]))
    pca.text <- pca.dist > quantile(pca.dist,0.95)
    pca.rst <- data.frame(PC1=pca.real$x[,1], PC2=pca.real$x[,2], names=colnames(Gene.cnt.scaled), type=type, pca.text = pca.text)
    p <- ggplot(pca.rst, aes(PC1, PC2, color = type)) + geom_point() + scale_colour_Rainbow() + theme_Publication() + 
      #  geom_text(data=subset(pca.rst, pca.text==T), aes(x=PC1, y=PC2, label=names)) + 
      xlab(paste0("PC1 (", signif(pca.explained[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pca.explained[2], 3), "%)"))
    print(p)
    p <- ggplot(pca.rst, aes(PC1, PC2, color = type)) + geom_text(aes(label=names), size=2) + scale_colour_Rainbow() + theme_Publication() 
    print(p)
  }
  if(plot & plot_ngene > 0) {
    ### for individual genes
    ngene <- plot_ngene
    nrowplot <- plot_nrow
    geneSel <- names(sort(pca.real$rotation[,1], decreasing = F)[1:ngene])
    pl <- lapply(1:ngene, function(x) 
      ggplot(pca.rst, aes(PC1, PC2, color = Gene.cnt.scaled[geneSel[x],])) + geom_point(size=1) + theme_graphOnly() + 
        scale_colour_gradientn(colors = topo.colors(10), guide=guide_colorbar(barheight=unit(3,"cm"))) + labs(color=geneSel[x]))
    p <- marrangeGrob(pl, nrow=nrowplot, ncol=nrowplot, top="major genes in PC1")
    print(p)
    geneSel <- names(sort(pca.real$rotation[,1], decreasing = T)[1:ngene])
    pl <- lapply(1:ngene, function(x) 
      ggplot(pca.rst, aes(PC1, PC2, color = Gene.cnt.scaled[geneSel[x],])) + geom_point(size=1) + theme_graphOnly() + 
        scale_colour_gradientn(colors = topo.colors(10), guide=guide_colorbar(barheight=unit(3,"cm"))) + labs(color=geneSel[x]))
    p <- marrangeGrob(pl, nrow=nrowplot, ncol=nrowplot, top="major genes in PC1")
    print(p)
    geneSel <- names(sort(pca.real$rotation[,2], decreasing = F)[1:ngene])
    pl <- lapply(1:ngene, function(x) 
      ggplot(pca.rst, aes(PC1, PC2, color = Gene.cnt.scaled[geneSel[x],])) + geom_point(size=1) + theme_graphOnly() + 
        scale_colour_gradientn(colors = topo.colors(10), guide=guide_colorbar(barheight=unit(3,"cm"))) + labs(color=geneSel[x]))
    p <- marrangeGrob(pl, nrow=nrowplot, ncol=nrowplot, top="major genes in PC2")
    print(p)
    geneSel <- names(sort(pca.real$rotation[,2], decreasing = T)[1:ngene])
    pl <- lapply(1:ngene, function(x) 
      ggplot(pca.rst, aes(PC1, PC2, color = Gene.cnt.scaled[geneSel[x],])) + geom_point(size=1) + theme_graphOnly() + 
        scale_colour_gradientn(colors = topo.colors(10), guide=guide_colorbar(barheight=unit(3,"cm"))) + labs(color=geneSel[x]))
    p <- marrangeGrob(pl, nrow=nrowplot, ncol=nrowplot, top="major genes in PC2")
    print(p)
  }
  if(permute) {
    #### PCA permutation
    permute.times <- pca.perm.n
    Gene.permute <- function(data) do.call(rbind, sapply(1:nrow(data), function(x) sample(data[x,]), simplify=F)) 
    pca.permuted <- sapply(1:permute.times, function(x){
      ev <- prcomp(t(Gene.permute(Gene.cnt.scaled)))$sdev ^ 2
      max( ev / sum(ev) )
    })
    pca.pval <- sapply(1:length(pca.real$sdev), function(x) {
      sum(pca.real$sdev[x]^2/sum(pca.real$sdev^2) < pca.permuted) / length(pca.permuted)
    })
    pca.padj <- p.adjust(pca.pval, method="BH")
    pca.padj.cutoff <- pca.padj.cutoff
    table(pca.padj < pca.padj.cutoff)
    pca_gene <- t(pca.real$rotation[, pca.padj < pca.padj.cutoff]) %*% Gene.cnt.scaled
    return(pca_gene) 
  }
  return(pca.real)
}

tSNE.analysis <- function(Gene.cnt.scaled, perplexity=30, max_iter=2000, try_times=100, plot=T, gene_expr=Gene.cnt.scaled, display=c(), plot_nrow=3, ...) {
  library(Rtsne)
  tsne.out <- sapply(1:try_times, function(x) Rtsne(t(Gene.cnt.scaled), dims=2, pca=F, perplexity=perplexity, max_iter=max_iter, ...))
  best <- which.min(sapply(1:try_times, function(x) tail(tsne.out[,x]$itercosts, 1)))
  type <- factor(celltypes(colnames(Gene.cnt.scaled)))
  tsne.rst <- data.frame(tSNE.1 = tsne.out[,best]$Y[,1], tSNE.2 = tsne.out[,best]$Y[,2], type=type, names=colnames(Gene.cnt.scaled))
  if(plot) {
    p <- ggplot(tsne.rst, aes(tSNE.1, tSNE.2, color = type)) + geom_point() + scale_colour_Rainbow() + theme_Publication()
    print(p)
    p <- ggplot(tsne.rst, aes(tSNE.1, tSNE.2, color = type)) + geom_text(aes(label=names), size=2) + scale_colour_Rainbow() + theme_Publication()
    print(p)
  }
  display <- intersect(display, rownames(gene_expr))
  if(plot & !is.null(display) & length(display) > 0) {
    ### for individual genes
    nrowplot <- plot_nrow
    ngene <- length(display)
    pl <- lapply(1:ngene, function(x) 
      ggplot(tsne.rst, aes(tSNE.1, tSNE.2, color = gene_expr[display[x],])) + geom_point(size=1) + theme_graphOnly() + 
        scale_colour_gradientn(colors = topo.colors(10), guide=guide_colorbar(barheight=unit(3,"cm"))) + labs(color=display[x]))
    p <- marrangeGrob(pl, nrow=nrowplot, ncol=nrowplot, top="")
    print(p)
  }
  return(tsne.rst)
}

diffusionmap.analysis <- function(Gene.cnt.scaled, dims = 1:3, dist.method=c("euclidean", "cosine", "rankcor"), 
                                  sigma = "local", n_local = 5:7, density_norm = TRUE, plot=T, 
                                  gene_expr=Gene.cnt.scaled, display=c(), plot_nrow=3, ...) { 
  dfmap <- DiffusionMap(t(Gene.cnt.scaled), sigma=sigma, density_norm=density_norm, distance=distance, n_local=n_local, ...)
  type <- factor(celltypes(colnames(Gene.cnt.scaled)))
  dfmap.rst <- data.frame(DC1=dfmap$DC1, DC2=dfmap$DC2, DC3=dfmap$DC3, type=type, names=colnames(Gene.cnt.scaled))
  if(plot) { 
    p <- plot_ly(dfmap.rst, x = ~DC1, y = ~DC2, z = ~DC3, color=~type, colors = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#a0a013","#dd4444")) %>%
      add_markers() %>%
      layout(scene = list(xaxis = list(title = 'DC1'),
                          yaxis = list(title = 'DC2'),
                          zaxis = list(title = 'DC3')))
    htmlwidgets::saveWidget(as_widget(p), "diffusionmap.html")
    plot(dfmap, col=type, pch=20)
  }
  display <- intersect(display, rownames(gene_expr))
  if(plot & !is.null(display) & length(display) > 0) {
    ### for individual genes
    nrowplot <- plot_nrow
    ngene <- length(display)
    pl <- lapply(1:ngene, function(x) {
      p <- plot_ly(dfmap.rst, x = ~DC1, y = ~DC2, z = ~DC3, color= gene_expr[display[x],], colors = rev(rainbow(3))) %>%
        add_markers() %>%
        layout(scene = list(xaxis = list(title = 'DC1'),
                            yaxis = list(title = 'DC2'),
                            zaxis = list(title = 'DC3')))
      htmlwidgets::saveWidget(as_widget(p), paste0("diffusionmap_", display[x],".html"))
      ### 
      ggplot(dfmap.rst, aes(DC1, DC2, color = gene_expr[display[x],])) + geom_point(size=1) + theme_graphOnly() + 
        scale_colour_gradientn(colors = topo.colors(10), guide=guide_colorbar(barheight=unit(3,"cm"))) + labs(color=display[x])
    }) 
    p <- marrangeGrob(pl, nrow=nrowplot, ncol=nrowplot, top="")
    print(p)
  }
  return(dfmap)
}
