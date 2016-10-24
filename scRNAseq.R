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

theme_Publication <- function(base_size=14, base_family="Helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90, vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.line.x = element_line(color="black"),
           axis.line.y = element_line(color="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.size= unit(0.4, "cm"),
           legend.margin = unit(0.1, "cm"),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
}

theme_graphOnly <- function(base_size=14, base_family="Helvetica") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
           text = element_text(),
           axis.line=element_blank(),
           axis.text.x=element_blank(),
           axis.text.y=element_blank(),
           axis.ticks=element_blank(),
           axis.title.x=element_blank(),
           axis.title.y=element_blank(),
           panel.grid.major=element_blank(),
           panel.grid.minor=element_blank(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = "#f0f0f0"),
           panel.border = element_rect(colour = NA),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           legend.direction = "vertical",
           legend.key.size= unit(0.4, "cm"),
           legend.margin = unit(0.1, "cm"),
           legend.title = element_text(face="italic")
   ))
}


scale_fill_Publication <- function(...){
  discrete_scale("fill","Publication", manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_Publication <- function(...){
  discrete_scale("colour","Publication", manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

scale_colour_Rainbow <- function(size = 10, ...){
  discrete_scale("colour","Rainbow", manual_pal(values = rainbow(size)), ...)
}


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

Count.norm <- function(counts) {
  counts.sf <- estimateSizeFactorsForMatrix( counts )
  t( t(counts) / counts.sf )
}

HVG.identifier <- function(ERCC.cnt, Gene.cnt, plot=T, minBiolDisp=0.5^2, padjcutoff=0.1, topN=NULL) {
  sf.ERCC <- estimateSizeFactorsForMatrix( ERCC.cnt )
  sf.Gene <- estimateSizeFactorsForMatrix( Gene.cnt )
  ERCC.cnt.norm <- t( t(ERCC.cnt) / sf.ERCC )
  Gene.cnt.norm <- t( t(Gene.cnt) / sf.Gene )
  means.ERCC <- rowMeans( ERCC.cnt.norm )
  vars.ERCC <- rowVars( ERCC.cnt.norm )
  cv2.ERCC <- vars.ERCC / means.ERCC^2
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
    legend("bottomleft", legend=paste0("explained variances: ", signif(1 - residual / total, 3)), bty="n")
  }
 
  ######
  means.Gene <- rowMeans( Gene.cnt.norm )
  vars.Gene <- rowVars( Gene.cnt.norm )
  cv2.Gene <- vars.Gene / means.Gene^2
  
  xi <- mean( 1 / sf.ERCC )
  m <- ncol(Gene.cnt.norm)
  psia1theta <- mean( 1 / sf.ERCC ) +
    #psia1theta <- mean( 1 / sf.Gene ) + ### %%% modified by me
    ( coefficients(fit)["a1tilde"] - xi ) * mean( sf.ERCC / sf.Gene)
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
  HVG
} 


PCA.analysis <- function(Gene.cnt.norm, HVG, plot=T, pca.perm.n=100, pca.padj.cutoff=0.01, plot_ngene=0, plot_nrow=3) {
  center.scale <- function(data)  ( data - rowMeans(data) ) / sqrt(rowVars(data)) 
  Gene.cnt.scaled <- center.scale( log2(Gene.cnt.norm[HVG,]+1) )
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
}

tSNE.analysis <- function(Gene.cnt.scaled, perplexity=10, max_iter=2000, try_times=100, plot=T, gene_expr=c(), plot_nrow=3, ...) {
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
  if(!is.null(gene_expr)) {
    ### for individual genes
    nrowplot <- plot_nrow
    ngene <- length(gene_expr)
    pl <- lapply(1:ngene, function(x) 
      ggplot(tsne.rst, aes(tSNE.1, tSNE.2, color = Gene.cnt.scaled[gene_expr[x],])) + geom_point(size=1) + theme_graphOnly() + 
        scale_colour_gradientn(colors = topo.colors(10), guide=guide_colorbar(barheight=unit(3,"cm"))) + labs(color=gene_expr[x]))
    p <- marrangeGrob(pl, nrow=nrowplot, ncol=nrowplot, top="")
    print(p)
  }
}
