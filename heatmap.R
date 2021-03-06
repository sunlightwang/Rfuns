# heatmap
hclust2 <- function(x, method="ward.D", ...)
  hclust(x, method=method, ...)
dist2 <- function(x, ...)
  as.dist(1-cor(t(x), method="spearman"))

heatmap_cluster <- function(data, design.2, n.row.class=4, min.var=1, scale="row",labRow=FALSE,cexRow=1, cexCol=1,margins=c(10,10),
                            distfun=dist2, hclustfun=hclust2,...) { 
#n.row.class cannot be too large (<8)
  
  require(gplots)

  n.cc <- length(levels(as.factor(design.2)))
  cc <- rainbow(n.cc, start=0.2, end=1)[as.numeric(as.factor(design.2))]
                          
  ### remove row with low variability
  idx <- apply(data, 1, var) > min.var
  
  ### perform clustering on rows
  cl.row <- hclustfun(distfun(data[idx,]))
  ### extract cluster assignments; i.e. k=4 (rows) 
  gr.row <- cutree(cl.row, n.row.class)
  ### require(RColorBrewer)
  require(RColorBrewer)
  cr <- brewer.pal(n.row.class, "Set1")
  heatmap.2(as.matrix(data[idx,]), distfun=distfun, hclustfun=hclustfun, col=bluered(64), sepwidth=c(0,0),
            reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
            scale=scale, ColSideColors=cc, RowSideColors=cr[gr.row], key=TRUE, symkey=FALSE, labRow=labRow,
            rowsep=0, density.info="none", cexRow=cexRow, cexCol=cexCol, margins=margins,  trace="none", srtCol=45,...)
}
# data <- matrix(rnorm(2000), ncol=10)
# design.2 <- rep(1:2,5)
# heatmap_cluster(data, design.2)

###
#  heatmap.2(cor(cpm.log2[gene.intersect,idx1], method="spearman"), dendrogram="both", distfun=distfun, hclustfun=hclustfun, 
#            col=colfunc(32), sepwidth=c(0,0),cexCol=1,cexRow=1,
#            reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean), ColSideColors=cc[gr.row], symkey=F, srtCol=45,
#            scale="none",density.info="none", trace="none", keysize=1, key.xlab="CC", key.title="", xlab="", margins=c(10,12))
#
