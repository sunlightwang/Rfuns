# heatmap

library(gplots)
hclust2 <- function(x, method="ward.D", ...)
  hclust(x, method=method, ...)
dist2 <- function(x, ...)
  as.dist(1-cor(t(x), method="spearman"))
cc <- rep(rainbow(ncol(log2.cnt.matrix)/2, start=0.2, end=1), 2)

### perform clustering on rows
cl.row <- hclust2(dist2(log2.cnt.matrix[idx,]))
### extract cluster assignments; i.e. k=8 (rows) 
gr.row <- cutree(cl.row, 6)
### require(RColorBrewer)
require(RColorBrewer)
cr <- brewer.pal(6, "Set1")
heatmap.2(log2.cnt.matrix[idx,], distfun=dist2, hclustfun=hclust2, col=bluered(64), sepwidth=c(0,0),
          reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
          scale="row", ColSideColors=cc, RowSideColors=cr[gr.row], key=TRUE, symkey=FALSE, labRow=FALSE,
          rowsep=0, density.info="none", cexRow=1, cexCol=1, margins=c(6,11),  trace="none", srtCol=45)

###
heatmap.2(cor(cpm.log2[gene.intersect,idx1], method="spearman"), dendrogram="both", distfun=dist, hclustfun=hclust2, col=colfunc(32), sepwidth=c(0,0),cexCol=1,cexRow=1,
          reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean), ColSideColors=cc[gr.row], symkey=F, srtCol=45,
          scale="none",density.info="none", trace="none", keysize=1, key.xlab="CC", key.title="", xlab="", margins=c(10,12))
#
