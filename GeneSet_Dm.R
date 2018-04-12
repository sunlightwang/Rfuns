source("https://raw.githubusercontent.com/sunlightwang/Rfuns/master/theme.R")
library(GO.db)
library(AnnotationDbi)
library(org.Dm.eg.db)

loadGO <- function(geneIDtype=c("FLYBASE", "SYMBOL","ENSEMBL"), 
                   type=c("BP", "MF", "CC")) {
  type <- match.arg(type, c("BP", "MF", "CC"))
  GoTerms <- as.list(GOTERM)
  GoVect <- unlist(lapply(GoTerms, function(X) GOID(X)))
  GoDesc <- unlist(lapply(GoTerms, function(X) Term(X)))
  annots <- select(org.Dm.eg.db, keys=GoVect,
                   columns=geneIDtype, keytype="GOALL") # GOALL includes GO child notes

  # rm NAs
  annots <- annots[!is.na(annots[,4]),]
  if(length(type) == 1) {
    annots <- annots[ annots$ONTOLOGY == type,]
  }
  data.frame(GeneSet=annots$GOALL, Gene=annots[,4], Desc=GoDesc[annots$GO], annots[,2:3])
}

GS_enrich <- function(GeneList, bgGeneList=NULL, annots, padj_cutoff=0.05, minHit=5, minBgHit=10, maxBgHit=3000, plot=T, maxPlotTerm=25) { 
  uniqGeneID <- unique(annots$Gene)
  if(is.null(bgGeneList)) {bgGeneList <- uniqGeneID}
  #totalN <- length(intersect(bgGeneList, uniqGeneID))
  totalN <- length(bgGeneList)
  SelM <- length(intersect(GeneList, uniqGeneID))
  nm <- do.call(rbind, 
                tapply(1:nrow(annots), as.factor(annots$GeneSet), 
                       function(rows) {
                         list_n <- length(intersect(bgGeneList, annots$Gene[rows]))
                         sel_m <- length(intersect(GeneList, annots$Gene[rows]))
                         enrich_fold <- sel_m / list_n / (SelM / totalN)
                         return(c(list_n, sel_m, enrich_fold))
                       }) )
  nm <- nm[ nm[,1] >= minBgHit & nm[,1] <= maxBgHit & nm[,2] >= minHit, , drop=F]
  if(nrow(nm) < 1) {return(NULL)}
  colnames(nm) <- c("bgHit", "Hit", "EnrichFold")
  ## white ball: in list
  ## black ball: out of list
  m <- SelM # total white balls
  n <- totalN - SelM # total black balls 
  q <- nm[,2] # selected white balls
  k <- nm[,1] # total selected bass
  pval <- phyper(q, m, n, k, lower.tail=FALSE) + dhyper(q, m, n, k)
  padj <- p.adjust(pval, method="fdr")
  result <- data.frame(nm, pval=pval,padj=padj)
  result <- result[order(padj),]
  desc <- unique(data.frame(annots$GeneSet, annots$Desc))
  rownames(desc) <- desc$annots.GeneSet
  result <- data.frame(result, desc=substr(as.character(noquote(desc[rownames(result),2])), 1, 51))
  result.sig <- result[result$padj < padj_cutoff, ]
  if(plot) { 
    if( nrow(result.sig) > maxPlotTerm ) { result.sig <- result.sig[1:maxPlotTerm,] }
    margin_add <- 0
    if( nrow(result.sig) < maxPlotTerm ) { margin_add <- (maxPlotTerm - nrow(result.sig)) * 5 / 3 }
    p <- ggplot(result.sig, aes(x=reorder(desc, -padj), y=-log10(padj))) + geom_bar(aes(fill=EnrichFold), stat="identity") +
          coord_flip() + ylab("-log10(FDR)") + xlab("") + theme_Publication() + scale_fill_continuous(low="yellow", high="red") + 
          geom_hline(yintercept=-log10(padj_cutoff), color="grey50", linetype="dashed") + 
          theme(plot.margin=unit(c(10+margin_add,5,5+margin_add,5),"mm"))
    return(p)
  } else { 
    return(result.sig)
  }
}               
