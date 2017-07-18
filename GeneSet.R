library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

loadGO <- function(organism = c("mouse", "human"), geneIDtype=c("SYMBOL","ENSEMBL","MGI","REFSEQ","UNIGENE","UNIPROT"), 
                   type=c("BP", "MF", "CC")) {
  org <- match.arg(organism, c("mouse", "human"))
  type <- match.arg(type, c("BP", "MF", "CC"))
  GoTerms <- as.list(GOTERM)
  GoVect <- unlist(lapply(GoTerms, function(X) GOID(X)))
  GoDesc <- unlist(lapply(GoTerms, function(X) Term(X)))
  if(org == "mouse") {
    annots <- select(org.Mm.eg.db, keys=GoVect,
                     columns=geneIDtype, keytype="GOALL") # GOALL includes GO child notes
  } else {
    annots <- select(org.Hs.eg.db, keys=GoVect,
                     columns=geneIDtype, keytype="GOALL") # GOALL includes GO child notes
  }
  # rm NAs
  annots <- annots[!is.na(annots[,4]),]
  if(length(type) == 1) {
    annots <- annots[ annots$ONTOLOGY == type,]
  }
  data.frame(GeneSet=annots$GOALL, Gene=annots[,4], Desc=GoDesc[annots$GO], annots[,2:3])
}

GS_enrich <- function(GeneList, bgGeneList=NULL, annots, padj_cutoff=0.05, minHit=5, minBgHit=10) { 
  uniqGeneID <- unique(annots$Gene)
  if(is.null(bgGeneList)) {bgGeneList <- uniqGeneID}
  totalN <- length(intersect(bgGeneList, uniqGeneID))
  SelM <- length(intersect(GeneList, uniqGeneID))
  nm <- do.call(rbind, 
                tapply(1:nrow(annots), as.factor(annots$GeneSet), 
                       function(rows) {
                         list_n <- length(intersect(bgGeneList, annots$Gene[rows]))
                         sel_m <- length(intersect(GeneList, annots$Gene[rows]))
                         enrich_fold <- sel_m / list_n / (SelM / totalN)
                         return(c(list_n, sel_m, enrich_fold))
                       }) )
  nm <- nm[ nm[,1] >= minBgHit & nm[,2] >= minHit, , drop=F]
  if(nrow(nm) < 1) {return(NULL)}
  colnames(nm) <- c("bgHit", "Hit", "EnrichFold")
  ## white ball: in list
  ## black ball: out of list
  m <- SelM # total white balls
  n <- totalN - SelM # total black balls 
  q <- nm[,2] # selected white balls
  k <- nm[,1] # total selected bass
  pval <- phyper(q, m, n, k, lower.tail=FALSE)
  padj <- p.adjust(pval, method="fdr")
  result <- data.frame(nm, pval=pval,padj=padj)
  result <- result[order(padj),]
  desc <- unique(data.frame(annots$GeneSet, annots$Desc))
  rownames(desc) <- desc$annots.GeneSet
  result <- data.frame(result, desc=as.character(noquote(desc[rownames(result),2])))
  result[result$padj < padj_cutoff, ]
}                          
