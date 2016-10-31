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
                     columns=geneIDtype, keytype="GO") 
  } else {
    annots <- select(org.Hs.eg.db, keys=GoVect,
                     columns=geneIDtype, keytype="GO")
  }
  # rm NAs
  annots <- annots[!is.na(annots[,4]),]
  if(length(type) == 1) {
    annots <- annots[ annots$ONTOLOGY == type,]
  }
  data.frame(GeneSet=annots$GO, Gene=annots[,4], Desc=GoDesc[annots$GO], annots[,2:3])
}

GS_enrich <- function(GeneList, bgGeneList, annots, padj_cutoff=0.05) { 
  uniqGeneID <- unique(annots$Gene)
  totalN <- length(intersect(bgGeneList, uniqGeneID))
  SelM <- length(intersect(GeneList, uniqGeneID))
  nm <- do.call(rbind, 
                tapply(1:nrow(annots), as.factor(annots$GeneSet), 
                       function(rows) {
                         list_n <- length(intersect(bgGeneList, annots$Gene))
                         sel_m <- length(intersect(GeneList, annots$Gene))
                         enrich_fold <- sel_m / list_n / (SelM / totalN)
                         return(c(list_n, sel_m, enrich_fold))
                       }) )
  nm <- nm[ nm[,1] >= 10 & nm[,2] >= 5, ]
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
  result <- data.frame(result, desc=as.character(noquote(desc[rownames(result),2)))
  result[result$padj < padj_cutoff, ]
}                          
                          
### OLD
quickGS <- function(selGeneList, geneIDtype = "ENSEMBL",  fullGeneList = "", 
                    organism = c("mouse", "human"), minGSsize = 5, topNSig = 10) {

  
  uniqGeneID <- unique(annots[,4])
  if (fullGeneList == "") {
    totalN <- length(uniqGeneID)
  } else {
    totalN <- length(intersect(fullGeneList, uniqGeneID))
  }
  SelM <- length(intersect(selGeneList, uniqGeneID))
  nm <- do.call(rbind, 
                tapply(1:nrow(annots), as.factor(annots$GO), 
                       function(rows) {
                         if(fullGeneList == "") {
                           list_n <- length(unique(annots[rows,4]))
                         } else { 
                           list_n <- length(intersect(fullGeneList, annots[rows,4]))
                         } 
                         sel_m <- length(intersect(selGeneList, annots[rows,4]))
                         return(c(list_n, sel_m))
                       }) )
  nm <- nm[ nm[,1] >= minGSsize, ]
  colnames(nm) <- c("GSsize", "Hit")
  ## white ball: in list
  ## black ball: out of list
  m <- SelM # total white balls
  n <- totalN - SelM # total black balls 
  q <- nm[,2] # selected white balls
  k <- nm[,1] # total selected bass
  pval <- phyper(q, m, n, k, lower.tail=FALSE)
  result <- cbind(nm, pval)
  result <- result[order(pval)[1:topNSig],]
  terms <- sapply(rownames(result), function(X) Term(GoTerms[[which(names(GoTerms)==X)]]) )
  cbind(result, noquote(terms))
}
