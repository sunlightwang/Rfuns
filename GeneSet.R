library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

loadGO <- function(organism = c("mouse", "human"), geneIDtype=c("SYMBOL","ENSEMBL","MGI","REFSEQ","UNIGENE","UNIPROT"), 
                   type=c("BP", "MF", "CC")) {
  org <- match.arg(organism, c("mouse", "human"))
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
    annots <- annots[ annots$ONTOLOGY == "BP",]
  }
  data.frame(GeneSet=annots$GO, Gene=annots[,4], Desc=GoDesc[annots$GO], annots[,2:3])
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
