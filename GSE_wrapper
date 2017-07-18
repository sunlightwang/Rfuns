library(limma)
library(goseq)

run_goseq <- function(DEgenelist, Allgenelist, genome=c("hg19", "mm10"), geneID=c("geneSymbol"), 
                      test.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"), geneLengthCorrect=FALSE) {
  DEgenes <- rep(0, length(Allgenelist))
  names(DEgenes) <- Allgenelist
  DEgenes[ names(DEgenes) %in% DEgenelist] <- 1
  pwf <- nullp(DEgenes, genome, geneID)
  if(geneLengthCorrect) {
    pvals <- goseq(pwf, genome, geneID, test.cats=test.cats, method="Wallenius")
  } else { 
    pvals <- goseq(pwf, genome, geneID,  test.cats=test.cats, method="Hypergeometric")
  }
}

# pvals.1 <-  run_goseq(DEgenelist,  Allgenelist, "hg19", "geneSymbol", "GO:BP", TRUE)
# pvals.2 <-  run_goseq(DEgenelist,  Allgenelist, "hg19", "geneSymbol", "GO:BP", FALSE)
