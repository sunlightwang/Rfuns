source("https://raw.githubusercontent.com/sunlightwang/Rfuns/master/theme.R")
library(limma)
library(goseq)

run_goseq <- function(DEgenelist, Allgenelist, genome=c("hg19", "mm10"), geneID=c("geneSymbol"),
                               test.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"), geneLengthCorrect=FALSE,
                               minSize=10, maxSize=3000, gs_enrich_plot=T, topN=20) {
  DEgenes <- rep(0, length(Allgenelist))
  names(DEgenes) <- Allgenelist
  DEgenes[ names(DEgenes) %in% DEgenelist] <- 1
  pwf <- nullp(DEgenes, genome, geneID, plot.fit=F)
  if(geneLengthCorrect) {
    pvals <- goseq(pwf, genome, geneID, test.cats=test.cats, method="Wallenius", use_genes_without_cat=T)
  } else {
    pvals <- goseq(pwf, genome, geneID,  test.cats=test.cats, method="Hypergeometric", use_genes_without_cat=T)
  }
  fg.n <- length(DEgenelist)
  bg.n <- length(Allgenelist)
  if(gs_enrich_plot) {
    valid <- pvals$numInCat >= minSize & pvals$numInCat <= maxSize
    pvals <- pvals[valid, ]
    df <- data.frame(term=pvals$term, enrichment=pvals$numDEInCat/pvals$numInCat/(fg.n/bg.n),
                     FDR=-log10(p.adjust(pvals$over_represented_pvalue)))[1:topN,]
    p <- ggplot(df, aes(x=reorder(term, FDR), y=FDR)) + geom_bar(aes(fill=enrichment),stat="identity") +
    coord_flip() + ylab("-log10(FDR)") + xlab("") + theme_Publication() + scale_fill_continuous(low="blue", high="red")
    return(p)
  }
  return(pvals)
}

# pvals.1 <-  run_goseq(DEgenelist,  Allgenelist, "hg19", "geneSymbol", "GO:BP", TRUE)
# pvals.2 <-  run_goseq(DEgenelist,  Allgenelist, "hg19", "geneSymbol", "GO:BP", FALSE)
