source("https://raw.githubusercontent.com/sunlightwang/Rfuns/master/theme.R")
library(limma)
library(goseq)

# GOSeq
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
    if(test.cats == "KEGG") {
      library(KEGG.db)
      xx <- as.list(KEGGPATHID2NAME)
      df <- data.frame(term=unlist(xx[pvals$category]), enrichment=pvals$numDEInCat/pvals$numInCat/(fg.n/bg.n),
                     FDR=-log10(p.adjust(pvals$over_represented_pvalue)))[1:topN,]
    } else {
      df <- data.frame(term=pvals$term, enrichment=pvals$numDEInCat/pvals$numInCat/(fg.n/bg.n),
                     FDR=-log10(p.adjust(pvals$over_represented_pvalue, method="fdr")))[1:topN,]
    }
    p <- ggplot(df, aes(x=reorder(term, FDR), y=FDR)) + geom_bar(aes(fill=enrichment),stat="identity") +
    coord_flip() + ylab("-log10(FDR)") + xlab("") + theme_Publication() + scale_fill_continuous(low="yellow", high="red")
    return(p)
  }
  return(pvals)
}

# pvals.1 <-  run_goseq(DEgenelist,  Allgenelist, "hg19", "geneSymbol", "GO:BP", TRUE)
# pvals.2 <-  run_goseq(DEgenelist,  Allgenelist, "hg19", "geneSymbol", "GO:BP", FALSE)

### CAMERA
run_camera <- function(expr_log2, no.samples=c(2,2), genome=c("hg19", "mm10"), geneID=c("geneSymbol"),
                       test.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"), inter_gene_cor=TRUE, 
                       minSize=10, maxSize=3000, gs_enrich_plot=T, topN=20) {
  
  design <- cbind(Intercept=1,Group=c(rep(0,no.samples[1]), rep(1,no.samples[2])))
  require(goseq)
  reversemapping <- function(map) {
    tmp=unlist(map,use.names=FALSE)
    names(tmp)=rep(names(map),times=as.numeric(summary(map)[,1]))
    return(split(names(tmp),as.vector(tmp)))
  }
  gene2cat <- getgo(rownames(expr_log2), genome, geneID, fetch.cats=test.cats)
  cat2gene=reversemapping(gene2cat)
  gene2cat=reversemapping(cat2gene)
  
  cat2genes.idx0 <- ids2indices(cat2gene, rownames(expr_log2), remove.empty=TRUE)
  geneset.size <- lapply(cat2genes.idx0, length)
  cat2genes.idx <- cat2genes.idx0[geneset.size >= minSize & geneset.size <= maxSize]

  inter.gene.cor <- ifelse(inter_gene_cor, 0.01, 0)

  camera.rst <- camera(expr_log2, cat2genes.idx, design, inter.gene.cor=inter.gene.cor)
  camera.rst.topN <- camera.rst[1:topN, ]
  camera.rst.topN.GS <- rownames(camera.rst.topN)
  cat2genes.idx.topN <- cat2genes.idx[camera.rst.topN.GS]
}
