source("https://raw.githubusercontent.com/sunlightwang/Rfuns/master/theme.R")
library(limma)
library(goseq)

# GOSeq
run_goseq <- function(DEgenelist, Allgenelist, genome=c("hg19", "mm10"), geneID=c("geneSymbol"),
                      test.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"), geneLengthCorrect=FALSE,
                      minSize=10, maxSize=3000, gs_enrich_plot=T, padj_cutoff=0.05, topN=20, 
                      enrichment.limit=c(1,8), gene.list=FALSE, FDR=TRUE, pval_cutoff=0.05) {
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

  valid <- pvals$numInCat >= minSize & pvals$numInCat <= maxSize
  pvals <- pvals[valid, ]
  if(test.cats == "KEGG") {
    library(KEGG.db)
    xx <- as.list(KEGGPATHID2NAME)
    df <- data.frame(term=unlist(xx[pvals$category]), enrichment=pvals$numDEInCat/pvals$numInCat/(fg.n/bg.n),
                     p.value = -log10(pvals$over_represented_pvalue), 
                     FDR=-log10(p.adjust(pvals$over_represented_pvalue, method="fdr")), geneSet=pvals$category)
  } else {
    df <- data.frame(term=substr(pvals$term, 1, 51), enrichment=pvals$numDEInCat/pvals$numInCat/(fg.n/bg.n),
                     p.value = -log10(pvals$over_represented_pvalue), 
                     FDR=-log10(p.adjust(pvals$over_represented_pvalue, method="fdr")), geneSet=pvals$category)
  }
  if(gene.list) {
    GS2gene <- do.call(rbind, lapply(GS.gene.map(DEgenelist, genome=genome, geneID=geneID, GS.cats=test.cats, return.ID=TRUE), 
                                     function(x) paste(x, collapse=",")))                                     
    gene.list <- GS2gene[match(as.character(df$geneSet), rownames(GS2gene)),]                                
    df <- data.frame(df, gene.list=gene.list)
  }
  if(gs_enrich_plot) {
    df$enrichment [ df$enrichment > enrichment.limit[2] ] <- enrichment.limit[2]
    df <- df[1:topN,]
    if(FDR) {
      p <- ggplot(df, aes(x=reorder(term, FDR), y=FDR)) + geom_bar(aes(fill=enrichment),stat="identity") +
        coord_flip() + ylab("-log10(FDR)") + xlab("") + theme_Publication() + 
        scale_fill_gradient2(low="yellow", mid="orange", high="red", midpoint=2, limits=enrichment.limit) + 
        geom_hline(yintercept=-log10(padj_cutoff), color="grey50", linetype="dashed") 
      return(p) 
    } else {
      p <- ggplot(df, aes(x=reorder(term, p.value), y=p.value)) + geom_bar(aes(fill=enrichment),stat="identity") +
        coord_flip() + ylab("-log10(P value)") + xlab("") + theme_Publication() + 
        scale_fill_gradient2(low="yellow", mid="orange", high="red", midpoint=2, limits=enrichment.limit) + 
        geom_hline(yintercept=-log10(pval_cutoff), color="grey50", linetype="dashed") 
      return(p)
    }
  }
  if(FDR) {return(subset(df, FDR > -log10(padj_cutoff)))}
  if(!FDR) {return(subset(df, p.value > -log10(pval_cutoff)))}
}

# pvals.1 <-  run_goseq(DEgenelist,  Allgenelist, "hg19", "geneSymbol", "GO:BP", TRUE)
# pvals.2 <-  run_goseq(DEgenelist,  Allgenelist, "hg19", "geneSymbol", "GO:BP", FALSE)

### CAMERA
run_camera <- function(expr_log2, sample.cat=c(1,1,2,2), genome=c("hg19", "mm10"), geneID=c("geneSymbol"),
                       test.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"), inter_gene_cor=TRUE, FDR=0.1, 
                       minSize=10, maxSize=3000, gs_enrich_plot=T, topN=20, plot.nrow=4, plot.ncol=5, 
                       design = cbind(Intercept=1,Group=as.numeric(as.factor(sample.cat))-1), ... ) {
  require(goseq)
  reversemapping <- function(map) {
    tmp=unlist(map,use.names=FALSE)
    names(tmp)=rep(names(map),times=as.numeric(summary(map)[,1]))
    return(split(names(tmp),as.vector(tmp)))
  }
  gene2cat <- getgo(rownames(expr_log2), genome, geneID, fetch.cats=test.cats)
  names(gene2cat)=rownames(expr_log2)
  cat2gene=reversemapping(gene2cat)
  gene2cat=reversemapping(cat2gene)
  
  cat2genes.idx0 <- ids2indices(cat2gene, rownames(expr_log2), remove.empty=TRUE)
  geneset.size <- lapply(cat2genes.idx0, length)
  cat2genes.idx <- cat2genes.idx0[geneset.size >= minSize & geneset.size <= maxSize]
  if(test.cats == "KEGG") {
    require(KEGG.db)
    xx <- as.list(KEGGPATHID2NAME)
    names(cat2genes.idx) <- unlist(xx[names(cat2genes.idx)])
  } else { 
    require(GO.db)
    names(cat2genes.idx) <- select(GO.db, keys = names(cat2genes.idx), columns = c("TERM"))[,2]
    cat2genes.idx <- cat2genes.idx[! is.na(names(cat2genes.idx))]
  }
  inter.gene.cor <- ifelse(inter_gene_cor, 0.01, 0)

  camera.rst <- camera(expr_log2, cat2genes.idx, design, inter.gene.cor=inter.gene.cor, ...)
  N <- sum(camera.rst$FDR < FDR)
  camera.rst.topN <- camera.rst[1:min(N,topN), ]
  write.table(camera.rst.topN, sep="\t", quote=F)
  
  if(gs_enrich_plot) {
    n_page <- plot.nrow * plot.ncol
    n_plot <- nrow(camera.rst.topN) 
    idx.split <- split(1:n_plot, floor((1:n_plot-1)/n_page))
    p.list <- lapply(idx.split, function(xx) {
      camera.rst.topN.GS <- rownames(camera.rst.topN[xx,])
      cat2genes.idx.topN <- cat2genes.idx[camera.rst.topN.GS]
      ## sample 
      sample.mean_expr <- do.call(rbind, lapply(cat2genes.idx.topN, function(x) colMeans(expr_log2[x,])))
      colnames(sample.mean_expr) <- sample.cat
      require(reshape2)
      df <- melt(sample.mean_expr)
      colnames(df) <- c("ID","type","expr")
      p <- ggplot(df, aes(x=type, y=expr)) + geom_boxplot(aes(color=type), size=0.5) + 
               facet_wrap(~ ID, scales="free",nrow=plot.nrow,ncol=plot.ncol, labeller=label_wrap_gen(width=30)) + 
               theme_Publication() + theme(strip.text = element_text(face="plain", size=rel(0.5)))
      if( mean(table(sample.cat)) < 30 ) {
        p <- p + geom_point(size=0.5) 
      } 
      return(p)
      ## genes       
      #s1.idx <- design[,2] == unique(design[,2])[1]
      #s2.idx <- design[,2] == unique(design[,2])[2]
      #lapply(cat2genes.idx.topN, function(x) cbind(rowMeans(expr_log2[x,s1.idx]),rowMeans(expr_log2[x,s2.idx])))
    })      
    return(p.list)
  } else { 
    return(camera.rst.topN)
  }
}

                                  
GS.gene.map <- function(allGenes, genome=c("hg19", "mm10"), geneID=c("geneSymbol"),
                     GS.cats=c("GO:CC", "GO:BP", "GO:MF", "KEGG"), return.ID=FALSE) {
  require(goseq)
  reversemapping <- function(map) {
    tmp=unlist(map,use.names=FALSE)
    names(tmp)=rep(names(map),times=as.numeric(summary(map)[,1]))
    return(split(names(tmp),as.vector(tmp)))
  }
  gene2cat <- getgo(allGenes, genome, geneID, fetch.cats=GS.cats)
  names(gene2cat)=allGenes
  cat2gene=reversemapping(gene2cat)
  gene2cat=reversemapping(cat2gene)
  if(return.ID) { 
    return(cat2gene)
  }
  if(GS.cats == "KEGG") {
    require(KEGG.db)
    xx <- as.list(KEGGPATHID2NAME)
    names(cat2gene) <- unlist(xx[names(cat2gene)])
  } else { 
    require(GO.db)
    names(cat2gene) <- select(GO.db, keys = names(cat2gene), columns = c("TERM"))[,2]
  }
  cat2gene <- cat2gene[! is.na(names(cat2gene))]
  cat2gene
  #list(cat2gene=cat2gene, gene2cat=gene2cat)
}
                                  
