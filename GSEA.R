## this is modified from GSEA.R of Broad Institute

GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. 
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #   arg.ES: Location in gene.list where the peak running enrichment occurs (peak of the "mountain") 
  #   RES: Numerical vector containing the running enrichment score for all locations in the gene list 
  #   tag.indicator: Binary vector indicating the location of the gene sets (1's) in the gene list 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
  no.tag.indicator <- 1 - tag.indicator 
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector**alpha)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)      
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > - min.ES) {
    #      ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    #      ES <- min.ES
    ES <- signif(min.ES, digits=5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))    
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in 
  # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the 
  # observed one.
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  
  loc.vector <- vector(length=N, mode="numeric")
  peak.res.vector <- vector(length=Nh, mode="numeric")
  valley.res.vector <- vector(length=Nh, mode="numeric")
  tag.correl.vector <- vector(length=Nh, mode="numeric")
  tag.diff.vector <- vector(length=Nh, mode="numeric")
  tag.loc.vector <- vector(length=Nh, mode="numeric")
  
  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]
  
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  
  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else {
    tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(list(ES = ES))
  
}

GSEA.EnrichmentScore2 <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {  
  #
  # Computes the weighted GSEA score of gene.set in gene.list. It is the same calculation as in 
  # GSEA.EnrichmentScore but faster (x8) without producing the RES, arg.RES and tag.indicator outputs.
  # This call is intended to be used to asses the enrichment of random permutations rather than the 
  # observed one.
  # The weighted score type is the exponent of the correlation 
  # weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted). When the score type is 1 or 2 it is 
  # necessary to input the correlation vector with the values in the same order as in the gene list.
  #
  # Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   weighted.score.type: Type of score: weight: 0 (unweighted = Kolmogorov-Smirnov), 1 (weighted), and 2 (over-weighted)  
  #  correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  #
  # Outputs:
  #   ES: Enrichment score (real number between -1 and +1) 
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  
  loc.vector <- vector(length=N, mode="numeric")
  peak.res.vector <- vector(length=Nh, mode="numeric")
  valley.res.vector <- vector(length=Nh, mode="numeric")
  tag.correl.vector <- vector(length=Nh, mode="numeric")
  tag.diff.vector <- vector(length=Nh, mode="numeric")
  tag.loc.vector <- vector(length=Nh, mode="numeric")
  
  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]
  
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  
  if (weighted.score.type == 0) {
    tag.correl.vector <- rep(1, Nh)
  } else if (weighted.score.type == 1) {
    tag.correl.vector <- correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else if (weighted.score.type == 2) {
    tag.correl.vector <- correl.vector[tag.loc.vector]*correl.vector[tag.loc.vector]
    tag.correl.vector <- abs(tag.correl.vector)
  } else {
    tag.correl.vector <- correl.vector[tag.loc.vector]**weighted.score.type
    tag.correl.vector <- abs(tag.correl.vector)
  }
  
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  valley.res.vector <- peak.res.vector - tag.correl.vector
  max.ES <- max(peak.res.vector)
  min.ES <- min(valley.res.vector)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  
  return(list(ES = ES))
  
}


GSEA.HeatMapPlot <- function(V, row.names = F, col.labels, col.classes, col.names = F, main = " ", xlab=" ", ylab=" ") {
  #
  # Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
  #
  # The Broad Institute
  # SOFTWARE COPYRIGHT NOTICE AGREEMENT
  # This software and its documentation are copyright 2003 by the
  # Broad Institute/Massachusetts Institute of Technology.
  # All rights are reserved.
  #
  # This software is supplied without any warranty or guaranteed support
  # whatsoever. Neither the Broad Institute nor MIT can be responsible for
  # its use, misuse, or functionality.
  
  n.rows <- length(V[,1])
  n.cols <- length(V[1,])
  row.mean <- apply(V, MARGIN=1, FUN=mean)
  row.sd <- apply(V, MARGIN=1, FUN=sd)
  row.n <- length(V[,1])
  for (i in 1:n.rows) {
    if (row.sd[i] == 0) {
      V[i,] <- 0
    } else {
      V[i,] <- (V[i,] - row.mean[i])/(0.5 * row.sd[i])
    }
    V[i,] <- ifelse(V[i,] < -6, -6, V[i,])
    V[i,] <- ifelse(V[i,] > 6, 6, V[i,])
  }
  
  mycol <- c("#0000FF", "#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA", "#FF9DB0", "#FF7080", "#FF5A5A", "#FF4040", "#FF0D1D", "#FF0000") # blue-pinkogram colors. The first and last are the colors to indicate the class vector (phenotype). This is the 1998-vintage, pre-gene cluster, original pinkogram color map
  
  mid.range.V <- mean(range(V)) - 0.1
  heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
  heatm[1:n.rows,] <- V[seq(n.rows, 1, -1),]
  heatm[n.rows + 1,] <- ifelse(col.labels == 0, 7, -7)
  image(1:n.cols, 1:(n.rows + 1), t(heatm), col=mycol, axes=FALSE, main=main, xlab= xlab, ylab=ylab)
  
  if (length(row.names) > 1) {
    numC <- nchar(row.names)
    size.row.char <- 35/(n.rows + 5)
    size.col.char <- 25/(n.cols + 5)
    maxl <- floor(n.rows/1.6)
    for (i in 1:n.rows) {
      row.names[i] <- substr(row.names[i], 1, maxl)
    }
    row.names <- c(row.names[seq(n.rows, 1, -1)], "Class")
    axis(2, at=1:(n.rows + 1), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, font.axis=2, line=-1)
  }
  
  if (length(col.names) > 1) {
    axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
  }
  
  C <- split(col.labels, col.labels)
  class1.size <- length(C[[1]])
  class2.size <- length(C[[2]])
  axis(3, at=c(floor(class1.size/2),class1.size + floor(class2.size/2)), labels=col.classes, tick=FALSE, las = 1, cex.axis=1.25, font.axis=2, line=-1)
  
  return()
}

GSEA <- function(geneScores, geneSets, nperm = 1000, weighted.score.type = 1, 
  nom.p.val.threshold = -1, fwer.p.val.threshold = -1, fdr.q.val.threshold = 0.25, adjust.FDR.q.val = F, topgs = 10, 
  gs.size.threshold.min = 25, gs.size.threshold.max = 1000, gs.size.percent.max = 35%, random.seed = 123456) {
  
  set.seed(seed=random.seed, kind = NULL)
  adjust.param <- 0.5

  # geneScores
  geneScores <- geneScores
  geneIDs <- rownames(geneScores) 
  # gene set database TODO: overlapSize(geneSets, geneIDs) returns gene set overlapping size
  geneSet.size <- overlapSize(geneSets, geneIDs) # geneSets: list
  geneSets.flt <- geneSets[[ geneSet.size >= gs.size.threshold.min & geneSet.size <= min(gs.size.threshold.max, length(geneIDs) * gs.size.percent.max)]]
  
  # Read gene and gene set annotations if gene annotation file was provided
  all.gene.descs <- vector(length = N, mode ="character") 
  all.gene.symbols <- vector(length = N, mode ="character") 
  all.gs.descs <- vector(length = Ng, mode ="character") 
  
  if (is.data.frame(gene.ann)) {
    temp <- gene.ann
    a.size <- length(temp[,1])
    print(c("Number of gene annotation file entries:", a.size))  
    accs <- as.character(temp[,1])
    locs <- match(gene.labels, accs)
    all.gene.descs <- as.character(temp[locs, "Gene.Title"])
    all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
    rm(temp)
  } else  if (gene.ann == "") {
    for (i in 1:N) {
      all.gene.descs[i] <- gene.labels[i]
      all.gene.symbols[i] <- gene.labels[i]
    }
  } else {
    temp <- read.delim(gene.ann, header=T, sep=",", comment.char="", as.is=T)
    a.size <- length(temp[,1])
    print(c("Number of gene annotation file entries:", a.size))  
    accs <- as.character(temp[,1])
    locs <- match(gene.labels, accs)
    all.gene.descs <- as.character(temp[locs, "Gene.Title"])
    all.gene.symbols <- as.character(temp[locs, "Gene.Symbol"])
    rm(temp)
  }
  
  if (is.data.frame(gs.ann)) {
    temp <- gs.ann
    a.size <- length(temp[,1])
    print(c("Number of gene set annotation file entries:", a.size))  
    accs <- as.character(temp[,1])
    locs <- match(gs.names, accs)
    all.gs.descs <- as.character(temp[locs, "SOURCE"])
    rm(temp)
  } else if (gs.ann == "") {
    for (i in 1:Ng) {
      all.gs.descs[i] <- gs.desc[i]
    }
  } else {
    temp <- read.delim(gs.ann, header=T, sep="\t", comment.char="", as.is=T)
    a.size <- length(temp[,1])
    print(c("Number of gene set annotation file entries:", a.size))  
    accs <- as.character(temp[,1])
    locs <- match(gs.names, accs)
    all.gs.descs <- as.character(temp[locs, "SOURCE"])
    rm(temp)
  }
  
  
  Obs.indicator <- matrix(nrow= Ng, ncol=N)
  Obs.RES <- matrix(nrow= Ng, ncol=N)
  
  Obs.ES <- vector(length = Ng, mode = "numeric")
  Obs.arg.ES <- vector(length = Ng, mode = "numeric")
  Obs.ES.norm <- vector(length = Ng, mode = "numeric")
  
  time2 <- proc.time()
  
  # GSEA methodology
  
  # Compute observed and random permutation gene rankings
  
  obs.s2n <- vector(length=N, mode="numeric")
  signal.strength <- vector(length=Ng, mode="numeric")
  tag.frac <- vector(length=Ng, mode="numeric")
  gene.frac <- vector(length=Ng, mode="numeric")
  coherence.ratio <- vector(length=Ng, mode="numeric")
  obs.phi.norm <- matrix(nrow = Ng, ncol = nperm)
  correl.matrix <- matrix(nrow = N, ncol = nperm)
  obs.correl.matrix <- matrix(nrow = N, ncol = nperm)
  order.matrix <- matrix(nrow = N, ncol = nperm)
  obs.order.matrix <- matrix(nrow = N, ncol = nperm)
  
  nperm.per.call <- 100
  n.groups <- nperm %/% nperm.per.call
  n.rem <- nperm %% nperm.per.call
  n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
  n.ends <- cumsum(n.perms)
  n.starts <- n.ends - n.perms + 1
  
  if (n.rem == 0) {
    n.tot <- n.groups
  } else {
    n.tot <- n.groups + 1
  }
  
  for (nk in 1:n.tot) {
    call.nperm <- n.perms[nk]
    
    print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", n.starts[nk], "--", n.ends[nk], sep=" "))
    
    O <- GSEA.GeneRanking(A, class.labels, gene.labels, call.nperm, permutation.type = perm.type, sigma.correction = "GeneCluster", fraction=fraction, replace=replace, reverse.sign = reverse.sign)
    gc()
    
    order.matrix[,n.starts[nk]:n.ends[nk]] <- O$order.matrix
    obs.order.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
    correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$s2n.matrix
    obs.correl.matrix[,n.starts[nk]:n.ends[nk]] <- O$obs.s2n.matrix
    rm(O)
  }
  
  ##median of correlation of each genes in the gene list
  obs.s2n <- apply(obs.correl.matrix, 1, median)  # using median to assign enrichment scores
  obs.index <- order(obs.s2n, decreasing=T)            
  obs.s2n   <- sort(obs.s2n, decreasing=T)            
  
  obs.gene.labels <- gene.labels[obs.index]       
  obs.gene.descs <- all.gene.descs[obs.index]       
  obs.gene.symbols <- all.gene.symbols[obs.index]       
  
  for (r in 1:nperm) {
    correl.matrix[, r] <- correl.matrix[order.matrix[,r], r]
  }
  for (r in 1:nperm) {
    obs.correl.matrix[, r] <- obs.correl.matrix[obs.order.matrix[,r], r]
  }
  
  
  ##gene list ranked!
  gene.list2 <- obs.index
  for (i in 1:Ng) {
    print(paste("Computing observed enrichment for gene set:", i, gs.names[i], sep=" ")) 
    if( gs.names[i]=='DOID:14330'){
      1+1                                                                                      
    }
    gene.set <- gs[i,gs[i,] != "null"]
    gene.set2 <- vector(length=length(gene.set), mode = "numeric")
    gene.set2 <- match(gene.set, gene.labels)
    GSEA.results <- GSEA.EnrichmentScore(gene.list=gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector = obs.s2n)
    
    Obs.ES[i] <- GSEA.results$ES
    Obs.arg.ES[i] <- GSEA.results$arg.ES
    Obs.RES[i,] <- GSEA.results$RES
    Obs.indicator[i,] <- GSEA.results$indicator
    if (Obs.ES[i] >= 0) {  # compute signal strength
      tag.frac[i] <- sum(Obs.indicator[i,1:Obs.arg.ES[i]])/size.G[i]
      gene.frac[i] <- Obs.arg.ES[i]/N
    } else {
      tag.frac[i] <- sum(Obs.indicator[i, Obs.arg.ES[i]:N])/size.G[i]
      gene.frac[i] <- (N - Obs.arg.ES[i] + 1)/N
    }
    signal.strength[i] <- tag.frac[i] * (1 - gene.frac[i]) * (N / (N - size.G[i]))
  }
  
  # Compute enrichment for random permutations 
  
  phi <- matrix(nrow = Ng, ncol = nperm)
  phi.norm <- matrix(nrow = Ng, ncol = nperm)
  obs.phi <- matrix(nrow = Ng, ncol = nperm)
  
  if (reshuffling.type == "sample.labels") { # reshuffling phenotype labels
    for (i in 1:Ng) {
      print(paste("Computing random permutations' enrichment for gene set:", i, gs.names[i], sep=" ")) 
      if( gs.names[i]=='DOID:3856'){
        1+1                                                                                      
      }
      gene.set <- gs[i,gs[i,] != "null"]
      gene.set2 <- vector(length=length(gene.set), mode = "numeric")
      gene.set2 <- match(gene.set, gene.labels)
      for (r in 1:nperm) {
        gene.list2 <- order.matrix[,r]
        if (use.fast.enrichment.routine == F) {
          GSEA.results <- GSEA.EnrichmentScore(gene.list=gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=correl.matrix[, r])   
        } else {
          GSEA.results <- GSEA.EnrichmentScore2(gene.list=gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=correl.matrix[, r])   
        }
        phi[i, r] <- GSEA.results$ES
      }
      if (fraction < 1.0) { # if resampling then compute ES for all observed rankings
        for (r in 1:nperm) {
          obs.gene.list2 <- obs.order.matrix[,r]
          if (use.fast.enrichment.routine == F) {
            GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r]) 
          } else {
            GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])
          }
          obs.phi[i, r] <- GSEA.results$ES
        }
      } else { # if no resampling then compute only one column (and fill the others with the same value)
        obs.gene.list2 <- obs.order.matrix[,1]
        if (use.fast.enrichment.routine == F) {
          GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r]) 
        } else {
          GSEA.results <- GSEA.EnrichmentScore2(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])
        }
        #this is the same as Obs.ES!
        obs.phi[i, 1] <- GSEA.results$ES
        for (r in 2:nperm) {
          obs.phi[i, r] <- obs.phi[i, 1]
        }
      }
      gc()
    }
    
  } else if (reshuffling.type == "gene.labels") { # reshuffling gene labels
    for (i in 1:Ng) {
      gene.set <- gs[i,gs[i,] != "null"]
      gene.set2 <- vector(length=length(gene.set), mode = "numeric")
      gene.set2 <- match(gene.set, gene.labels)
      for (r in 1:nperm) {
        reshuffled.gene.labels <- sample(1:rows)
        GSEA.results <- GSEA.EnrichmentScore(gene.list=reshuffled.gene.labels, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.s2n)   
        phi[i, r] <- GSEA.results$ES
      }
      if (fraction < 1.0) { # if resampling then compute ES for all observed rankings
        for (r in 1:nperm) {
          obs.gene.list2 <- obs.order.matrix[,r]
          GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
          obs.phi[i, r] <- GSEA.results$ES
        }
      } else { # if no resampling then compute only one column (and fill the others with the same value)
        obs.gene.list2 <- obs.order.matrix[,1]
        GSEA.results <- GSEA.EnrichmentScore(gene.list=obs.gene.list2, gene.set=gene.set2, weighted.score.type=weighted.score.type, correl.vector=obs.correl.matrix[, r])   
        obs.phi[i, 1] <- GSEA.results$ES
        for (r in 2:nperm) {
          obs.phi[i, r] <- obs.phi[i, 1]
        }
      }
      gc()
    }
  }
  
  # Compute 3 types of p-values
  
  # Find nominal p-values       
  
  print("Computing nominal p-values...")
  
  p.vals <- matrix(0, nrow = Ng, ncol = 2)
  
  if (OLD.GSEA == F) {
    for (i in 1:Ng) {
      if(gs.names[i] == 'DOID:3856'){
        1+1
      }
      pos.phi <- NULL
      neg.phi <- NULL
      for (j in 1:nperm) {
        if (phi[i, j] >= 0) {
          pos.phi <- c(pos.phi, phi[i, j]) 
        } else {
          neg.phi <- c(neg.phi, phi[i, j]) 
        }
      }
      ES.value <- Obs.ES[i]
      if (ES.value >= 0) {
        p.vals[i, 1] <- signif(sum(pos.phi >= ES.value)/length(pos.phi), digits=10)
      } else {
        p.vals[i, 1] <- signif(sum(neg.phi <= ES.value)/length(neg.phi), digits=10)
      }
    }
  } else {  # For OLD GSEA compute the p-val using positive and negative values in the same histogram
    for (i in 1:Ng) {
      if (Obs.ES[i] >= 0) {
        p.vals[i, 1] <-  sum(phi[i,] >= Obs.ES[i])/length(phi[i,])
        p.vals[i, 1] <-  signif(p.vals[i, 1], digits=5)
      } else {
        p.vals[i, 1] <-  sum(phi[i,] <= Obs.ES[i])/length(phi[i,])
        p.vals[i, 1] <-  signif(p.vals[i, 1], digits=5)
      }
    }
  }
  
  # Find effective size 
  
  erf <- function (x) 
  {
    2 * pnorm(sqrt(2) * x)
  }
  
  KS.mean <- function(N) { # KS mean as a function of set size N
    S <- 0
    for (k in -100:100) {
      if (k == 0) next
      S <- S + 4 * (-1)**(k + 1) * (0.25 * exp(-2 * k * k * N) - sqrt(2 * pi) *  erf(sqrt(2 * N) * k)/(16 * k * sqrt(N)))
    }
    return(abs(S))
  }
  
  # KS.mean.table <- vector(length=5000, mode="numeric")
  
  # for (i in 1:5000) {
  #    KS.mean.table[i] <- KS.mean(i)
  # }
  
  # KS.size <-  vector(length=Ng, mode="numeric")
  
  # Rescaling normalization for each gene set null
  
  print("Computing rescaling normalization for each gene set null...")
  
  if (OLD.GSEA == F) {
    for (i in 1:Ng) {
      if(gs.names[i] == 'DOID:9074'){
        1+1
      }
      pos.phi <- NULL
      neg.phi <- NULL
      for (j in 1:nperm) {
        if (phi[i, j] >= 0) {
          pos.phi <- c(pos.phi, phi[i, j]) 
        } else {
          neg.phi <- c(neg.phi, phi[i, j]) 
        }
      }
      pos.m <- mean(pos.phi)
      neg.m <- mean(abs(neg.phi))
      
      #         if (Obs.ES[i] >= 0) {
      #            KS.size[i] <- which.min(abs(KS.mean.table - pos.m))
      #         } else {
      #            KS.size[i] <- which.min(abs(KS.mean.table - neg.m))
      #         }
      
      #pos.phi <- pos.phi/pos.m
      #neg.phi <- neg.phi/neg.m
      for (j in 1:nperm) {
        if (phi[i, j] >= 0) {
          phi.norm[i, j] <- phi[i, j]/pos.m
        } else {
          phi.norm[i, j] <- phi[i, j]/neg.m
        }
      }
      for (j in 1:nperm) {
        if (obs.phi[i, j] >= 0) {
          obs.phi.norm[i, j] <- obs.phi[i, j]/pos.m
        } else {
          obs.phi.norm[i, j] <- obs.phi[i, j]/neg.m
        }
      }
      if (Obs.ES[i] >= 0) {
        Obs.ES.norm[i] <- Obs.ES[i]/pos.m
      } else {
        Obs.ES.norm[i] <- Obs.ES[i]/neg.m
      }
    }
  } else {  # For OLD GSEA does not normalize using empirical scaling 
    for (i in 1:Ng) {
      for (j in 1:nperm) {
        phi.norm[i, j] <- phi[i, j]/400
      }
      for (j in 1:nperm) {
        obs.phi.norm[i, j] <- obs.phi[i, j]/400
      }
      Obs.ES.norm[i] <- Obs.ES[i]/400
    }
  }
  
  # Save intermedite results
  
  if (save.intermediate.results == T) {
    
    filename <- paste(output.directory, doc.string, ".phi.txt", sep="", collapse="")
    write.table(phi, file = filename, quote=F, col.names= F, row.names=F, sep = "\t")
    
    filename <- paste(output.directory, doc.string, ".obs.phi.txt", sep="", collapse="")
    write.table(obs.phi, file = filename, quote=F, col.names= F, row.names=F, sep = "\t")
    
    filename <- paste(output.directory, doc.string, ".phi.norm.txt", sep="", collapse="")
    write.table(phi.norm, file = filename, quote=F, col.names= F, row.names=F, sep = "\t")
    
    filename <- paste(output.directory, doc.string, ".obs.phi.norm.txt", sep="", collapse="")
    write.table(obs.phi.norm, file = filename, quote=F, col.names= F, row.names=F, sep = "\t")
    
    filename <- paste(output.directory, doc.string, ".Obs.ES.txt", sep="", collapse="")
    write.table(Obs.ES, file = filename, quote=F, col.names= F, row.names=F, sep = "\t")
    
    filename <- paste(output.directory, doc.string, ".Obs.ES.norm.txt", sep="", collapse="")
    write.table(Obs.ES.norm, file = filename, quote=F, col.names= F, row.names=F, sep = "\t")
  }
  #-------------------------------------------------------------------------------------------------------------------  
  #-------------------------------------------------------------------------------------------------------------------  
  #-------------------------------------------------------------------------------------------------------------------  
  #-------------------------------------------------------------------------------------------------------------------  
  # Compute FWER p-vals
  
  print("Computing FWER p-values...")
  
  if (OLD.GSEA == F) {
    max.ES.vals.p <- NULL
    max.ES.vals.n <- NULL
    for (j in 1:nperm) {
      pos.phi <- NULL
      neg.phi <- NULL
      for (i in 1:Ng) {
        if (phi.norm[i, j] >= 0) {
          pos.phi <- c(pos.phi, phi.norm[i, j]) 
        } else {
          neg.phi <- c(neg.phi, phi.norm[i, j]) 
        }
      }
      if (length(pos.phi) > 0) {
        max.ES.vals.p <- c(max.ES.vals.p, max(pos.phi))
      }
      if (length(neg.phi) > 0) {
        max.ES.vals.n <- c(max.ES.vals.n, min(neg.phi))
      }
    }
    for (i in 1:Ng) {
      ES.value <- Obs.ES.norm[i]
      if (Obs.ES.norm[i] >= 0) {
        p.vals[i, 2] <- signif(sum(max.ES.vals.p >= ES.value)/length(max.ES.vals.p), digits=5)
      } else {
        p.vals[i, 2] <- signif(sum(max.ES.vals.n <= ES.value)/length(max.ES.vals.n), digits=5)
      }
    }
  } else {  # For OLD GSEA compute the FWER using positive and negative values in the same histogram
    max.ES.vals <- NULL
    for (j in 1:nperm) {
      max.NES <- max(phi.norm[,j])
      min.NES <- min(phi.norm[,j])
      if (max.NES > - min.NES) {
        max.val <- max.NES
      } else {
        max.val <- min.NES
      }
      max.ES.vals <- c(max.ES.vals, max.val)
    }
    for (i in 1:Ng) {
      if (Obs.ES.norm[i] >= 0) {
        p.vals[i, 2] <- sum(max.ES.vals >= Obs.ES.norm[i])/length(max.ES.vals)
      } else {
        p.vals[i, 2] <- sum(max.ES.vals <= Obs.ES.norm[i])/length(max.ES.vals)
      }
      p.vals[i, 2] <-  signif(p.vals[i, 2], digits=4)
    }
  }
  
  # Compute FDRs 
  
  print("Computing FDR q-values...")
  
  NES <- vector(length=Ng, mode="numeric")
  phi.norm.mean  <- vector(length=Ng, mode="numeric")
  obs.phi.norm.mean  <- vector(length=Ng, mode="numeric")
  phi.norm.median  <- vector(length=Ng, mode="numeric")
  obs.phi.norm.median  <- vector(length=Ng, mode="numeric")
  phi.norm.mean  <- vector(length=Ng, mode="numeric")
  obs.phi.mean  <- vector(length=Ng, mode="numeric")
  FDR.mean <- vector(length=Ng, mode="numeric")
  FDR.median <- vector(length=Ng, mode="numeric")
  phi.norm.median.d <- vector(length=Ng, mode="numeric")
  obs.phi.norm.median.d <- vector(length=Ng, mode="numeric")
  
  Obs.ES.index <- order(Obs.ES.norm, decreasing=T)
  Orig.index <- seq(1, Ng)
  Orig.index <- Orig.index[Obs.ES.index]
  Orig.index <- order(Orig.index, decreasing=F)
  Obs.ES.norm.sorted <- Obs.ES.norm[Obs.ES.index]
  gs.names.sorted <- gs.names[Obs.ES.index]
  
  for (k in 1:Ng) {
    NES[k] <- Obs.ES.norm.sorted[k]
    ES.value <- NES[k]
    count.col <- vector(length=nperm, mode="numeric")
    obs.count.col <- vector(length=nperm, mode="numeric")
    for (i in 1:nperm) {
      phi.vec <- phi.norm[,i]
      obs.phi.vec <- obs.phi.norm[,i]
      if (ES.value >= 0) {
        count.col.norm <- sum(phi.vec >= 0)
        obs.count.col.norm <- sum(obs.phi.vec >= 0)
        count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec >= ES.value)/count.col.norm, 0)
        obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec >= ES.value)/obs.count.col.norm, 0)
      } else {
        count.col.norm <- sum(phi.vec < 0)
        obs.count.col.norm <- sum(obs.phi.vec < 0)
        count.col[i] <- ifelse(count.col.norm > 0, sum(phi.vec <= ES.value)/count.col.norm, 0)
        obs.count.col[i] <- ifelse(obs.count.col.norm > 0, sum(obs.phi.vec <= ES.value)/obs.count.col.norm, 0)
      }
    }
    phi.norm.mean[k] <- mean(count.col)
    obs.phi.norm.mean[k] <- mean(obs.count.col)
    phi.norm.median[k] <- median(count.col)
    obs.phi.norm.median[k] <- median(obs.count.col)
    FDR.mean[k] <- ifelse(phi.norm.mean[k]/obs.phi.norm.mean[k] < 1, phi.norm.mean[k]/obs.phi.norm.mean[k], 1)
    FDR.median[k] <- ifelse(phi.norm.median[k]/obs.phi.norm.median[k] < 1, phi.norm.median[k]/obs.phi.norm.median[k], 1)
  }
  
  # adjust q-values
  
  if (adjust.FDR.q.val == T) {
    pos.nes <- length(NES[NES >= 0])
    min.FDR.mean <- FDR.mean[pos.nes]
    min.FDR.median <- FDR.median[pos.nes]
    for (k in seq(pos.nes - 1, 1, -1)) {
      if (FDR.mean[k] < min.FDR.mean) {
        min.FDR.mean <- FDR.mean[k]
      }
      if (min.FDR.mean < FDR.mean[k]) {
        FDR.mean[k] <- min.FDR.mean
      }
    }
    
    neg.nes <- pos.nes + 1
    min.FDR.mean <- FDR.mean[neg.nes]
    min.FDR.median <- FDR.median[neg.nes]
    for (k in seq(neg.nes + 1, Ng)) {
      if (FDR.mean[k] < min.FDR.mean) {
        min.FDR.mean <- FDR.mean[k]
      }
      if (min.FDR.mean < FDR.mean[k]) {
        FDR.mean[k] <- min.FDR.mean
      }
    }
  }
  
  obs.phi.norm.mean.sorted <- obs.phi.norm.mean[Orig.index]
  phi.norm.mean.sorted <- phi.norm.mean[Orig.index]
  FDR.mean.sorted <- FDR.mean[Orig.index]
  FDR.median.sorted <- FDR.median[Orig.index]
  
  #   Compute global statistic
  
  glob.p.vals <- vector(length=Ng, mode="numeric")
  NULL.pass <- vector(length=nperm, mode="numeric")
  OBS.pass <- vector(length=nperm, mode="numeric")
  
  for (k in 1:Ng) {
    NES[k] <- Obs.ES.norm.sorted[k]
    if (NES[k] >= 0) {
      for (i in 1:nperm) {
        NULL.pos <- sum(phi.norm[,i] >= 0)            
        NULL.pass[i] <- ifelse(NULL.pos > 0, sum(phi.norm[,i] >= NES[k])/NULL.pos, 0)
        OBS.pos <- sum(obs.phi.norm[,i] >= 0)
        OBS.pass[i] <- ifelse(OBS.pos > 0, sum(obs.phi.norm[,i] >= NES[k])/OBS.pos, 0)
      }
    } else {
      for (i in 1:nperm) {
        NULL.neg <- sum(phi.norm[,i] < 0)
        NULL.pass[i] <- ifelse(NULL.neg > 0, sum(phi.norm[,i] <= NES[k])/NULL.neg, 0)
        OBS.neg <- sum(obs.phi.norm[,i] < 0)
        OBS.pass[i] <- ifelse(OBS.neg > 0, sum(obs.phi.norm[,i] <= NES[k])/OBS.neg, 0)
      }
    }
    glob.p.vals[k] <- sum(NULL.pass >= mean(OBS.pass))/nperm
  }
  glob.p.vals.sorted <- glob.p.vals[Orig.index]
  
  # Produce results report
  
  print("Producing result tables and plots...")
  
  Obs.ES <- signif(Obs.ES, digits=5)
  Obs.ES.norm <- signif(Obs.ES.norm, digits=5)
  p.vals <- signif(p.vals, digits=4)
  signal.strength <- signif(signal.strength, digits=3)
  tag.frac <- signif(tag.frac, digits=3)
  gene.frac <- signif(gene.frac, digits=3)
  FDR.mean.sorted <- signif(FDR.mean.sorted, digits=5)
  FDR.median.sorted <-  signif(FDR.median.sorted, digits=5)
  glob.p.vals.sorted <- signif(glob.p.vals.sorted, digits=5)
  
  report <- data.frame(cbind(gs.names, size.G, all.gs.descs, Obs.ES, Obs.ES.norm, p.vals[,1], FDR.mean.sorted, p.vals[,2], tag.frac, gene.frac, signal.strength, FDR.median.sorted, glob.p.vals.sorted))
  names(report) <- c("GS", "SIZE", "SOURCE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "Tag %", "Gene %", "Signal", "FDR (median)", "glob.p.val")
  #       print(report)
  report2 <- report
  report.index2 <- order(Obs.ES.norm, decreasing=T)
  for (i in 1:Ng) {
    report2[i,] <- report[report.index2[i],]
  }   
  report3 <- report
  report.index3 <- order(Obs.ES.norm, decreasing=F)
  for (i in 1:Ng) {
    report3[i,] <- report[report.index3[i],]
  }   
  phen1.rows <- length(Obs.ES.norm[Obs.ES.norm >= 0])
  phen2.rows <- length(Obs.ES.norm[Obs.ES.norm < 0])
  report.phen1 <- report2[1:phen1.rows,]
  report.phen2 <- report3[1:phen2.rows,]
  
  if (output.directory != "")  {
    if (phen1.rows > 0) {
      filename <- paste(output.directory, doc.string, ".SUMMARY.RESULTS.REPORT.", phen1,".txt", sep="", collapse="")
      write.table(report.phen1, file = filename, quote=F, row.names=F, sep = "\t")
    }
    if (phen2.rows > 0) {
      filename <- paste(output.directory, doc.string, ".SUMMARY.RESULTS.REPORT.", phen2,".txt", sep="", collapse="")
      write.table(report.phen2, file = filename, quote=F, row.names=F, sep = "\t")
    }
  }
  
  # Global plots
  
  if (output.directory != "")  {
    if (non.interactive.run == F) {
      if (.Platform$OS.type == "windows") {
        glob.filename <- paste(output.directory, doc.string, ".global.plots", sep="", collapse="")
        windows(width = 10, height = 10)
      } else if (.Platform$OS.type == "unix") {
        glob.filename <- paste(output.directory, doc.string, ".global.plots.pdf", sep="", collapse="")
        pdf(file=glob.filename, height = 10, width = 10)
      }
    } else {
      if (.Platform$OS.type == "unix") {
        glob.filename <- paste(output.directory, doc.string, ".global.plots.pdf", sep="", collapse="")
        pdf(file=glob.filename, height = 10, width = 10)
      } else if (.Platform$OS.type == "windows") {
        glob.filename <- paste(output.directory, doc.string, ".global.plots.pdf", sep="", collapse="")
        pdf(file=glob.filename, height = 10, width = 10)
      }
    }
  }
  
  nf <- layout(matrix(c(1,2,3,4), 2, 2, byrow=T), c(1,1), c(1,1), TRUE)
  
  # plot S2N correlation profile
  
  location <- 1:N
  max.corr <- max(obs.s2n)
  min.corr <- min(obs.s2n)
  
  x <- plot(location, obs.s2n, ylab = "Signal to Noise Ratio (S2N)", xlab = "Gene List Location", main = "Gene List Correlation (S2N) Profile", type = "l", lwd = 2, cex = 0.9, col = 1)            
  for (i in seq(1, N, 20)) {
    lines(c(i, i), c(0, obs.s2n[i]), lwd = 3, cex = 0.9, col = colors()[12]) # shading of correlation plot
  }
  x <- points(location, obs.s2n, type = "l", lwd = 2, cex = 0.9, col = 1)            
  lines(c(1, N), c(0, 0), lwd = 2, lty = 1, cex = 0.9, col = 1) # zero correlation horizontal line
  temp <- order(abs(obs.s2n), decreasing=T)
  arg.correl <- temp[N]
  lines(c(arg.correl, arg.correl), c(min.corr, 0.7*max.corr), lwd = 2, lty = 3, cex = 0.9, col = 1) # zero correlation vertical line
  
  area.bias <- signif(100*(sum(obs.s2n[1:arg.correl]) + sum(obs.s2n[arg.correl:N]))/sum(abs(obs.s2n[1:N])), digits=3)
  area.phen <- ifelse(area.bias >= 0, phen1, phen2)
  delta.string <- paste("Corr. Area Bias to \"", area.phen, "\" =", abs(area.bias), "\\%", sep="", collapse="")
  zero.crossing.string <- paste("Zero Crossing at location ", arg.correl, " (",  signif(100*arg.correl/N, digits=3), " \\%)")
  leg.txt <- c(delta.string, zero.crossing.string)
  legend(x=N/10, y=max.corr, bty="n", bg = "white", legend=leg.txt, cex = 0.9)
  
  leg.txt <- paste("\"", phen1, "\" ", sep="", collapse="")
  text(x=1, y=-0.05*max.corr, adj = c(0, 1), labels=leg.txt, cex = 0.9)
  
  leg.txt <- paste("\"", phen2, "\" ", sep="", collapse="")
  text(x=N, y=0.05*max.corr, adj = c(1, 0), labels=leg.txt, cex = 0.9)
  
  if (Ng > 1) { # make these plots only if there are multiple gene sets.
    
    # compute plots of actual (weighted) null and observed
    
    phi.densities.pos <- matrix(0, nrow=512, ncol=nperm)
    phi.densities.neg <- matrix(0, nrow=512, ncol=nperm)
    obs.phi.densities.pos <- matrix(0, nrow=512, ncol=nperm)
    obs.phi.densities.neg <- matrix(0, nrow=512, ncol=nperm)
    phi.density.mean.pos <- vector(length=512, mode = "numeric")
    phi.density.mean.neg <- vector(length=512, mode = "numeric")
    obs.phi.density.mean.pos <- vector(length=512, mode = "numeric")
    obs.phi.density.mean.neg <- vector(length=512, mode = "numeric")
    phi.density.median.pos <- vector(length=512, mode = "numeric")
    phi.density.median.neg <- vector(length=512, mode = "numeric")
    obs.phi.density.median.pos <- vector(length=512, mode = "numeric")
    obs.phi.density.median.neg <- vector(length=512, mode = "numeric")
    x.coor.pos <-  vector(length=512, mode = "numeric")
    x.coor.neg <-  vector(length=512, mode = "numeric")
    
    for (i in 1:nperm) {
      pos.phi <- phi.norm[phi.norm[, i] >= 0, i]
      if (length(pos.phi) > 2) {
        temp <- density(pos.phi, adjust=adjust.param, n = 512, from=0, to=3.5)
      } else {
        temp <- list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      phi.densities.pos[, i] <- temp$y
      norm.factor <- sum(phi.densities.pos[, i])
      phi.densities.pos[, i] <- phi.densities.pos[, i]/norm.factor
      if (i == 1) {
        x.coor.pos <- temp$x
      }
      
      neg.phi <- phi.norm[phi.norm[, i] < 0, i]
      if (length(neg.phi) > 2) {
        temp <- density(neg.phi, adjust=adjust.param, n = 512, from=-3.5, to=0)
      } else {
        temp <- list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      phi.densities.neg[, i] <- temp$y
      norm.factor <- sum(phi.densities.neg[, i])
      phi.densities.neg[, i] <- phi.densities.neg[, i]/norm.factor
      if (i == 1) {
        x.coor.neg <- temp$x
      }
      pos.phi <- obs.phi.norm[obs.phi.norm[, i] >= 0, i]
      if (length(pos.phi) > 2) {
        temp <- density(pos.phi, adjust=adjust.param, n = 512, from=0, to=3.5)
      } else {
        temp <- list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      obs.phi.densities.pos[, i] <- temp$y
      norm.factor <- sum(obs.phi.densities.pos[, i])
      obs.phi.densities.pos[, i] <- obs.phi.densities.pos[, i]/norm.factor
      
      neg.phi <- obs.phi.norm[obs.phi.norm[, i] < 0, i]
      if (length(neg.phi)> 2) {  
        temp <- density(neg.phi, adjust=adjust.param, n = 512, from=-3.5, to=0)
      } else {
        temp <- list(x = 3.5*(seq(1, 512) - 1)/512, y = rep(0.001, 512))
      }
      obs.phi.densities.neg[, i] <- temp$y
      norm.factor <- sum(obs.phi.densities.neg[, i])
      obs.phi.densities.neg[, i] <- obs.phi.densities.neg[, i]/norm.factor
      
    }
    phi.density.mean.pos <- apply(phi.densities.pos, 1, mean)
    phi.density.mean.neg <- apply(phi.densities.neg, 1, mean)
    
    obs.phi.density.mean.pos <- apply(obs.phi.densities.pos, 1, mean)
    obs.phi.density.mean.neg <- apply(obs.phi.densities.neg, 1, mean)
    
    phi.density.median.pos <- apply(phi.densities.pos, 1, median)
    phi.density.median.neg <- apply(phi.densities.neg, 1, median)
    
    obs.phi.density.median.pos <- apply(obs.phi.densities.pos, 1, median)
    obs.phi.density.median.neg <- apply(obs.phi.densities.neg, 1, median)
    
    x <- c(x.coor.neg, x.coor.pos)
    x.plot.range <- range(x)
    y1 <- c(phi.density.mean.neg, phi.density.mean.pos)
    y2 <- c(obs.phi.density.mean.neg, obs.phi.density.mean.pos)
    y.plot.range <- c(-0.3*max(c(y1, y2)),  max(c(y1, y2)))
    
    print(c(y.plot.range, max(c(y1, y2)), max(y1), max(y2)))
    
    plot(x, y1, xlim = x.plot.range, ylim = 1.5*y.plot.range, type = "l", lwd = 2, col = 2, xlab = "NES", ylab = "P(NES)", main = "Global Observed and Null Densities (Area Normalized)")
    
    y1.point <- y1[seq(1, length(x), 2)]
    y2.point <- y2[seq(2, length(x), 2)]
    x1.point <- x[seq(1, length(x), 2)]
    x2.point <- x[seq(2, length(x), 2)]
    
    #     for (i in 1:length(x1.point)) {
    #       lines(c(x1.point[i], x1.point[i]), c(0, y1.point[i]), lwd = 3, cex = 0.9, col = colors()[555]) # shading 
    #     }
    #
    #     for (i in 1:length(x2.point)) {
    #       lines(c(x2.point[i], x2.point[i]), c(0, y2.point[i]), lwd = 3, cex = 0.9, col = colors()[29]) # shading 
    #     }
    
    points(x, y1, type = "l", lwd = 2, col = colors()[555])
    points(x, y2, type = "l", lwd = 2, col = colors()[29])
    
    for (i in 1:Ng) {
      col <- ifelse(Obs.ES.norm[i] > 0, 2, 3) 
      lines(c(Obs.ES.norm[i], Obs.ES.norm[i]), c(-0.2*max(c(y1, y2)), 0), lwd = 1, lty = 1, col = 1)
    }
    leg.txt <- paste("Neg. ES: \"", phen2, " \" ", sep="", collapse="")
    text(x=x.plot.range[1], y=-0.25*max(c(y1, y2)), adj = c(0, 1), labels=leg.txt, cex = 0.9)
    leg.txt <- paste(" Pos. ES: \"", phen1, "\" ", sep="", collapse="")
    text(x=x.plot.range[2], y=-0.25*max(c(y1, y2)), adj = c(1, 1), labels=leg.txt, cex = 0.9)
    
    leg.txt <- c("Null Density", "Observed Density", "Observed NES values")
    c.vec <- c(colors()[555], colors()[29], 1)
    lty.vec <- c(1, 1, 1)
    lwd.vec <- c(2, 2, 2)
    legend(x=0, y=1.5*y.plot.range[2], bty="n", bg = "white", legend=leg.txt, lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 0.9)
    
    B <- A[obs.index,]
    if (N > 300) {
      C <- rbind(B[1:100,], rep(0, Ns), rep(0, Ns), B[(floor(N/2) - 50 + 1):(floor(N/2) + 50),], rep(0, Ns), rep(0, Ns), B[(N - 100 + 1):N,])
    } 
    rm(B)
    GSEA.HeatMapPlot(V = C, col.labels = class.labels, col.classes = class.phen, main = "Heat Map for Genes in Dataset")
    
    # p-vals plot
    nom.p.vals <- p.vals[Obs.ES.index,1]
    FWER.p.vals <- p.vals[Obs.ES.index,2]
    plot.range <- 1.25*range(NES)
    plot(NES, FDR.mean, ylim = c(0, 1), xlim = plot.range, col = 1, bg = 1, type="p", pch = 22, cex = 0.75, xlab = "NES", main = "p-values vs. NES", ylab ="p-val/q-val")
    points(NES, nom.p.vals, type = "p", col = 2, bg = 2, pch = 22, cex = 0.75)
    points(NES, FWER.p.vals, type = "p", col = colors()[577], bg = colors()[577], pch = 22, cex = 0.75)
    leg.txt <- c("Nominal p-value", "FWER p-value", "FDR q-value")
    c.vec <- c(2, colors()[577], 1)
    pch.vec <- c(22, 22, 22)
    legend(x=-0.5, y=0.5, bty="n", bg = "white", legend=leg.txt, pch = pch.vec, col = c.vec, pt.bg = c.vec, cex = 0.9)
    lines(c(min(NES), max(NES)), c(nom.p.val.threshold, nom.p.val.threshold), lwd = 1, lty = 2, col = 2) 
    lines(c(min(NES), max(NES)), c(fwer.p.val.threshold, fwer.p.val.threshold), lwd = 1, lty = 2, col = colors()[577]) 
    lines(c(min(NES), max(NES)), c(fdr.q.val.threshold, fdr.q.val.threshold), lwd = 1, lty = 2, col = 1) 
    
    if (non.interactive.run == F) {  
      if (.Platform$OS.type == "windows") {
        savePlot(filename = glob.filename, type ="jpeg", device = dev.cur())
      } else if (.Platform$OS.type == "unix") {
        dev.off()
      }
    } else {
      dev.off()
    }
    
  } # if Ng > 1
  
  #----------------------------------------------------------------------------
  # Produce report for each gene set passing the nominal, FWER or FDR test or the top topgs in each side
  
  if (topgs > floor(Ng/2)) {
    topgs <- floor(Ng/2)
  }
  
  for (i in 1:Ng) {
    if ((p.vals[i, 1] <= nom.p.val.threshold) ||
          (p.vals[i, 2] <= fwer.p.val.threshold) ||
          (FDR.mean.sorted[i] <= fdr.q.val.threshold) || 
          (is.element(i, c(Obs.ES.index[1:topgs], Obs.ES.index[(Ng - topgs + 1): Ng])))) {
      
      #  produce report per gene set
      
      kk <- 1
      gene.number <- vector(length = size.G[i], mode = "character")
      gene.names <- vector(length = size.G[i], mode = "character")
      gene.symbols <- vector(length = size.G[i], mode = "character")
      gene.descs <- vector(length = size.G[i], mode = "character")
      gene.list.loc <- vector(length = size.G[i], mode = "numeric")
      core.enrichment <- vector(length = size.G[i], mode = "character")
      gene.s2n <- vector(length = size.G[i], mode = "numeric")
      gene.RES <- vector(length = size.G[i], mode = "numeric")
      rank.list <- seq(1, N)
      
      if (Obs.ES[i] >= 0) {
        set.k <- seq(1, N, 1)
        phen.tag <- phen1
        loc <- match(i, Obs.ES.index)
      } else {
        set.k <- seq(N, 1, -1)
        phen.tag <- phen2
        loc <- Ng - match(i, Obs.ES.index) + 1
      }
      
      for (k in set.k) {
        if (Obs.indicator[i, k] == 1) {
          gene.number[kk] <- kk
          gene.names[kk] <- obs.gene.labels[k]
          gene.symbols[kk] <- substr(obs.gene.symbols[k], 1, 15)
          gene.descs[kk] <- substr(obs.gene.descs[k], 1, 40)
          gene.list.loc[kk] <- k
          gene.s2n[kk] <- signif(obs.s2n[k], digits=3)
          gene.RES[kk] <- signif(Obs.RES[i, k], digits = 3)
          if (Obs.ES[i] >= 0) {
            core.enrichment[kk] <- ifelse(gene.list.loc[kk] <= Obs.arg.ES[i], "YES", "NO")
          } else {
            core.enrichment[kk] <- ifelse(gene.list.loc[kk] > Obs.arg.ES[i], "YES", "NO")
          }
          kk <- kk + 1
        }
      }
      
      gene.report <- data.frame(cbind(gene.number, gene.names, gene.symbols, gene.descs, gene.list.loc, gene.s2n, gene.RES, core.enrichment))
      names(gene.report) <- c("#", "GENE", "SYMBOL", "DESC", "LIST LOC", "S2N", "RES", "CORE_ENRICHMENT")
      
      #       print(gene.report)
      
      if (output.directory != "")  {
        
        filename <- paste(output.directory, doc.string, ".", gs.names[i], ".report.", phen.tag, ".", loc, ".txt", sep="", collapse="")
        write.table(gene.report, file = filename, quote=F, row.names=F, sep = "\t")
        
        if (non.interactive.run == F) {
          if (.Platform$OS.type == "windows") {
            gs.filename <- paste(output.directory, doc.string, ".", gs.names[i], ".plot.", phen.tag, ".", loc, sep="", collapse="")
            windows(width = 14, height = 6)
          } else if (.Platform$OS.type == "unix") {
            gs.filename <- paste(output.directory, doc.string, ".", gs.names[i], ".plot.", phen.tag, ".", loc, ".pdf", sep="", collapse="")
            pdf(file=gs.filename, height = 6, width = 14)
          }
        } else {
          if (.Platform$OS.type == "unix") {
            gs.filename <- paste(output.directory, doc.string, ".", gs.names[i], ".plot.", phen.tag, ".", loc, ".pdf", sep="", collapse="")
            pdf(file=gs.filename, height = 6, width = 14)
          } else if (.Platform$OS.type == "windows") {
            gs.filename <- paste(output.directory, doc.string, ".", gs.names[i], ".plot.", phen.tag, ".", loc, ".pdf", sep="", collapse="")
            pdf(file=gs.filename, height = 6, width = 14)
          }
        }
        
      }
      
      nf <- layout(matrix(c(1,2,3), 1, 3, byrow=T), 1, c(1, 1, 1), TRUE)
      ind <- 1:N
      min.RES <- min(Obs.RES[i,])
      max.RES <- max(Obs.RES[i,])
      if (max.RES < 0.3) max.RES <- 0.3
      if (min.RES > -0.3) min.RES <- -0.3
      delta <- (max.RES - min.RES)*0.50
      min.plot <- min.RES - 2*delta
      max.plot <- max.RES
      max.corr <- max(obs.s2n)
      min.corr <- min(obs.s2n)
      Obs.correl.vector.norm <- (obs.s2n - min.corr)/(max.corr - min.corr)*1.25*delta + min.plot
      zero.corr.line <- (- min.corr/(max.corr - min.corr))*1.25*delta + min.plot
      col <- ifelse(Obs.ES[i] > 0, 2, 4)
      
      # Running enrichment plot
      
      sub.string <- paste("Number of genes: ", N, " (in list), ", size.G[i], " (in gene set)", sep = "", collapse="")
      
      main.string <- paste("Gene Set ", i, ":", gs.names[i])
      plot(ind, Obs.RES[i,], main = main.string, sub = sub.string, xlab = "Gene List Index", ylab = "Running Enrichment Score (RES)", xlim=c(1, N), ylim=c(min.plot, max.plot), type = "l", lwd = 2, cex = 1, col = col)
      for (j in seq(1, N, 20)) {
        lines(c(j, j), c(zero.corr.line, Obs.correl.vector.norm[j]), lwd = 1, cex = 1, col = colors()[12]) # shading of correlation plot
      }
      lines(c(1, N), c(0, 0), lwd = 1, lty = 2, cex = 1, col = 1) # zero RES line
      lines(c(Obs.arg.ES[i], Obs.arg.ES[i]), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = col) # max enrichment vertical line
      for (j in 1:N) {
        if (Obs.indicator[i, j] == 1) {
          lines(c(j, j), c(min.plot + 1.25*delta, min.plot + 1.75*delta), lwd = 1, lty = 1, cex = 1, col = 1)  # enrichment tags
        }
      }
      lines(ind, Obs.correl.vector.norm, type = "l", lwd = 1, cex = 1, col = 1)
      lines(c(1, N), c(zero.corr.line, zero.corr.line), lwd = 1, lty = 1, cex = 1, col = 1) # zero correlation horizontal line
      temp <- order(abs(obs.s2n), decreasing=T)
      arg.correl <- temp[N]
      lines(c(arg.correl, arg.correl), c(min.plot, max.plot), lwd = 1, lty = 3, cex = 1, col = 3) # zero crossing correlation vertical line
      
      leg.txt <- paste("\"", phen1, "\" ", sep="", collapse="")
      text(x=1, y=min.plot, adj = c(0, 0), labels=leg.txt, cex = 1.0)
      
      leg.txt <- paste("\"", phen2, "\" ", sep="", collapse="")
      text(x=N, y=min.plot, adj = c(1, 0), labels=leg.txt, cex = 1.0)
      
      adjx <- ifelse(Obs.ES[i] > 0, 0, 1)
      
      leg.txt <- paste("Peak at ", Obs.arg.ES[i], sep="", collapse="")
      text(x=Obs.arg.ES[i], y=min.plot + 1.8*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
      
      leg.txt <- paste("Zero crossing at ", arg.correl, sep="", collapse="")
      text(x=arg.correl, y=min.plot + 1.95*delta, adj = c(adjx, 0), labels=leg.txt, cex = 1.0)
      
      # nominal p-val histogram
      
      sub.string <- paste("ES =", signif(Obs.ES[i], digits = 3), " NES =", signif(Obs.ES.norm[i], digits=3), "Nom. p-val=", signif(p.vals[i, 1], digits = 3),"FWER=", signif(p.vals[i, 2], digits = 3), "FDR=", signif(FDR.mean.sorted[i], digits = 3))
      temp <- density(phi[i,], adjust=adjust.param)
      x.plot.range <- range(temp$x)
      y.plot.range <- c(-0.125*max(temp$y), 1.5*max(temp$y))
      plot(temp$x, temp$y, type = "l", sub = sub.string, xlim = x.plot.range, ylim = y.plot.range, lwd = 2, col = 2, main = "Gene Set Null Distribution", xlab = "ES", ylab="P(ES)")
      x.loc <- which.min(abs(temp$x - Obs.ES[i]))
      lines(c(Obs.ES[i], Obs.ES[i]), c(0, temp$y[x.loc]), lwd = 2, lty = 1, cex = 1, col = 1)
      lines(x.plot.range, c(0, 0), lwd = 1, lty = 1, cex = 1, col = 1)
      
      leg.txt <- c("Gene Set Null Density", "Observed Gene Set ES value")
      c.vec <- c(2, 1)
      lty.vec <- c(1, 1)
      lwd.vec <- c(2, 2)
      legend(x=-0.2, y=y.plot.range[2], bty="n", bg = "white", legend=leg.txt, lty = lty.vec, lwd = lwd.vec, col = c.vec, cex = 1.0)
      
      leg.txt <- paste("Neg. ES \"", phen2, "\" ", sep="", collapse="")
      text(x=x.plot.range[1], y=-0.1*max(temp$y), adj = c(0, 0), labels=leg.txt, cex = 1.0)
      leg.txt <- paste(" Pos. ES: \"", phen1, "\" ", sep="", collapse="")
      text(x=x.plot.range[2], y=-0.1*max(temp$y), adj = c(1, 0), labels=leg.txt, cex = 1.0)
      
      # create pinkogram for each gene set
      
      kk <- 1
      
      pinko <- matrix(0, nrow = size.G[i], ncol = cols)
      pinko.gene.names <- vector(length = size.G[i], mode = "character")
      for (k in 1:rows) {
        if (Obs.indicator[i, k] == 1) {
          pinko[kk,] <- A[obs.index[k],]
          pinko.gene.names[kk] <- obs.gene.symbols[k]
          kk <- kk + 1
        }
      }
      GSEA.HeatMapPlot(V = pinko, row.names = pinko.gene.names, col.labels = class.labels, col.classes = class.phen, col.names = sample.names, main =" Heat Map for Genes in Gene Set", xlab=" ", ylab=" ")
      
      if (non.interactive.run == F) {  
        if (.Platform$OS.type == "windows") {
          savePlot(filename = gs.filename, type ="jpeg", device = dev.cur())
        } else if (.Platform$OS.type == "unix") {
          dev.off()
        }
      } else {
        dev.off()
      }
      
    } # if p.vals thres
    
  } # loop over gene sets
  
  
  return(list(report1 = report.phen1, report2 = report.phen2))
  
}  # end of definition of GSEA.analysis

