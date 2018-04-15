require("biomaRt")

### here only provides human and mouse gene symbol converstion functions 

MouseGeneToHuman <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genelist, 
    mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genes)
}

PigGeneToHuman <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  pig = useMart("ensembl", dataset = "sscrofa_gene_ensembl")
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genelist, 
    mart = pig, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genes)
}

MonkeyGeneToHuman <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  monkey = useMart("ensembl", dataset = "mfascicularis_gene_ensembl")
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genelist, 
    mart = monkey, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genes)
}

HumanGeneToMouse <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genelist, 
    mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  return(genes)
}

HumanGeneToPig <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  pig = useMart("ensembl", dataset = "sscrofa_gene_ensembl")
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genelist, 
    mart = human, attributesL = c("hgnc_symbol"), martL = pig, uniqueRows=T)
  return(genes)
}

MouseGeneToHuman.ens <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genes = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = genelist, 
    mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
  return(genes)
}

HumanGeneToMouse.ens <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genes = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = genelist, 
    mart = human, attributesL = c("ensembl_gene_id"), martL = mouse, uniqueRows=T)
  return(genes)
}
