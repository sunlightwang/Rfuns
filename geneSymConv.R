require("biomaRt")

### here only provides human and mouse gene symbol converstion functions 

MouseGeneToHuman <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genelist, 
    mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}

HumanGeneToMouse <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genelist, 
    mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  return(genesV2)
}
