require("biomaRt")

##### show ensembl data sets
# ensembl = useEnsembl(biomart="ensembl")
# listDatasets(ensembl)

### here only provides human and rat gene symbol converstion functions 

RatGeneToHuman <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  genes = getLDS(attributes = c("rgd_symbol"), filters = "mgi_symbol", values = genelist, 
    mart = rat, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genes)
}

HumanGeneToRat <- function(genelist){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  genes = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = genelist, 
    mart = human, attributesL = c("rgd_symbol"), martL = rat, uniqueRows=T)
  return(genes)
}
