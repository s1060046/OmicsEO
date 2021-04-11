#' getGOterms
#'
#' sets up gsea
#' @param species species of the sample
#' @param mirror mirror for useEnsembl
#' @import biomaRt
#' @export
#'
#'
getGOterms <- function(species, mirror = "useast"){
  map_mart  <- useEnsembl(biomart="ensembl", dataset=species, mirror = mirror)
  GO_terms <- getBM(attributes = c("ensembl_gene_id", "name_1006", "external_gene_name"), mart = map_mart)
  return(GO_terms)
}
