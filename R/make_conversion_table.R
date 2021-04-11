#' make_conversion_table
#'
#' grabs biomart info
#' @param species name of biomart species mouse = mmusculus_gene_ensembl, human = hsapiens_gene_ensembl
#' @param species_1 name of biomart species mouse = mmusculus_gene_ensembl, human = hsapiens_gene_ensembl
#' @import biomaRt
#' @export

make_conversion_table <- function(species,species_1){
  map_mart  <- useEnsembl(biomart="ensembl", dataset=species)
  map_mart_1  <- useEnsembl(biomart="ensembl", dataset=species_1)

  exp_data <-getLDS(attributes = c("external_gene_name","ensembl_gene_id"),
                              mart = map_mart,
                              attributesL = c("external_gene_name","ensembl_gene_id"),
                              martL = map_mart_1)

  return(exp_data)
}
