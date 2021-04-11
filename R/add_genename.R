#' add_genename
#'
#' sets up gsea
#' @param species species of the sample
#' @param mirror mirror for useEnsembl
#' @param data data object from get_deseq_res
#' @import biomaRt dplyr
#' @export
#'
#'
add_genename <- function(species, mirror = "useast", data){
  map_mart  <- useEnsembl(biomart="ensembl", dataset=species, mirror = mirror)
  data = left_join(data, getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                               filter = "ensembl_gene_id",
                               value = data$ensembl_gene_id,
                               mart = map_mart),
                   by = c("ensembl_gene_id" = "ensembl_gene_id"))
  return(data)
}
