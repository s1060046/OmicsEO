#' get_biomart
#'
#' grabs biomart info
#' @param species name of biomart species mouse = mmusculus_gene_ensembl, human = hsapiens_gene_ensembl
#' @import biomaRt
#' @export

get_biomart <- function(species){
  map_mart  <- useEnsembl(biomart="ensembl", dataset=species)

  exp_data <- getBM(attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                                   "external_gene_name","description",
                                   "5_utr_start", "5_utr_end",
                                   "transcript_length",
                                   "cds_length"),
                    mart = map_mart)

  return(exp_data)
}
