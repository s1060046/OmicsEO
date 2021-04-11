#' get_deseq_res
#'
#' exports DEseq2 results
#' @param dds deseq object
#' @param shrinkage shrinkage method
#' @export

get_deseq_res <- function(dds, shrinkage){
  res <- lfcShrink(dds = dds, type = shrinkage, coef = "condition_exp_vs_cont")
  res <- as.data.frame(res)
  res$ensembl_gene_id <- row.names(res)
  return(as.data.frame(res))
}
