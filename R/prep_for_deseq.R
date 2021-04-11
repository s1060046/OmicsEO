#' prep_for_deseq
#'
#' prepares count data for DESeq
#' @param input needs to be an output from input_featurecounts
#' @param patter pattern to grab. defult is NULL and will return the whole input, if pattern is supplied it will return all columns with the pattern
#' @export
#'
#'
prep_for_deseq <- function(input, pattern=NULL){
  if(is.null(pattern)){
    row.names(input) <- input$gene.name
    input$gene.name <- NULL
    return(input)
  } else {
    gene_name <- input$gene.name

    data <- input[,grep(pattern, names(input))]

    row.names(data) <- gene_name
    return(data)
  }

}
