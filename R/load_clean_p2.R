#' load_clean_p2
#'
#' loads and cleans p2 data pre normalisation
#' @param input file from miguel
#' @param pep_filter filters peptide count less than given number
#' @import readr
#' @import DEP
#' @export
#'
#'
load_clean_p2 <- function(input, pep_filter){
  P2_data_noNorm <- read_csv(input)
  names(P2_data_noNorm) <- gsub("roll-up","LFQ", names(P2_data_noNorm))

  clean_DEP_p2 <- function(x){
    keep <- c("protGroups.input", "gene", "ID","gene.label", "pepcount", "name")
    info <- x[,keep]
    lfq <- x[,grep("LFQ", colnames(x))]
    exp <- cbind(info,lfq)
    names(exp) <- gsub(" ", ".", names(exp))
    return(exp)
  }
  P2_data_noNorm_pep <- subset(P2_data_noNorm, P2_data_noNorm$pepcount > pep_filter)
  P2_data_noNorm_pep <- make_unique(P2_data_noNorm_pep, "ID", "gene.label", delim = ";")
  P2_data_noNorm_pep_unique <- clean_DEP_p2(P2_data_noNorm_pep)

  return(P2_data_noNorm_pep_unique)
}
