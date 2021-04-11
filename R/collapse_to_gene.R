#' collapse_to_gene
#'
#' collapses transcript level info to gene level info my taking the most abundant transcript
#' @param tpm output from input_kallisto
#' @param trans_info output from get_biomart
#' @param selection pattern selection for TPM calculation. defult = NULL. set it to "WT" to only use WT data for
#' @export
#'
#'

collapse_to_gene <- function(tpm, trans_info, selection = NULL){
  if(is.null(selection)){
    tpm <- tpm
  } else {
    trans.name <-  tpm$trans.name
    dat <- tpm[,grep(selection, names(tpm))]
    tpm <- cbind(trans.name, dat)
  }
  detect_trans_info <- subset(trans_info, trans_info$ensembl_transcript_id %in% tpm$trans.name)
  message(paste("total gene",
              length(unique(trans_info$ensembl_gene_id)),
              ",",
              length(unique(detect_trans_info$ensembl_gene_id)),
              "detected", sep = " "))

  loop = 0
  Exp_info <- data.frame()
  for(gene in unique(detect_trans_info$ensembl_gene_id)){
    info_dat <- subset(detect_trans_info, detect_trans_info$ensembl_gene_id == gene)
    dat <- subset(tpm, tpm$trans.name %in% info_dat$ensembl_transcript_id)
    dat$av <- rowMeans(dat[-1])
    exp_data <- subset(dat, dat$av == max(dat$av))
    if (dim(exp_data)[1] > 1){
      exp_info <- subset(info_dat, info_dat$ensembl_transcript_id %in% exp_data$trans.name)
      exp_info$is.unique = "no"
      Exp_info <- rbind(Exp_info, exp_info)
      loop = loop +1
    } else{
      exp_info <- subset(info_dat, info_dat$ensembl_transcript_id == exp_data$trans.name)
      exp_info$is.unique = "yes"
      Exp_info <- rbind(Exp_info, exp_info)
      loop = loop +1
    }
    if(loop%%1000 == 0){message(paste("processed ", loop, " genes", sep = ""))}
  }
  return(Exp_info)
}
