#' deseq_plot_qc
#'
#' plots PCA
#' @param directory directory where raw output from feature counts are kept
#' @import limma
#' @import ggpubr
#' @import SummarizedExperiment
#' @export

deseq_plot_qc <- function(dds){
  vsd <- vst(dds, blind=FALSE)
  rld <- rlog(dds, blind=FALSE)
  mat_vsd <- assay(vsd)
  mat_vsd <- limma::removeBatchEffect(mat_vsd, vsd$pairs)
  assay(vsd) <- mat_vsd
  vsd_plot <- plotPCA(vsd, intgroup=c("condition")) + ggtitle("Variance stabilizing transformations")

  mat_rld <- assay(rld)
  mat_rld <- limma::removeBatchEffect(mat_rld, rld$pairs)
  assay(rld) <- mat_rld
  rld_plot <- plotPCA(rld, intgroup=c("condition")) + ggtitle("Regularized logarithm")

  ggpubr::ggarrange(vsd_plot, rld_plot)
}
