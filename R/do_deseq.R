#' do_deseq
#'
#' makes DESeq2 object
#' @param count_data countdata, from prep_for_deseq function
#' @param sample_data meta data
#' @param desing_formula design formula for comparison
#' @param low_count_filter filter for low counts, will remove genes with average counts < filter. Defult = 0
#' @import DESeq2
#' @export

do_deseq <- function(count_data, sample_data, design_formula,low_count_filter=0) {
  dds <- DESeqDataSetFromMatrix(countData =count_data, colData=sample_data, design=design_formula)

  count_filter_func <- function(x) { rowMeans(x) > low_count_filter}

  dds <- dds[count_filter_func(counts(dds)),]

  if("condition" %in% colnames(sample_data)) {
    dds$condition <- relevel(dds$condition, ref="cont")
  }

  results <- DESeq(dds)

  return(results)

}
