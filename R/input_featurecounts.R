#' input_featurecounts
#'
#' takes directory where feature count results are kept and loads them. add _<sample name> to the file names
#' @param directory directory where raw output from feature counts are kept
#' @import readr dplyr
#' @export
#'
#'
#'
#'
input_featurecounts <- function(directory){
  #files need to be raw output from featurecounts, in rawdata folder and must contain sample ID after _ in file name
  file_name <-paste(directory , list.files(directory), sep="")
  for(i in file_name){
    if(i == file_name[1]){
      res <- read_table2(i, skip = 1)[c(1,7)]
      samp_name <- strsplit(i , "_")[[1]][length(strsplit(i , "_")[[1]])]
      names(res) <- c("gene.name", samp_name)
      raw_counts <- res
    } else{
      res <- read_table2(i, skip = 1)[c(1,7)]
      samp_name <- strsplit(i , "_")[[1]][length(strsplit(i , "_")[[1]])]
      names(res) <- c("gene.name", samp_name)
      raw_counts <- full_join(raw_counts, res)
    }
  }
  return(raw_counts)
}
