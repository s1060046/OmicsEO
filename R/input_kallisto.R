#' input_kallisto
#'
#' takes directory where feature count results are kept and loads them. add _<sample name> to the file names
#' @param directory directory where tsv output from kallisto is kept
#' @import readr
#' @import dplyr
#' @export

input_kallisto <- function(directory){
  #files need to be raw output from featurecounts, in rawdata folder and must contain sample ID after _ in file name
  file_name <-paste(directory , list.files(directory), sep="")
  for(i in file_name){
    if(i == file_name[1]){
      res <- read_table2(i)[c(1,5)]
      samp_name <- strsplit(i , "_")[[1]][length(strsplit(i , "_")[[1]])]
      names(res) <- c("trans.name", samp_name)
      raw_counts <- res
    } else{
      res <- read_table2(i)[c(1,5)]
      samp_name <- strsplit(i , "_")[[1]][length(strsplit(i , "_")[[1]])]
      names(res) <- c("trans.name", samp_name)
      raw_counts <- full_join(raw_counts, res)
    }
  }
  raw_counts$trans.name <- apply(raw_counts, 1, FUN = function(x){strsplit(x[1], "\\.")[[1]][1]})

  return(raw_counts)
}
