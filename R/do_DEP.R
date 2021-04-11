#' do_DEP
#'
#' perform DEP
#' @param dataframe proteomics dataset
#' @param metadata metadata for the dataframe
#' @param permutation Yes or No for performing 0 permutation
#' @param DE_type Differential expression type
#' @param control set control if DE type is control
#' @param test test conditions if DE type is manual
#' @param design_formula design formula
#' @import DEP
#' @export
#'
#'


do_DEP <- function (dataframe, metadata, permutation = c("Yes", "No"),
                    DE_type = c("control", "all", "manual"),
                    control = NULL,
                    test = NULL,
                    design_formula =formula(~0 + condition))
{
  se <- make_se(dataframe, grep("LFQ.", names(dataframe)), metadata)
  se <- normalize_vsn(se)
  if(permutation == "Yes"){
    proteins_MNAR <- get_df_long(se) %>%
      group_by(name, condition) %>%
      summarize(NAs = all(is.na(intensity))) %>%
      filter(NAs) %>%
      pull(name) %>%
      unique()

    # Get a logical vector
    MNAR <- names(se) %in% proteins_MNAR

    # Perform a mixed imputation
    se <- impute(
      se,
      fun = "mixed",
      randna = !MNAR, # we have to define MAR which is the opposite of MNAR
      mar = "bpca", # imputation function for MAR
      mnar = "QRILC") # imputation function for MNAR
  }


  data_diff <- test_diff_sangedit(se, type = DE_type,
                                  control = control, test = test, design_formula = design_formula )
  data_diff <- as.data.frame(rowData(data_diff))

  return(list(data_diff, se))

}
