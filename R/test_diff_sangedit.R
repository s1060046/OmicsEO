#' test_diff_sangedit
#'
#' this is a modified version of test_diff function from DEP. It reports t statistics for gsea
#' @param se DEP object
#' @param type either control all manual
#' @param control if type == control set control
#' @param Character, The contrasts that will be tested if type = "manual". These should be formatted as "SampleA_vs_SampleB" or c("SampleA_vs_SampleC", "SampleB_vs_SampleC").
#' @param design_formula Formula, Used to create the design matrix.
#' @import limma
#' @import purrr
#' @import tibble
#' @import tidyr
#' @import dplyr
#' @import SummarizedExperiment
#' @export
#'
#'


test_diff_sangedit <- function (se, type = c("control", "all", "manual"), control = NULL,
                            test = NULL, design_formula = formula(~0 + condition))
{
  assertthat::assert_that(inherits(se, "SummarizedExperiment"),
                          is.character(type), class(design_formula) == "formula")

  type <- match.arg(type)
  col_data <- colData(se)
  raw <- assay(se)
  if (any(!c("name", "ID") %in% colnames(rowData(se, use.names = FALSE)))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_unique() and make_se() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(col_data))) {
    stop("'label', 'condition' and/or 'replicate' columns are not present in '",
         deparse(substitute(se)), "'\nRun make_se() or make_se_parse() to obtain the required columns",
         call. = FALSE)
  }
  if (any(is.na(raw))) {
    warning("Missing values in '", deparse(substitute(se)),
            "'")
  }
  if (!is.null(control)) {
    assertthat::assert_that(is.character(control), length(control) ==
                              1)
    if (!control %in% unique(col_data$condition)) {
      stop("run test_diff() with a valid control.\nValid controls are: '",
           paste0(unique(col_data$condition), collapse = "', '"),
           "'", call. = FALSE)
    }
  }
  variables <- terms.formula(design_formula) %>% attr(., "variables") %>%
    as.character() %>% .[-1]
  if (any(!variables %in% colnames(col_data))) {
    stop("run make_diff() with an appropriate 'design_formula'")
  }
  if (variables[1] != "condition") {
    stop("first factor of 'design_formula' should be 'condition'")
  }
  for (var in variables) {
    temp <- factor(col_data[[var]])
    assign(var, temp)
  }
  design <- model.matrix(design_formula, data = environment())
  colnames(design) <- gsub("condition", "", colnames(design))
  conditions <- as.character(unique(condition))
  if (type == "all") {
    cntrst <- apply(utils::combn(conditions, 2), 2, paste,
                    collapse = " - ")
    if (!is.null(control)) {
      flip <- grep(paste("^", control, sep = ""), cntrst)
      if (length(flip) >= 1) {
        cntrst[flip] <- cntrst[flip] %>% gsub(paste(control,
                                                    "- ", sep = " "), "", .) %>% paste(" - ",
                                                                                       control, sep = "")
      }
    }
  }
  if (type == "control") {
    if (is.null(control))
      stop("run test_diff(type = 'control') with a 'control' argument")
    cntrst <- paste(conditions[!conditions %in% control],
                    control, sep = " - ")
  }
  if (type == "manual") {
    if (is.null(test)) {
      stop("run test_diff(type = 'manual') with a 'test' argument")
    }
    assertthat::assert_that(is.character(test))
    if (any(!unlist(strsplit(test, "_vs_")) %in% conditions)) {
      stop("run test_diff() with valid contrasts in 'test'",
           ".\nValid contrasts should contain combinations of: '",
           paste0(conditions, collapse = "', '"), "', for example '",
           paste0(conditions[1], "_vs_", conditions[2]),
           "'.", call. = FALSE)
    }
    cntrst <- gsub("_vs_", " - ", test)
  }
  message("Tested contrasts: ", paste(gsub(" - ", "_vs_",
                                           cntrst), collapse = ", "))
  fit <- lmFit(raw, design = design)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- contrasts.fit(fit, made_contrasts)
  if (any(is.na(raw))) {
    for (i in cntrst) {
      covariates <- strsplit(i, " - ") %>% unlist
      single_contrast <- makeContrasts(contrasts = i,
                                       levels = design[, covariates])
      single_contrast_fit <- contrasts.fit(fit[, covariates],
                                           single_contrast)
      contrast_fit$coefficients[, i] <- single_contrast_fit$coefficients[,
                                                                         1]
      contrast_fit$stdev.unscaled[, i] <- single_contrast_fit$stdev.unscaled[,
                                                                             1]
    }
  }
  eB_fit <- eBayes(contrast_fit)
  retrieve_fun <- function(comp, fit = eB_fit) {
    res <- topTable(fit, sort.by = "t", coef = comp, number = Inf,
                    confint = TRUE)
    res <- res[!is.na(res$t), ]
    fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
    res$qval <- fdr_res$qval
    res$lfdr <- fdr_res$lfdr
    res$comparison <- rep(comp, dim(res)[1])
    res <- rownames_to_column(res)
    return(res)
  }
  limma_res <- map_df(cntrst, retrieve_fun)
  table <- limma_res %>% select(rowname, logFC, CI.L, CI.R, t,
                                P.Value, qval, comparison) %>% mutate(comparison = gsub(" - ",
                                                                                        "_vs_", comparison)) %>% gather(variable, value, -c(rowname,
                                                                                                                                            comparison)) %>% mutate(variable = recode(variable,
                                                                                                                                                                                      logFC = "log2FoldChange", P.Value = "pvalue", qval = "padj")) %>%
    unite(temp, comparison, variable) %>% spread(temp, value)
  rowData(se) <- merge(rowData(se, use.names = FALSE), table,
                       by.x = "name", by.y = "rowname", all.x = TRUE, sort = FALSE)
  return(se)
}
