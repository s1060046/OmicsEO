# OmicsEO
Collection of R scripts used for TRAP-Seq analysis

## Setup

Package can be downloaded by 
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
install.packages("devtools")
devtools::install_github(“s1060046/OmicsEO")
```

## Required files
1. need to have a folder containing all Featurecount outputs for DESeq2

## Example code for DESeq
1. Load in the package/
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
require(OmicsEO)
```

2. load in featurecount data to R environment.This code reads Featurecount files and compiles them to R object with gene.name column and each sample on each column.
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
counts <- input_featurecounts(directory = <your_directory>)
```
3. If you need to subset the count data you can use code below to subset your count data. Even if you don't use this to prepare the count data for DESeq2
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
counts_prep <- prep_for_deseq(<count_data>, <Unique identifier for the comparison>)
```

4.perform DESeq2 analysis. Make design dataframe as appropriate (Makae sure to name the control group "cont" which will be used as the reference point). For DESeq2 you can select low_count_filter. genes with average count less than the filter will be removed. For get_deseq_res select one of the lfcShrink options.
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
counts_prep <- data.frame(
  row.names = colnames(counts_prep),
  condition = c(rep("exp",3), rep("cont",3)), pairs = c("a","b","c","a","b","c")
)

veh_deseq <- do_deseq(count_data = counts_veh,
                                sample_data = counts_veh_design,
                                design_formula = ~ pairs + condition,
                                low_count_filter = 5)
veh_deseq_res <- get_deseq_res(veh_deseq, shrinkage = "normal")
```
