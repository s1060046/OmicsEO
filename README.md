# OmicsEO
Collection of R scripts used for TRAP-Seq analysis. Run alignment and counting on Eddie before using this.

## Setup

Package can be downloaded by 
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
install.packages("devtools")
devtools::install_github(“s1060046/OmicsEO")
```

## Required files
1. need to have a folder containing all Featurecount outputs for DESeq2
2. need to have a separate folder containing all kallisto output for TPM

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

4.perform DESeq2 analysis. Make design dataframe as appropriate (Makae sure to name the control group "cont" which will be used as the reference point). For DESeq2 you can select low_count_filter. genes with average count less than the filter will be removed. For get_deseq_res select one of the lfcShrink options. lfc shrinkage is important for downstream GSEA analysis
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
counts_prep_design <- data.frame(
  row.names = colnames(counts_prep),
  condition = c(rep("exp",3), rep("cont",3)), pairs = c("a","b","c","a","b","c")
)

deseq <- do_deseq(count_data = counts_prep,
                                sample_data = counts_prep_design,
                                design_formula = ~ pairs + condition,
                                low_count_filter = 5)
deseq_res <- get_deseq_res(deseq, shrinkage = "normal")
```
5.For exporting data we are going to attach gene names and description. write csv function exports the data as csv file
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
deseq_res <- add_genename(deseq_res)
write.csv(deseq_res, "deseq_result.csv")
```
6.perform GSEA for GO level analysis. First, we need to generate a table with all the genes and their GO term information. getGOterm function takes species as input and generates this 
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
GO_term_list <- getGOterms("mmusculus_gene_ensembl")
gsea_res <- do_gsea(GO_terms = GO_term_list, data = deseq_res)
write.csv(gsea_res, "gsea_result.csv")
```

## Example code for TPM
1. load in featurecount data to R environment.This code reads Featurecount files and compiles them to R object with gene.name column and each sample on each column.
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
TPM <- input_kallisto(directory = <your_directory>)
```
2. generate table of transcript information and collapse TPMs to gene level
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
gene_transcript_info <- get_biomart("mmusculus_gene_ensembl")

TPM_gene <- collapse_to_gene(tpm = TPM, trans_info = gene_transcript_info)
```
