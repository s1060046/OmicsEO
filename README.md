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
1. Load in the package
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
require(OmicsEO)
```

2. load in featurecount data to R environment
this code reads Featurecount files and compiles them to R object with gene.name column and each sample on each column
```{r, echo = TRUE, eval = TRUE, collapse = TRUE}
counts <- input_featurecounts(directory = <your_directory>)
```
