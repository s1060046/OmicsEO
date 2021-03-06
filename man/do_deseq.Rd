% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_deseq.R
\name{do_deseq}
\alias{do_deseq}
\title{do_deseq}
\usage{
do_deseq(count_data, sample_data, design_formula, low_count_filter = 0)
}
\arguments{
\item{count_data}{countdata, from prep_for_deseq function}

\item{sample_data}{meta data}

\item{low_count_filter}{filter for low counts, will remove genes with average counts < filter. Defult = 0}

\item{desing_formula}{design formula for comparison}
}
\description{
makes DESeq2 object
}
