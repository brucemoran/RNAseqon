% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_RNAseq_output.R
\name{nf_core_rnaseq_star_salmon_parser}
\alias{nf_core_rnaseq_star_salmon_parser}
\title{Parse featurecounts input as a follow-on from standard nf-core/rnaseq pipeline using STAR-Salmon (v3.3 currently)}
\usage{
nf_core_rnaseq_star_salmon_parser(tag = NULL)
}
\arguments{
\item{tag}{string to prefix output}
}
\value{
nz_fc_co raw count object with no lines summing to zeros
}
\description{
Parse featurecounts input as a follow-on from standard nf-core/rnaseq pipeline using STAR-Salmon (v3.3 currently)
}
