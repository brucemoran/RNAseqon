% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aracne-ap_viper.R
\name{highest_exp_wides}
\alias{highest_exp_wides}
\title{The situation arises where a single gene name maps to multiple gene identifiers
e.g. ENSG00000271503, ENSG00000274233 both map to CCL5
this function takes wide format input with at least 3 rows:
1: gene name
2: gene identifier
3: expression value(s) \link{can be index of columns}
function returns a wide format object with highest expression name+id
also returns full mapping info and mean expression value per name, id pairs}
\usage{
highest_exp_wides(wide_object, name_col, id_col, value_cols)
}
\arguments{
\item{wide_object}{a wide tibble or df}

\item{name_col}{naming column (gene name)}

\item{id_col}{identification column (gene ID, i.e. ENSID)}

\item{value_col}{expression value column (can be index)}
}
\value{
wide format with higehst exp. when multiple genes same named
}
\description{
The situation arises where a single gene name maps to multiple gene identifiers
e.g. ENSG00000271503, ENSG00000274233 both map to CCL5
this function takes wide format input with at least 3 rows:
1: gene name
2: gene identifier
3: expression value(s) \link{can be index of columns}
function returns a wide format object with highest expression name+id
also returns full mapping info and mean expression value per name, id pairs
}
