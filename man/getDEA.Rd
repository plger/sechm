% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{getDEA}
\alias{getDEA}
\title{getDEA}
\usage{
getDEA(se, dea = NULL, homogenize = FALSE, sort = TRUE)
}
\arguments{
\item{se}{A \code{\link[SummarizedExperiment]{SummarizedExperiment-class}},
with DEAs each saved as a rowData column of `se`, with the column name
prefixed with "DEA."}

\item{dea}{The optional name of the DEA to extract}

\item{homogenize}{Logical; whether to homogenize the DEA}

\item{sort}{Logical; whether to return the table sorted by significance}
}
\value{
The DEA data.frame if `dea` is given, otherwise a named list of
  data.frames.
}
\description{
Extracts (standardized) DEA results from the rowData of an SE object.
}
\examples{
# loading example SE
data("Chen2017", package="sechm")
# this ones doesn't have saved DEAs in the standard format:
getDEA(Chen2017)
}
