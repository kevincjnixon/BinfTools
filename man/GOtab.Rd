% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GO_Tables.R
\name{GOtab}
\alias{GOtab}
\title{Pander table of top n GO results}
\usage{
GOtab(GOres, ts = c(10, 500), sig = TRUE, title = "", n = 10)
}
\arguments{
\item{GOres}{A single results table from GO_GEM() function with returnRes = T.}

\item{ts}{Numeric vector of length 2 indicating the minimum and maximum term sizes to filter GO results. Defaults to c(10,500).}

\item{sig}{Boolean indicating if the top significant terms should be shown. If FALSE, top enriched terms are shown. Default is TRUE.}

\item{title}{Character for the table legend description}

\item{n}{number of top terms to show. Default is 10. If there are fewer than n terms after filtering, only those terms will be shown.}
}
\value{
Prints a pandoc.table of results to stdout.
}
\description{
Make a Pandoc table of the top n results from a GO_GEM() results object. This is useful for making markdown reports.
}
