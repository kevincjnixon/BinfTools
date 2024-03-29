% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Correlations.R
\name{avgExp}
\alias{avgExp}
\title{Calculate average gene expression between replicates}
\usage{
avgExp(counts, cond, method = c("mean", "median", "geoMean"))
}
\arguments{
\item{counts}{Data.frame of counts (rows=genes, columns=samples)}

\item{cond}{character vector of length ncol(counts) assigning each column to a condition}

\item{method}{character of either "mean", "median", or "geoMean" of how to calculate average expression. Default is "mean"}
}
\value{
Data frame with rows=genes and columns=condition with indicated average of gene expression beween replicates
}
\description{
Calculate average gene expression between replicates
}
