% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_plot.R
\name{count_plot}
\alias{count_plot}
\title{A function to make a plot of normalized counts}
\usage{
count_plot(
  counts,
  scaling = "zscore",
  genes,
  condition,
  con = NULL,
  title = "expression",
  compare = NULL,
  col = "Dark2",
  method = "ind",
  pair = F,
  pc = 1,
  yax = NULL,
  showStat = T,
  style = "violin"
)
}
\arguments{
\item{counts}{Normalized counts from a DESeq2 object - use 'counts(dds, normalized=T)'}

\item{scaling}{Method used to scale counts per gene across samples. 'zscore', 'log10', or 'none'. Default is 'zscore'}

\item{genes}{Character vector of genes to subset from counts. Must correspond with rownames(counts).}

\item{condition}{Character vector of conditions in DESeq2 object. Must be in order of columns (counts).}

\item{con}{Character indicating the control condition. (To be displayed first in plot). Default NULL.}

\item{title}{Character vector indicating title of plot. Defaults to "expression"}

\item{compare}{List of character vectors (each of length 2) indicating pairwise comparisons. If NULL, all possible comparisons will be made. Default is NULL}

\item{col}{Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used. Default is "Dark2"}

\item{method}{Character indicating what to plot. One of "ind", "mean", "geoMean", or "median", or "perMean". Defaults to "ind" for individual data points (one point per sample).}

\item{pair}{Boolean indicating if t-test should be independent (F; default) or paired (T).}

\item{pc}{Numeric indicating the pseudocount to be added when scaling="log10". Default=1.}

\item{yax}{Character indicating the y-axis label. Leave NULL if going with default axis label.}

\item{showStat}{Boolean indicating if statistics should be plotted.}

\item{style}{Character indicating the style of plot ("violin" or "box"). Defaults to "violin".}
}
\value{
Generates a violin or box plot
}
\description{
This function takes normalized counts of specific genes from a DESeq2 counts
object, scales them, and creates a plot with pairwise t-tests by condition
}
