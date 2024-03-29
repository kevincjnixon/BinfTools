% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zheat.R
\name{heatClus}
\alias{heatClus}
\title{Generate clusters from hierarchical heatmap clustering}
\usage{
heatClus(out, level = round(max(out$tree_row$height)))
}
\arguments{
\item{out}{output from zheat() when retClus=TRUE}

\item{level}{Level to cut the hierarchy tree to genreate clusters. Default is the rounded max tree height.}
}
\value{
Data frame with rownames as genes and a column named "Cluster" indicating the clusters generated. And a figure showing the cut hierarchy tree.
}
\description{
Generate clusters from hierarchical heatmap clustering
}
