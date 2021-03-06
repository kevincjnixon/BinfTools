% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenerateGSEA.R
\name{GenerateGSEA}
\alias{GenerateGSEA}
\title{Create a rnk file for Gene Set Enrichment Analysis}
\usage{
GenerateGSEA(res, filename = "GSEA.rnk", bystat = T, byFC = F)
}
\arguments{
\item{res}{A DESeq2 results object obtained from 'results(dds)' or a data.frame with the same column name values as a DESeq2 results object and rownames as genes}

\item{filename}{Path to the output .rnk file. Default is "./GSEA.rnk"}

\item{bystat}{Boolean values determining if genes should be ranked by column 'stat'. If TRUE, and no column 'stat', genes will be ranked using the -log10 of column 'pvalue'. Default is TRUE.}

\item{byFC}{Boolean values determining if genes should be ranked by column 'log2FoldChange'Default is FALSE.}
}
\value{
Exports a ranked gene list. If both bystat and byFC are true, ranking will be abs(stat)*log2FoldChange.
}
\description{
This function creates a ranked list of genes for use with a PreRanked
Gene Set Enrichment Analysis (GSEA) (www.gsea-msigdb.org/gsea)
}
