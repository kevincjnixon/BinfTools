#'Convert outputs from tools other than DESeq2 for compatibility with BinfTools
#'
#'These functions convert results from Limma or EdgeR for compaitibility with BinfTools
#'
#'@param topTable A data frame originating from the output of limma::topTable()
#'@param geneNames Acharacter vector specifying the names of the genes. Defaults to rownames(topTable) or rownames(topTags$table).
#'@return A data frame of results compatible with BinfTools commands
#' @export

fromLimma<-function(topTable, geneNames=rownames(topTable)){
  res<-data.frame(row.names=geneNames,
                     baseMean=2^(topTable$AveExpr),
                  log2FoldChange=topTable$logFC,
                     pvalue=topTable$P.Value,
                     padj=topTable$adj.P.Val)
  return(res)
}

#'Convert outputs from tools other than DESeq2 for compatibility with BinfTools
#'
#'These functions convert results from Limma or EdgeR for compaitibility with BinfTools
#'
#'@param topTags A TopTags object originating from the output of edgeR::topTags()
#'@param geneNames Acharacter vector specifying the names of the genes. Defaults to rownames(topTable) or rownames(topTags$table).
#'@return A data frame of results compatible with BinfTools commands
#'@export

fromEdgeR<-function(topTags, geneNames=rownames(topTags$table)){
  print(geneNames)
  res<-data.frame(row.names=geneNames,
                      baseMean=2^(topTags$table$logCPM),
                  log2FoldChange=topTags$table$logFC,
                     pvalue=topTags$table$PValue,
                     padj=topTags$table$FDR)
  return(res)
}
