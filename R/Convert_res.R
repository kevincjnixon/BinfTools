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
                     baseMean=topTable$AveExpr,
                  log2FoldChange=topTable$logFC,
                     stat=topTable$t,
                     pvalue=topTable$P.Value,
                     padj=topTable$adj.P.Val)
  return(res)
}

#'Convert outputs from tools other than DESeq2 for compatibility with BinfTools
#'
#'These functions convert results from Limma or EdgeR for compaitibility with BinfTools
#'
#'@param topTags A TopTags object originating from the output of edgeR::topTags()
#'@param geneNames A character vector specifying the names of the genes. Defaults to rownames(topTable) or rownames(topTags$table).
#'@param statCol A character indicating the column name of topTags$table that contains the calculated stats (e.g. 'LR', 'F', etc) - default NULL uses -log10(PValue)
#'@return A data frame of results compatible with BinfTools commands
#'@export

fromEdgeR<-function(topTags, geneNames=rownames(topTags$table), statCol=NULL){
  print(geneNames)
  res<-data.frame(row.names=geneNames,
                      baseMean=2^(topTags$table$logCPM),
                  log2FoldChange=topTags$table$logFC,
                     pvalue=topTags$table$PValue,
                     padj=topTags$table$FDR)
  if(is.null(statCol)){
    res$stat<- -log10(res$pvalue)
  } else {
    res$stat<- topTags$table[,which(colnames(topTags$table) %in% statCol)]
  }
  res<-res[,c(1,2,5,3,4)]
  return(res)
}
