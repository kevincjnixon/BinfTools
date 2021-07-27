GSEAclus<-function(gsFile, term, clusName, dir){
  if(length(term)!=length(clusName)){
    stop("term must be of same length as clusName")
  }
  gs<-names(qusage::read.gmt(gsFile))
  clus<-rep("", length(gs))
  message("Identifying GSEA names with terms...")
  pb<-txtProgressBar(min=0, max=length(term), style=3)
  for(i in 1:length(term)){
    clus[grep(term[i], gs, ignore.case=T)]<-clusName[i]
    setTxtProgressBar(pb, i)
  }
  comb<-cbind(gs, clus)
  colnames(comb)<-c("NAME","Cluster")
  write.table(comb, paste0(dir,"GSEA_clusters.tsv"),
              quote=F, row.names=F, sep="\t")
  message("\nDone")
  #return(comb)
}

#' Create a rnk file for Gene Set Enrichment Analysis
#'
#' This function creates a ranked list of genes for use with a PreRanked
#' Gene Set Enrichment Analysis (GSEA) (www.gsea-msigdb.org/gsea)
#'
#' @param res A DESeq2 results object obtained from 'results(dds)' or a data.frame with the same column name values as a DESeq2 results object and rownames as genes
#' @param filename Path to the output .rnk file. Default is "./GSEA.rnk"
#' @param bystat Boolean values determining if genes should be ranked by column 'stat'. If TRUE, and no column 'stat', genes will be ranked using the -log10 of column 'pvalue'. Default is TRUE.
#' @param byFC Boolean values determining if genes should be ranked by column 'log2FoldChange'Default is FALSE.
#' @return Exports a ranked gene list. If both bystat and byFC are true, ranking will be abs(stat)*log2FoldChange.
#' @export

GenerateGSEA<-function(res, filename="GSEA.rnk", bystat=T, byFC=F){
  GSEA<-data.frame(NAME=rownames(res), Rank=rep(NA,nrow(res)))
  #res<-res[complete.cases(res),]
  #Check to see if there is a column named 'stat', if not, rename 'pvalue' to stat
  if(!any(colnames(res) %in% "stat")){
	print("No column named 'stat'... using -log10 'pvalue'")
	colnames(res)[which(colnames(res) %in% "pvalue")] <- "stat"
	res$stat <- -log(res$stat, 10)
  }
  if(isTRUE(bystat) && isFALSE(byFC)){
	  for(i in 1:nrow(res)){
		if(res$log2FoldChange[i] > 0){
		  GSEA$Rank[i]<-(abs(res$stat[i]))
		}
		if(res$log2FoldChange[i] < 0 ){
		  GSEA$Rank[i]<-(abs(res$stat[i]))*-1
		}
	  }
  }
  if(isTRUE(byFC) && isFALSE(bystat)){
	  for(i in 1:nrow(res)){
		if(res$log2FoldChange[i] > 0){
		  GSEA$Rank[i]<-res$log2FoldChange[i]
		}
		if(res$log2FoldChange[i] < 0 ){
		  GSEA$Rank[i]<-res$log2FoldChange[i]
		}
	  }
  }
  if(isTRUE(byFC) && isTRUE(bystat)){
	  for(i in 1:nrow(res)){
		if(res$log2FoldChange[i] > 0){
		  GSEA$Rank[i]<-abs(res$stat[i])*abs(res$log2FoldChange[i])
		}
		if(res$log2FoldChange[i] < 0 ){
		  GSEA$Rank[i]<-abs(res$stat[i])*abs(res$log2FoldChange[i])*-1
		}
	  }
  }
  write.table(GSEA[complete.cases(GSEA),], filename, quote=FALSE, row.names=FALSE, sep="\t")
}
