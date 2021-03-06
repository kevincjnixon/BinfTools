#' Pathway Enrichment Analysis on Clusters
#'
#' Run Pathway Enrichment Analysis (same as GO_GEM()) on output from clusFigs()
#'
#' @param clusRes data.frame object that is output from clusFigs()
#' @param writeGene Boolean indicating if the genes in each cluster should be written out to '*prefix*_cluster*n*.genes.txt". Default is FALSE.
#' @param ... Arguments for GO_GEM() see ?GO_GEM for details. *bg* is automatically set as rownames(clusRes).
#' @return For each cluster, output from GO_GEM() is generated and output to "*prefix*_cluster_*n*"
#' @export

clusGO<-function(clusRes, species="hsapiens", bg=rownames(clusRes), source=NULL, corr="fdr", iea=FALSE, prefix="ClusGO", ts=c(10,500),
                 pdf=T, fig=T, returnGost=F, writeRes=T, writeGem=T, returnRes=F, writeGene=F){
  clusRes<-clusRes[order(clusRes$cluster),]
  clusters<-unique(clusRes$cluster)
  res<-list()
  for(i in 1:length(clusters)){
    genes<-rownames(subset(clusRes, cluster==clusters[i]))
    if(isTRUE(writeGene)){
      write.table(genes, file=paste0(prefix,"_cluster",clusters[i],".genes.txt"), quote=F, col.names=F, row.names=F, sep="\t")
    }
    print(paste("Analyzing",length(genes),"genes in cluster", clusters[i]))
    if(isFALSE(returnGost) & isFALSE(returnRes)){
      tryCatch({BinfTools::GO_GEM(geneList=genes, species=species, bg=bg, source=source, corr=corr,
                  iea=iea, prefix=paste0(prefix,"_cluster",clusters[i]), ts=ts, pdf=pdf, fig=fig, returnGost=returnGost,
                  writeRes=writeRes, writeGem=writeGem, returnRes=returnRes)},
               error=function(e){
                 message(paste0("No enriched terms in cluster ",clusters[i],"."))
                 message(paste0("Moving on to cluster ",clusters[i+1],"."))
               })
    }
    if(isTRUE(returnGost)){
      res[[i]]<-tryCatch({BinfTools::GO_GEM(geneList=genes, species=species, bg=bg, source=source, corr=corr,
                  iea=iea, prefix=paste0(prefix,"_cluster",clusters[i]), ts=ts, pdf=pdf, fig=fig, returnGost=T,
                  writeRes=writeRes, writeGem=writeGem, returnRes=returnRes)},
               error=function(e){
                 message(paste0("No enriched terms in cluster ",clusters[i],"."))
                 message(paste0("Moving on to cluster ",clusters[i+1],"."))
                 return(NA)
               })
      names(res)[i]<-paste0("cluster",clusters[i])
    }
    if(isTRUE(returnRes)){
      res[[i]]<-tryCatch({BinfTools::GO_GEM(geneList=genes, species=species, bg=bg, source=source, corr=corr,
                  iea=iea, prefix=paste0(prefix,"_cluster",clusters[i]), ts=ts, pdf=pdf, fig=fig, returnGost=F,
                  writeRes=writeRes, writeGem=writeGem, returnRes=T)},
               error=function(e){
                 message(paste0("No enriched terms in cluster ",clusters[i],"."))
                 message(paste0("Moving on to cluster ",clusters[i+1],"."))
                 return(NA)
               })
      names(res)[i]<-paste0("cluster",clusters[i])
    }
  }
  if(isTRUE(returnGost) | isTRUE(returnRes)){
    return(res)
  }
}
