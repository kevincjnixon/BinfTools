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
                 pdf=TRUE, fig=TRUE, returnGost=FALSE, writeRes=TRUE, writeGem=FALSE, returnRes=FALSE, writeGene=FALSE){

  clusGenes<-lapply(split(clusRes, clusRes$cluster), rownames)
  names(clusGenes)<-paste0("cluster",names(clusGenes))

  return(BinfTools::GO_GEM(clusGenes, species=species, bg=bg, source=source, corr=corr,
                           iea=iea, prefix=prefix, ts=ts, pdf=pdf, fig=fig, returnRes=returnRes, returnGost=returnGost, writeGene=writeGene, writeGem=writeGem, writeRes=writeRes))
}
