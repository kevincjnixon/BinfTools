clusHeatmap<-function(mat, gaps, title, annotdf, hmcol, labgenes){
  if(is.null(hmcol)){
    hmcol<-colorRampPalette(c("blue","grey","red"))(100)
  }
  mat<-as.matrix(mat)
  lim<-max(abs(mat[is.finite(mat)]))
  if(!is.null(labgenes)){
    tmp<-rep(" ", nrow(mat))
    for(i in 1:length(labgenes)){
      tmp[which(rownames(mat) %in% labgenes[i], arr.ind=T)]<- labgenes[i]
    }
    labgenes<-tmp
  }
  pheatmap::pheatmap(mat, scale="none", color=hmcol, cluster_rows=F, cluster_cols=F,
                     legend=T, labels_row = labgenes, annotation_row=annotdf, gaps_col=seq(1:ncol(mat-1)),
                     show_colnames = T, gaps_row=gaps, main=title, breaks=seq(from=-lim, to=lim, length.out=100))
}
clusBar<-function(mat, title, col, avgExp, cluslev=NULL){
  ylab<-"Average log2 Fold-Change (+/- sd)"
  if(isTRUE(avgExp)){
    ylab<-"Average Expression (+/- sd)"
  }
  clus<-mat$cluster
  if(class(clus)!="factor"){
    clus<-factor(mat$cluster)
  }
  mat<-mat[,-(which(colnames(mat) %in% "cluster"))]
  clusList<-split(mat, clus)
  avgFC<-as.data.frame(do.call("rbind", lapply(clusList, colMeans)))
  colse<-function(x){
    y<-c()
    for(i in 1:ncol(x)){
      #y<-c(y, (sd(x[,i])/(sqrt(length(x[,i])))))
      y<-c(y, sd(x[,i]))
    }
    names(y)<-colnames(x)
    return(y)
  }
  seFC<-as.data.frame(do.call("rbind", lapply(clusList, colse)))
  clusters<-as.factor(c(rep(levels(clus), dim(avgFC)[2])))
  maxclus<-max(levels(clusters))
  if(is.null(cluslev)){
    cluslev<-c(1:maxclus)
  }
  y<- seFC %>% tidyr::gather(key="Comparison", value="SE")
  x<- avgFC %>% tidyr::gather(key="Comparison", value="AvgFC") %>% dplyr::mutate(cluster=clusters) %>%
    dplyr::mutate(SE=y$SE) %>% dplyr::group_by(cluster)
  x <- x %>% dplyr::mutate(cluster=forcats::fct_relevel(cluster, as.character(cluslev)))
  x<-as.data.frame(x)
  #print(levels(x$cluster))
  dodge<-ggplot2::position_dodge(width=0.9)
  p<-ggplot2::ggplot(x, ggplot2::aes(x=cluster, y=AvgFC, fill=factor(Comparison)))+
    ggplot2::geom_bar(stat="identity", position= ggplot2::position_dodge(), color="black") +
    ggplot2::geom_errorbar(ggplot2::aes(ymax=AvgFC+SE, ymin=AvgFC-SE), position=dodge,
                           width=0.2) +
    ggplot2::labs(title=title, y=ylab, x="Cluster") +
    ggplot2::guides(fill=ggplot2::guide_legend(title="Comparison")) +
    ggplot2::theme_minimal() + ggplot2::scale_fill_manual(values=colPal(col))
  print(p)
  return(x)
}
#' k-means clustered figures
#'
#' This function is for use with more than one results object. Combine the results
#' from multiple comparisons by clustering genes based on their log2 fold changes
#' between conditions for each comparison. This function will create a clustered
#' heatmap and bar plot of the average log2 fold-change in each cluster.
#' @param resList A list of results data frames. names(resList) will be used
#' @param numClus Number of k-means clusters to cluster the results
#' @param title Character indicating the titles of the plots to be made
#' @param labgenes Character vector corresponding to rownames in resList of genes to label in heatmap
#' @param col Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used for the bar plot. Default is "Dark2"
#' @param hmcol Colour Ramp Palette of length 100 indicating the colour palette of the heatmap. Leave NULL for default.
#' @param avgExp Boolean. If set to TRUE, resList should actually be an average expression matrix instead of list of DESeq2 results objects. (See AvgExp())
#' @param con Character indicating the control condition (equal to one of names(resList), or the colnames of the average expression matrix) to order the clusters from greatest to smallest average. Default is NULL for no ordering. Can also be releveled after using BinfTools:::clusRelev().
#' @return A data frame of log2FoldChanges for each comparison as columns (rows are genes) and a column named "cluster" indicating the cluster each gene belongs to. Two figures: a Heatmap of log2FoldChange of each gene ordered into clusters and a bar plot of the average log2FoldChange in each cluster (+/- SD) by comparison.
#' @export
clusFigs<-function(resList, numClus, title="Clustered Results", labgenes="", col="Dark2", hmcol=NULL, avgExp=F, con=NULL){
  mat<-NULL
  if(isTRUE(avgExp)){
    mat<-resList
  } else {
    #Make sure there is more than one entry in the list
    if(length(resList)<2){
      stop("resList must have at least 2 entries...")
    }
    #Get the names of the comparisons
    compNames<-names(resList)
    #Create a matrix for clustering based on the log2FoldChanges of each comparison
    mat<-data.frame(genes=rownames(resList[[1]]), comp1=resList[[1]]$log2FoldChange)
    for(i in 2:length(resList)){
      tmp<-data.frame(genes=rownames(resList[[i]]), tmp=resList[[i]]$log2FoldChange)
      mat<-merge(mat, tmp, by.x="genes", by.y="genes")
    }
    #Set genes to rownames
    rownames(mat)<-mat$genes
    #remove 'genes' column
    mat<-mat[,-(which(colnames(mat) %in% "genes"))]
    colnames(mat)<-compNames
    #Remove any NAs
    mat<-as.matrix(mat[complete.cases(mat),])
  }
  #Now we can cluster!
  print(paste("Clustering", length(resList),"results with",numClus,"clusters ..."))
  #Set seed
  set.seed(1234)
  clusRes<-kmeans(mat,centers=numClus, iter.max=1000)
  #Add the cluster to mat
  mat<-as.data.frame(mat)
  mat$cluster<-clusRes$cluster
  mat<-mat[order(mat$cluster),]
  if(!is.null(con)){
    newClusLev<-order(unlist(lapply(split(mat[,which(colnames(mat) %in% con)], mat$cluster), mean)), decreasing=T)
    mat <- mat %>% dplyr::mutate(cluster= forcats::fct_relevel(as.character(cluster), as.character(newClusLev)))
    newClus<-mat$cluster
    for(i in 1:numClus){
      message("changing cluster ",newClusLev[i]," to ",i,"...")
      newClus[which(mat$cluster==newClusLev[i])]<-i
    }
    newClus<-forcats::fct_relevel(as.character(newClus), as.character(1:numClus))
    mat$cluster<-newClus
    #mat <- mat%>% dplyr::mutate(cluster=forcats::fct_relevel(as.character(cluster), as.character(1:numClus)))
    mat<-mat[order(mat$cluster),]
  }
  #Calculate the gaps for the heatmap
  gaps=c()
  for(i in 1:(numClus-1)){
    gaps<-c(gaps, sum(gaps[length(gaps)],length(which(mat$cluster == i))))
  }
  #Create the annotation data frame
  annotdf<-data.frame(row.names=rownames(mat), cluster=mat$cluster)
  #Pass it through to the heatmap function:
  clusHeatmap(mat[,-(which(colnames(mat) %in% "cluster"))], gaps, title, annotdf, hmcol, labgenes)
  #Pass it through to the barplot function:
  tmp<-clusBar(mat, title, col=col, avgExp)
  #return the cluster matrix
  return(mat)
}

clusRelev<-function(clusRes, cluslev, rename=T, title="Releveled Clusters", col="Dark2", hmcol=NULL, labgenes=""){
  numClus<-max(clusRes$cluster)
  clusRes$cluster<-factor(clusRes$cluster)
  if(length(cluslev) == length(levels(factor(clusRes$cluster)))){
    clusRes <- clusRes %>% dplyr::mutate(cluster= forcats::fct_relevel(cluster, as.character(cluslev)))
  } else {
    newlev<-c(cluslev, levels(factor(clusRes$cluster))[!which(levels(factor(clusRes$cluster)) %in% cluslev)])
    clusRes <- clusRes %>% dplyr::mutate(cluster = forcats::fct_relevel(cluster, as.character(newlev)))
    cluslev<-as.numeric(newlev)
  }
  clusRes<-as.data.frame(clusRes)
  #clusRes$cluster<-as.numeric(clusRes$cluster)
  if(isTRUE(rename)){
    newClus<-clusRes$cluster
    for(i in 1:numClus){
      message("changing cluster ",clusLev[i]," to ",i,"...")
      newClus[which(clusRes$cluster==clusLev[i])]<-i
    }
    newClus<-forcats::fct_relevel(as.character(newClus), as.character(1:numClus))
    clusRes$cluster<-newClus
    #mat <- mat%>% dplyr::mutate(cluster=forcats::fct_relevel(as.character(cluster), as.character(1:numClus)))
    clusRes<-clusRes[order(clusRes$cluster),]
    
    #newClus<-clusRes$cluster
    #for(i in 1:numClus){
      ##message("Cluster ",as.numeric(cluslev[i])," is now cluster ",i,".")
      #newClus[clusRes$cluster==as.numeric(cluslev[i])]<-i
    #}
    #clusRes$cluster<-newClus
    #clusRes<-clusRes[order(clusRes$cluster),]
    #cluslev<-c(1:numClus)
  }
  #clusRes$cluster<-as.numeric(clusRes$cluster)
  clusRes<-clusRes[order(clusRes$cluster),]
  #Calculate the gaps for the heatmap
  gaps=c()
  for(i in 1:(numClus-1)){
    gaps<-c(gaps, sum(gaps[length(gaps)],length(which(clusRes$cluster == i))))
  }
  #Create the annotation data frame
  annotdf<-data.frame(row.names=rownames(clusRes), cluster=clusRes$cluster)
  #Pass it through to the heatmap function:
  clusHeatmap(clusRes[,-(which(colnames(clusRes) %in% "cluster"))], gaps, title, annotdf, hmcol, labgenes)
  #Pass it through to the barplot function:
  #print(as.numeric(cluslev))
  tmp<-clusBar(clusRes, title, col=col, cluslev=as.numeric(cluslev), avgExp=F)
  #return the cluster clusResrix
  return(clusRes)
}

#Subset results objects using various parameters

subRes<-function(res, genes=NULL, p=NULL, pval=NULL, FC=NULL){
  if(!is.null(genes)){
    message("subsetting using ", length(genes)," genes.")
    res<-res[which(rownames(res) %in% genes),]
  }
  if(!is.null(p)){
    message("subsetting using padj<",p,".")
    res<-subset(res, padj<p)
  }
  if(!is.null(pval)){
    message("subsetting using pvalue<",pval,".")
    res<-subset(res, pvalue<pval)
  }
  if(!is.null(FC)){
    message("subsetting using absolute log2FoldChange>", FC,".")
    res<-subset(res, abs(log2FoldChange)>FC)
  }
  return(res)
}
