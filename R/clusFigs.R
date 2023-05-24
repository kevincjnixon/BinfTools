.clusHeatmap<-function(mat, gaps, title, annotdf, hmcol, labgenes){
  if(is.null(hmcol)){
    hmcol<-colorRampPalette(c("blue","grey","red"))(100)
  }
  mat<-as.matrix(mat)
  lim<-max(abs(mat[is.finite(mat)]))
  if(!is.null(labgenes)){
    tmp<-rep(" ", nrow(mat))
    for(i in 1:length(labgenes)){
      tmp[which(rownames(mat) %in% labgenes[i], arr.ind=TRUE)]<- labgenes[i]
    }
    labgenes<-tmp
  }
  pheatmap::pheatmap(mat, scale="none", color=hmcol, cluster_rows=FALSE, cluster_cols=FALSE,
                     legend=TRUE, labels_row = labgenes, annotation_row=annotdf, gaps_col=seq(1:ncol(mat-1)),
                     show_colnames = TRUE, gaps_row=gaps, main=title, breaks=seq(from=-lim, to=lim, length.out=100))
}
.clusBar<-function(mat, title, col, avgExp, cluslev=NULL){
  ylab<-"Average log2 Fold-Change (+/- sd)"
  if(isTRUE(avgExp)){
    ylab<-"Average Expression (+/- sd)"
  }
  clus<-mat$cluster
  if(!is(clus, "factor")){
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
  x <- x %>% dplyr::mutate(Comparison=forcats::fct_relevel(Comparison, as.character(colnames(mat))))
  x<-as.data.frame(x)
  #print(levels(x$cluster))
  #print(levels(x$Comparison))
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

#' Make a boxplot grouped by gene sets or clusters
#'
#' @param x Data frame or named list of data frames of average gene expression or log2FoldChanges of a subset of genes, where each column is a condition or contrast (in the order you want displayed), and rows are genes.
#' @param yax Character indicating the y-axis label of the boxplot. Default="Log2FoldChange".
#' @param col Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used for the bar plot. Default is "Dark2"
#' @param title Character indicating the title of the plot
#' @param showStat Boolean indicating if significance of pairwise contrasts within each group should be shown. Default=TRUE
#' @return Boxplot or grouped boxplot and tibble containing statistics from pairwise t-tests.
#' @export
#' @examples
#' groups<-gmtHeat(counts, cond, gmt, retGroups=T)
#' plotClusBox(x=groups, yax="Log2FoldChange", col="Dark2", title="", showStat=TRUE)

plotClusBox<-function (x, yax = "Log2FoldChange", col = "Dark2", title = "",
                       showStat = TRUE)
{
  require(dplyr, quietly = TRUE)
  res <- NULL
  compare <- NULL
  pwc <- NULL
  isList = FALSE
  if (is.list(x)) {
    isList = TRUE
    ord <- colnames(x[[1]])
    res <- lapply(x, tidyr::gather, key = "group", value = "Expression")
    res <- dplyr::bind_rows(res, .id = "cluster")
    res$group <- forcats::fct_relevel(as.factor(res$group),
                                      ord)
    res$group1 <- paste(res$group, res$cluster, sep = ".")
  }
  else {
    res <- tidyr::gather(x, key = "group", value = "Expression")
  }
  res <- as.data.frame(res)
  #return(res)
  if (isTRUE(isList)) {
    pwc <- res %>% dplyr::group_by(cluster) %>% rstatix::wilcox_test(Expression ~
                                                                       group) %>% rstatix::adjust_pvalue(method = "BH") %>%
      rstatix::add_significance("p.adj")
    #split pwc by cluster, then reorder and paste back into tible
    pwc<-split(pwc, f=pwc$cluster)
    pwc<-do.call("rbind",pwc[match(names(x), names(pwc))])
    pwc <- pwc %>% rstatix::add_xy_position(x = "cluster")
    pwc$x<-pwc$x[order(pwc$x)]
    pwc$xmin<-pwc$xmin[order(pwc$xmin)]
    pwc$xmax<-pwc$xmax[order(pwc$xmax)]

  }
  else {
    pwc <- res %>% rstatix::pairwise_wilcox_test(Expression ~
                                                   group, p.adjust.method = "BH")
    #split pwc by cluster, then reorder and paste back into tible
    pwc<-split(pwc, f=pwc$group)
    pwc<-do.call("rbind",pwc[match(colnames(x), names(pwc))])
    pwc <- pwc %>% rstatix::add_xy_position(x = "group")
    pwc$x<-pwc$x[order(pwc$x)]
    pwc$xmin<-pwc$xmin[order(pwc$xmin)]
    pwc$xmax<-pwc$xmax[order(pwc$xmax)]
  }
  p <- NULL
  if (isTRUE(isList)) {
    p <- ggpubr::ggboxplot(res, x = "cluster", y = "Expression",
                           fill = "group") + ggplot2::labs(title = title, y = yax,
                                                           x = "Cluster") + ggplot2::theme_minimal() + ggplot2::scale_fill_manual(values = BinfTools::colPal(col))
  }
  else {
    p <- ggpubr::ggboxplot(res, x = "group", y = "Expression",
                           fill = "group") + ggplot2::labs(title = title, y = yax,
                                                           x = "Condition") + ggplot2::theme_minimal() + ggplot2::scale_fill_manual(values = BinfTools::colPal(col))
  }
  if (isTRUE(showStat)) {
    if(isTRUE(isList)){
      p <- p + ggpubr::stat_pvalue_manual(pwc, label = "p.adj.signif",
                                          tip.length = 0, step.increase = 0.1, step.group.by = "cluster", hide.ns = TRUE)
    } else {
      p <- p + ggpubr::stat_pvalue_manual(pwc, label = "p.adj.signif",
                                          tip.length = 0, step.increase = 0.1, step.group.by = "group", hide.ns = TRUE)
    }
  }
  print(p)
  return(pwc)
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
#' @param showStat Boolean. If set to TRUE, significance symbols will be shown for significant contrasts in boxplot.
#' @param retStat Boolean. If set to TRUE, stats for contrasts (p-values and FDR) will be returned for each cluster along with results in a list object.
#' @return A data frame of log2FoldChanges for each comparison as columns (rows are genes) and a column named "cluster" indicating the cluster each gene belongs to. Two figures: a Heatmap of log2FoldChange of each gene ordered into clusters and a bar plot of the average log2FoldChange in each cluster (+/- SD) by comparison.
#' @export
clusFigs<-function(resList, numClus, title="Clustered Results", labgenes="", col="Dark2", hmcol=NULL, avgExp=FALSE, con=NULL, showStat=TRUE, retStat=FALSE){
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
    newClusLev<-order(unlist(lapply(split(mat[,which(colnames(mat) %in% con)], mat$cluster), mean)), decreasing=TRUE)
    mat <- mat %>% dplyr::mutate(cluster= forcats::fct_relevel(as.character(cluster), as.character(newClusLev)))
    newClus<-mat$cluster
    for(i in 1:numClus){
      message("changing cluster ",newClusLev[i]," to ",i,"...")
      newClus[which(mat$cluster == newClusLev[i])]<-i
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
  .clusHeatmap(mat[,-(which(colnames(mat) %in% "cluster"))], gaps, title, annotdf, hmcol, labgenes)
  #Plot the boxplot
  clusList<-split(mat[,-which(colnames(mat) %in% "cluster")], mat$cluster)
  yax="log2FoldChange"
  if(isTRUE(avgExp)){
    yax<-"Average Normalized Expression"
  }
  stats<-plotClusBox(clusList, yax, col, title, showStat)
  #Pass it through to the barplot function:
  tmp<-.clusBar(mat, title, col=col, avgExp)
  if(isTRUE(retStat)){
    return(list(clusRes=mat, stat=stats))
  }
  #return the cluster matrix
  return(mat)
}

#' Reorder clusters generated from clusFigs() and replot.
#'
#' @param clusRes Results dataframe from clusRes()
#' @param cluslev Numeric vector indicating the new order of clusters for clusRes
#' @param rename Boolean indicating if clusteres should be renamed based on their new order. Default=TRUE
#' @param title Character indicating the title for the plots. Default= "Releveled Clsuters"
#' @param col Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used for the bar plot. Default is "Dark2"
#' @param hmcol Colour Ramp Palette of length 100 indicating the colour palette of the heatmap. Leave NULL for default.
#' @param labgenes Character vector corresponding to rownames in clusRes of genes to label in heatmap. Set to NULL to label all genes. Default is "" which labeles no genes.
#' @param avgExp Boolean. Set to TRUE if clusRes was generated with 'avgExp=TRUE' in clusFigs().
#' @param showStat Boolean. If set to TRUE, significance symbols will be shown for significant contrasts in boxplot.
#' @param retStat Boolean. If set to TRUE, stats for contrasts (p-values and FDR) will be returned for each cluster along with results in a list object.
#' @return A data frame of log2FoldChanges for each comparison as columns (rows are genes) and a column named "cluster" indicating the cluster each gene belongs to. Two figures: a Heatmap of log2FoldChange of each gene ordered into clusters and a bar plot of the average log2FoldChange in each cluster (+/- SD) by comparison.
#' @export

clusRelev<-function(clusRes, cluslev, rename=TRUE, title="Releveled Clusters", col="Dark2", hmcol=NULL, labgenes="", avgExp=FALSE, showStat=TRUE, retStat=FALSE){
  numClus<-length(levels(factor(clusRes$cluster)))
  clusRes$cluster<-factor(clusRes$cluster)
  if(length(cluslev) == length(levels(factor(clusRes$cluster)))){
    clusRes <- clusRes %>% dplyr::mutate(cluster= forcats::fct_relevel(cluster, as.character(cluslev)))
  } else {
    newlev<-c(cluslev, levels(factor(clusRes$cluster))[!which(levels(factor(clusRes$cluster)) %in% cluslev)])
    clusRes <- clusRes %>% dplyr::mutate(cluster = forcats::fct_relevel(cluster, as.character(newlev)))
    cluslev<-as.numeric(newlev)
  }
  clusRes<-as.data.frame(clusRes)
  if(isTRUE(rename)){
    newClus<-clusRes$cluster
    for(i in 1:numClus){
      message("changing cluster ",cluslev[i]," to ",i,"...")
      newClus[which(clusRes$cluster == cluslev[i])]<-i
    }
    newClus<-forcats::fct_relevel(as.character(newClus), as.character(1:numClus))
    clusRes$cluster<-newClus
    clusRes<-clusRes[order(clusRes$cluster),]
    cluslev<-1:numClus
  }
  clusRes<-clusRes[order(clusRes$cluster),]
  ##Calculate the gaps for the heatmap
  gaps=c()
  for(i in 1:(length(cluslev)-1)){
    gaps<-c(gaps, sum(gaps[length(gaps)],length(which(clusRes$cluster == cluslev[i]))))
  }
  ##Create the annotation data frame
  annotdf<-data.frame(row.names=rownames(clusRes), cluster=clusRes$cluster)
  ##Pass it through to the heatmap function:
  .clusHeatmap(clusRes[,-(which(colnames(clusRes) %in% "cluster"))], gaps, title, annotdf, hmcol, labgenes)
  clusList<-split(clusRes[,-which(colnames(clusRes) %in% "cluster")], clusRes$cluster)
  yax="log2FoldChange"
  if(isTRUE(avgExp)){
    yax<-"Average Normalized Expression"
  }
  stats<-plotClusBox(clusList, yax, col, title, showStat)
  ##Pass it through to the barplot function:
  tmp<-.clusBar(clusRes, title, col=col, avgExp, cluslev=as.numeric(cluslev))
  if(isTRUE(retStat)){
    return(list(clusRes=clusRes, stat=stats))
  }
  ##return the cluster clusResrix
  return(clusRes)
}

#' Subset results objects using various parameters
#'
#' @param results data frame
#' @param genes character vector of genes in rownames(res) to subset. Default=NULL.
#' @param p adjusted p-value threshold. Genes must have padj<p to be returned. Default=NULL.
#' @param pval raw p-value threshold. Genes must have pvalue<pval to be returned. Default=NULL.
#' @param FC absolute log2FoldChange threshold. Genes must have |log2FoldChange|>FC to be returned. Default=NULL.
#' @return Data frame subset to genes matching any and all of the provided parametes.
#' @export

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
