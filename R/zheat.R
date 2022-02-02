heatClus<-function(out, level=round(max(out$tree_row$height))){
  plot(out$tree_row)
  abline(h=level, col="red", lty=2, lwd=2)
  x<-sort(cutree(out$tree_row, h=level))
  return(data.frame(row.names=names(x),
                    Cluster=as.factor(x)))
}
#' Z-score normalized heatmap of gene expression
#'
#' A function to generate a z-score normalized heatmap of gene expression from
#' a DESeq2 counts object.
#'
#' @param genes A character vector of genes to subset from rownames(counts) to create the heatmap
#' @param counts Normalized count matrix with rows as genes and columns as samples. Use 'counts(dds, normalized=T)'.
#' @param conditions Character vector indicating conditions belonging to each sample (same order as colnames(counts))
#' @param con Character vector indicating the control condition level in 'conditions'. Default is "WT".
#' @param title Character vector indicating the title of the figure
#' @param labgenes Character vector corresponding to rownames(counts) of genes to be labelled on the side of the heatmap. Leave as NULL to label all genes. Use "" to label no genes.
#' @param zscore Boolean indicating if counts should be z-score normalized (default=T)
#' @param avgExp Boolean indicating if values should be averaged within each condition (i.e. plot one value per condition). Default=F.
#' @param rclus Boolean indicating if rows (genes) should be clustered. Leave FALSE for genes to be ordered on decreasing expression in 'con'.
#' @param hmcol Color ramp palette of length 100 for custom heatmap colouring. Leave NULL for defaults.
#' @param retClus Boolean indicating if clustered row names should be returnd. Default=FALSE.
#' @param annoCols Character vector indicating RColorBrewer palette name or colours for annotation colours if rclus is provided a data.frame.
#' @return Heatmap of z-score normalized gene expression with control condition samples appearing first.
#' @export

zheat<-function(genes=NULL, counts, conditions, con="WT", title="DEGs", labgenes=NULL, zscore=T, avgExp=F, rclus=F, hmcol=NULL, retClus=F, annoCols="Dark2"){
  if(is.null(hmcol)){
    hmcol<-colorRampPalette(c("blue","grey","red"))(100)
    if(isFALSE(zscore)){
      hmcol<-colorRampPalette(c("green","yellow","orange","red"))(100)
    }
  }
  zmat<-c()
  conditions<-as.factor(conditions)
  normcounts<-as.matrix(counts)
  #If bnorm is true, zscore normalize the counts BEFORE pulling genes, otherwise, do it after
  if(isTRUE(zscore)){
    message("scaling to all genes...")
    #genes points to the rows of counts to use in the heatmap - if null, use all genes in counts
    normcounts<-t(scale(t(counts)))
  }
	if(!is.null(genes)){
		message("pulling certain genes...")
		zmat<-normcounts[which(rownames(normcounts) %in% genes),]
	} else {
		zmat<-normcounts
	}
  #Order the heatmap rows based on average control z-score (control condition is dicated by con
  #Start by releveling the conditions so that con is first:
  #print(levels(conditions))
  if(length(con)==length(levels(factor(conditions)))){
    conditions<-forcats::fct_relevel(conditions, con)
  } else {
    ind<-which(levels(conditions)==con)
    if(ind != 1){
      x<-1:length(levels(conditions))
      x<-x[-ind]
      x<-c(ind,x)
      y<-levels(conditions)[x]
      #levels(conditions)<-y
      conditions<-forcats::fct_relevel(conditions, y)
    }
  }
  #print(levels(conditions)[1])
  zmat<-zmat[order(rowMeans(zmat[,which(conditions==con[1])]), decreasing=T),]
  #Now reorder the zmat to match the order of conditions (control then RNF)
  tmp.zmat<-zmat[,which(conditions==levels(conditions)[1])]
  tmp.conditions <- as.character(conditions[which(conditions == levels(conditions)[1])])
  #print(head(tmp.zmat))
  for(i in 2:length(levels(conditions))){
    tmp.zmat<-cbind(tmp.zmat, zmat[,which(conditions==levels(conditions)[i])])
    tmp.conditions <- c(tmp.conditions, as.character(conditions[which(conditions == levels(conditions)[i])]))
    #print(head(tmp.zmat))
  }
  zmat<-tmp.zmat
  conditions <- tmp.conditions
  if(isTRUE(avgExp)){
    #print(head(zmat))
    #print(conditions)
    #print(levels(conditions))
    message("Averaging values within each condition...")
    zmat<-as.matrix(avgExp(zmat, conditions, "mean"))
    #zmat<-zmat[,levels(conditions)]
    #print(head(zmat))
  }
  #print(head(zmat))
  tmp=NULL
  #Get the absolute max value to centre the scale around 0:
  lim<-c(0, max(zmat[is.finite(zmat)]))
  if(isTRUE(zscore)){
    lim<-c(max(abs(zmat[is.finite(zmat)]))*-1,max(abs(zmat[is.finite(zmat)])))
  }
  if(!is.logical(zscore) && zscore=="keepMin"){
    lim<-c(max(abs(zmat[is.finite(zmat)]))*-1,max(abs(zmat[is.finite(zmat)])))
  }
  #print(lim)
  #Check to see if 'genes' is null and if it's not, replace it with a vector of length rownames
  #Where genes specified in 'genes' are in the correct place, and all other values are " "
  if(!is.null(labgenes)){
    tmp<-rep(" ", nrow(zmat))
    for(i in 1:length(labgenes)){
      tmp[which(rownames(zmat) %in% labgenes[i], arr.ind=T)]<- labgenes[i]
    }
    labgenes<-tmp
  }
  if(!is.data.frame(rclus)){
    annoCols=NULL
  } else {
    cols<-colPal(annoCols)
    annoCols<-list()
    for(i in 1:ncol(rclus)){
      rclus[,i]<-factor(rclus[,i])
      annoCols[[i]]<-cols[1:length(levels(rclus[,i]))]
      names(annoCols[[i]])<-levels(rclus[,i])
      names(annoCols)[i]<-colnames(rclus)[i]
    }
  }
  #Now make a heatmap
  out<-NULL
  if(is.logical(rclus)){
    out<-pheatmap::pheatmap(zmat, color=hmcol, show_colnames=T, cluster_cols=F, cluster_rows=rclus, main=title, labels_row=labgenes,breaks=seq(from=lim[1], to=lim[2], length.out=100))
  } else {
    #rclus<-rclus[which(rownames(rclus) %in% rownames(zmat)),]
    rclus<-rclus[order(rclus[,1]),,drop=F]
	  #print(head(zmat))
    zmat<-zmat[match(rownames(rclus), rownames(zmat)),]
	  #print(head(zmat))
    out<-pheatmap::pheatmap(zmat, color=hmcol, show_colnames=T, cluster_cols=F, cluster_rows=F, main=title, annotation_row=rclus, annotation_colors=annoCols, labels_row=labgenes, breaks=seq(from=lim[1], to=lim[2], length.out=100))
  }
  if(isTRUE(retClus)){
    return(out)
  }
}
######################################################################
######################################################################
#' Z-score normalized heatmap of gene expression
#'
#' A function to generate a z-score normalized heatmap of gene expression from
#' a DESeq2 counts object.
#'
#' @param genes A character vector of genes to subset from rownames(counts) to create the heatmap
#' @param counts Normalized count matrix with rows as genes and columns as samples. Use 'counts(dds, normalized=T)'.
#' @param conditions Character vector indicating conditions belonging to each sample (same order as colnames(counts))
#' @param order.by Character vector indicating the arbitrary appearance order of the conditions level. Default is the default order of levels(conditions).
#' @param sort.genes.by Character vector indicating the samples for sorting heatmap rows (genes).
#' @param title Character vector indicating the title of the figure
#' @param labgenes Character vector corresponding to rownames(counts) of genes to be labelled on the side of the heatmap. Leave as NULL to label all genes. Use "" to label no genes.
#' @param zscore Boolean indicating if counts should be z-score normalized (default=T)
#' @param avgExp Boolean indicating if values should be averaged within each condition (i.e. plot one value per condition). Default=F.
#' @param rclus Boolean indicating if rows (genes) should be clustered. Leave FALSE for genes to be ordered on decreasing expression in 'con'.
#' @param hmcol Color ramp palette of length 100 for custom heatmap colouring. Leave NULL for defaults.
#' @param retClus Boolean indicating if clustered row names should be returnd. Default=FALSE.
#' @param annoCols Character vector indicating RColorBrewer palette name or colours for annotation colours if rclus is provided a data.frame.
#' @return Heatmap of z-score normalized gene expression with control condition samples appearing first.
#' @export

zheat_v2 <- function (genes = NULL, counts, conditions, order.by = NULL, sort.genes.by = NULL, title = "DEGs",
                      labgenes = NULL, zscore = T, avgExp = F, rclus = F, hmcol = NULL,
                      retClus = F, annoCols = "Dark2")
{
  if (is.null(hmcol)) {
    hmcol <- colorRampPalette(c("blue", "grey", "red"))(100)
    if (isFALSE(zscore)) {
      hmcol <- colorRampPalette(c("green", "yellow", "orange",
                                  "red"))(100)
    }
  }
  zmat <- c()
  conditions <- as.factor(conditions)
  normcounts <- as.matrix(counts)
  if (isTRUE(zscore)) {
    message("scaling to all genes...")
    normcounts <- t(scale(t(counts)))
  }
  if (!is.null(genes)) {
    message("pulling certain genes...")
    zmat <- normcounts[which(rownames(normcounts) %in% genes),
    ]
  }
  else {
    zmat <- normcounts
  }

  if (is.null(order.by)){
    order.by <- levels(conditions)
  }else{
    if (length(order.by) != length(levels(factor(conditions)))){
      stop("The order.by vector must be a permutation of the levels of the conditions vector")
    }
    conditions <- forcats::fct_relevel(conditions, order.by)
  }

  if (is.null(sort.genes.by)){
    zmat <- zmat[order(rowMeans(zmat[, which(conditions == order.by[1])]), decreasing = T), ]
  }else{
    zmat <- zmat[order(rowMeans(zmat[, which(conditions %in% sort.genes.by)]), decreasing = T), ]
  }

  tmp.zmat <- zmat[, which(conditions == levels(conditions)[1])]
  tmp.conditions <- conditions[which(conditions == levels(conditions)[1])]
  # Creat a new vector
  for (i in 2:length(levels(conditions))) {
    tmp.zmat <- cbind(tmp.zmat, zmat[, which(conditions == levels(conditions)[i])])
    tmp.conditions <- c(tmp.conditions, conditions[which(conditions == levels(conditions)[i])])
  }
  zmat <- tmp.zmat
  conditions <- tmp.conditions
  if (isTRUE(avgExp)) {
    message("Averaging values within each condition...")
    zmat <- as.matrix(avgExp(zmat, conditions, "mean"))
    zmat <- zmat[, levels(conditions)]
  }
  tmp = NULL
  lim <- c(0, max(zmat[is.finite(zmat)]))
  if (isTRUE(zscore)) {
    lim <- c(max(abs(zmat[is.finite(zmat)])) * -1, max(abs(zmat[is.finite(zmat)])))
  }
  if (!is.logical(zscore) && zscore == "keepMin") {
    lim <- c(max(abs(zmat[is.finite(zmat)])) * -1, max(abs(zmat[is.finite(zmat)])))
  }
  if (!is.null(labgenes)) {
    tmp <- rep(" ", nrow(zmat))
    for (i in 1:length(labgenes)) {
      tmp[which(rownames(zmat) %in% labgenes[i], arr.ind = T)] <- labgenes[i]
    }
    labgenes <- tmp
  }
  if (!is.data.frame(rclus)) {
    annoCols = NULL
  }
  else {
    cols <- colPal(annoCols)
    annoCols <- list()
    for (i in 1:ncol(rclus)) {
      rclus[, i] <- factor(rclus[, i])
      annoCols[[i]] <- cols[1:length(levels(rclus[, i]))]
      names(annoCols[[i]]) <- levels(rclus[, i])
      names(annoCols)[i] <- colnames(rclus)[i]
    }
  }
  out <- NULL
  if (is.logical(rclus)) {
    out <- pheatmap::pheatmap(zmat, color = hmcol, show_colnames = T,
                              cluster_cols = F, cluster_rows = rclus, main = title,
                              labels_row = labgenes, breaks = seq(from = lim[1],
                                                                  to = lim[2], length.out = 100))
  }
  else {
    rclus <- rclus[order(rclus[, 1]), , drop = F]
    zmat <- zmat[match(rownames(rclus), rownames(zmat)),
    ]
    out <- pheatmap::pheatmap(zmat, color = hmcol, show_colnames = T,
                              cluster_cols = F, cluster_rows = F, main = title,
                              annotation_row = rclus, annotation_colors = annoCols,
                              labels_row = labgenes, breaks = seq(from = lim[1],
                                                                  to = lim[2], length.out = 100))
  }
  if (isTRUE(retClus)) {
    return(out)
  }
}
