#' Generate clusters from hierarchical heatmap clustering
#'
#' @param out output from zheat() when retClus=TRUE
#' @param level Level to cut the hierarchy tree to genreate clusters. Default is the rounded max tree height.
#' @return Data frame with rownames as genes and a column named "Cluster" indicating the clusters generated. And a figure showing the cut hierarchy tree.
#' @export
heatClus<-function(out, level=round(max(out$tree_row$height))){
  plot(out$tree_row)
  abline(h=level, col="red", lty=2, lwd=2)
  x<-sort(cutree(out$tree_row, h=level))
  return(data.frame(row.names=names(x),
                    Cluster=as.factor(x)))
}

#' Add clusters to results or count matrix
#'
#' Make a data frame of DESeq2 results or gene expression counts that shows the cluster each gene belongs to
#' @param resList A list or single data frame of DESeq2 results or gene expression counts
#' @param heatClus The data frame with rownames=genes and a column named Cluster. The output of heatClus
#' @return The original object(s) with a column named 'cluster' - for use with clusRelev()
#' @export

heatClusRes<-function(resList, heatClus){
  mat<-NULL
  if(is.data.frame(resList)){
    message("input is average expression count matrix.")
    mat<-resList
  } else {
    message("input is list of DESeq2 results objects.")
    if (length(resList) < 2) {
      stop("resList must have at least 2 entries...")
    }
    compNames <- names(resList)
    mat <- data.frame(genes = rownames(resList[[1]]), comp1 = resList[[1]]$log2FoldChange)
    for (i in 2:length(resList)) {
      tmp <- data.frame(genes = rownames(resList[[i]]),
        tmp = resList[[i]]$log2FoldChange)
      mat <- merge(mat, tmp, by.x = "genes", by.y = "genes")
    }
    rownames(mat) <- mat$genes
    mat <- mat[, -(which(colnames(mat) %in% "genes"))]
    colnames(mat) <- compNames
  }
  heatClus<-heatClus[which(rownames(heatClus) %in% rownames(mat)),, drop=FALSE]
  mat<-mat[which(rownames(mat) %in% rownames(heatClus)),]
  heatClus<-heatClus[match(rownames(mat), rownames(heatClus)),,drop=FALSE]
  mat$cluster<-heatClus$Cluster
  mat$cluster<-as.numeric(mat$cluster)
  mat <- mat[complete.cases(mat), ]
  mat<-mat[order(mat$cluster),]
  return(mat)
}


#' Z-score normalized heatmap of gene expression
#'
#' A function to generate a z-score normalized heatmap of gene expression from
#' a DESeq2 counts object.
#'
#' @param genes A character vector of genes to subset from rownames(counts) to create the heatmap
#' @param counts Normalized count matrix with rows as genes and columns as samples. Use 'counts(dds, normalized=TRUE)'.
#' @param conditions Character vector indicating conditions belonging to each sample (same order as colnames(counts))
#' @param con Character vector indicating the control condition level in 'conditions'. Default is "WT".
#' @param title Character vector indicating the title of the figure
#' @param labgenes Character vector corresponding to rownames(counts) of genes to be labelled on the side of the heatmap. Leave as NULL to label all genes. Use "" to label no genes.
#' @param zscore Boolean indicating if counts should be z-score normalized (default=TRUE)
#' @param avgExp Boolean indicating if values should be averaged within each condition (i.e. plot one value per condition). Default=FALSE.
#' @param rclus Boolean indicating if rows (genes) should be clustered. Leave FALSE for genes to be ordered on decreasing expression in 'con'.
#' @param hmcol Color ramp palette of length 100 for custom heatmap colouring. Leave NULL for defaults.
#' @param retClus Boolean indicating if clustered row names should be returnd. Default=FALSE.
#' @param annoCols Character vector indicating RColorBrewer palette name or colours for annotation colours if rclus is provided a data.frame.
#' @param fs Numeric indicating font size for row names. Default is 10.
#' @return Heatmap of z-score normalized gene expression with control condition samples appearing first.
#' @export

zheat<-function(genes=NULL, counts, conditions, con="WT", title="DEGs", labgenes=NULL, zscore=TRUE, avgExp=FALSE, rclus=FALSE, hmcol=NULL, retClus=FALSE, annoCols="Dark2", fs=10){
  if(is.null(hmcol)){
    hmcol<-colorRampPalette(c("blue","grey","red"))(100)
    if(isFALSE(zscore)){
      hmcol<-colorRampPalette(c("green","yellow","orange","red"))(100)
    }
  }
  zmat<-c()
  conditions<-as.factor(conditions)
  normcounts<-as.matrix(counts)
  ##If norm is true, zscore normalize the counts BEFORE pulling genes, otherwise, do it after
  if(isTRUE(zscore)){
    message("scaling to all genes...")
    normcounts<-t(scale(t(counts)))
  }
  ##genes points to the rows of counts to use in the heatmap - if null, use all genes in counts
	if(!is.null(genes)){
		message("pulling certain genes...")
		zmat<-normcounts[which(rownames(normcounts) %in% genes),]
	} else {
		zmat<-normcounts
	}
  ##Order the heatmap rows based on average control z-score (control condition is dicated by con
  ##Start by releveling the conditions so that con is first:
  if(length(con) == length(levels(factor(conditions)))){
    conditions<-forcats::fct_relevel(conditions, con)
  } else {
    ind<-which(levels(conditions) == con)
    if(ind != 1){
      x<-1:length(levels(conditions))
      x<-x[-ind]
      x<-c(ind,x)
      y<-levels(conditions)[x]
      conditions<-forcats::fct_relevel(conditions, y)
    }
  }
  zmat<-zmat[order(rowMeans(zmat[,which(conditions == con[1])]), decreasing=TRUE),]
  ##Now reorder the zmat to match the order of conditions (control then RNF)
  tmp.zmat<-zmat[,which(conditions == levels(conditions)[1])]
  tmp.conditions <- as.character(conditions[which(conditions == levels(conditions)[1])])
  for(i in 2:length(levels(conditions))){
    tmp.zmat<-cbind(tmp.zmat, zmat[,which(conditions == levels(conditions)[i])])
    tmp.conditions <- c(tmp.conditions, as.character(conditions[which(conditions == levels(conditions)[i])]))
  }
  zmat<-tmp.zmat
  conditions <- tmp.conditions
  if(isTRUE(avgExp)){
    message("Averaging values within each condition...")
    zmat<-as.matrix(avgExp(zmat, conditions, "mean"))
  }
  tmp=NULL
  ##Get the absolute max value to centre the scale around 0:
  lim<-c(0, max(zmat[is.finite(zmat)]))
  if(isTRUE(zscore)){
    lim<-c(max(abs(zmat[is.finite(zmat)]))*-1,max(abs(zmat[is.finite(zmat)])))
  }
  if(!is.logical(zscore) && zscore == "keepMin"){
    lim<-c(max(abs(zmat[is.finite(zmat)]))*-1,max(abs(zmat[is.finite(zmat)])))
  }
  ##Check to see if 'genes' is null and if it's not, replace it with a vector of length rownames
  ##Where genes specified in 'genes' are in the correct place, and all other values are " "
  if(!is.null(labgenes)){
    tmp<-rep(" ", nrow(zmat))
    for(i in 1:length(labgenes)){
      tmp[which(rownames(zmat) %in% labgenes[i], arr.ind=TRUE)]<- labgenes[i]
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
  ##Now make a heatmap
  out<-NULL
  if(is.logical(rclus)){
    out<-pheatmap::pheatmap(zmat, color=hmcol, show_colnames=TRUE, cluster_cols=FALSE, cluster_rows=rclus, main=title, labels_row=labgenes,breaks=seq(from=lim[1], to=lim[2], length.out=100), fontsize_row=fs)
  } else {
    rclus<-rclus[order(rclus[,1]),,drop=FALSE]
    zmat<-zmat[match(rownames(rclus), rownames(zmat)),]
    out<-pheatmap::pheatmap(zmat, color=hmcol, show_colnames=TRUE, cluster_cols=FALSE, cluster_rows=FALSE, main=title, annotation_row=rclus, annotation_colors=annoCols, labels_row=labgenes, breaks=seq(from=lim[1], to=lim[2], length.out=100), fontsize_row=10)
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
#' @param counts Normalized count matrix with rows as genes and columns as samples. Use 'counts(dds, normalized=TRUE)'.
#' @param conditions Character vector indicating conditions belonging to each sample (same order as colnames(counts))
#' @param order.by Character vector indicating the arbitrary appearance order of the conditions level. Default is the default order of levels(conditions).
#' @param sort.genes.by Character vector indicating the samples for sorting heatmap rows (genes).
#' @param title Character vector indicating the title of the figure
#' @param labgenes Character vector corresponding to rownames(counts) of genes to be labelled on the side of the heatmap. Leave as NULL to label all genes. Use "" to label no genes.
#' @param zscore Boolean indicating if counts should be z-score normalized (default=TRUE)
#' @param avgExp Boolean indicating if values should be averaged within each condition (i.e. plot one value per condition). Default=FALSE.
#' @param rclus Boolean indicating if rows (genes) should be clustered. Leave FALSE for genes to be ordered on decreasing expression in 'con'.
#' @param hmcol Color ramp palette of length 100 for custom heatmap colouring. Leave NULL for defaults.
#' @param retClus Boolean indicating if clustered row names should be returnd. Default=FALSE.
#' @param annoCols Character vector indicating RColorBrewer palette name or colours for annotation colours if rclus is provided a data.frame.
#' @param border_color color of cell borders on heatmap, use NA if no border should be drawn.
#' @return Heatmap of z-score normalized gene expression with control condition samples appearing first.
#' @export

zheat_v2 <- function (genes = NULL, counts, conditions, order.by = NULL, sort.genes.by = NULL, title = "DEGs",
                      labgenes = NULL, zscore = TRUE, avgExp = FALSE, rclus = FALSE, hmcol = NULL,
                      retClus = FALSE, annoCols = "Dark2", border_color = "grey60",
                      cellheight = NA, fontsize = 10)
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
    zmat <- zmat[order(rowMeans(zmat[, which(conditions == order.by[1])]), decreasing = TRUE), ]
  }else{
    zmat <- zmat[order(rowMeans(zmat[, which(conditions %in% sort.genes.by)]), decreasing = TRUE), ]
  }

  tmp.zmat <- zmat[, which(conditions == levels(conditions)[1])]
  tmp.conditions <- conditions[which(conditions == levels(conditions)[1])]
  ## Create a new vector
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
      tmp[which(rownames(zmat) %in% labgenes[i], arr.ind = TRUE)] <- labgenes[i]
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
    out <- pheatmap::pheatmap(zmat, color = hmcol, show_colnames = TRUE,
                              cluster_cols = FALSE, cluster_rows = rclus, main = title, border_color = border_color,
                              fontsize = fontsize, cellheight = cellheight,
                              labels_row = labgenes, breaks = seq(from = lim[1],
                                                                  to = lim[2], length.out = 100))
  }
  else {
    rclus <- rclus[order(rclus[, 1]), , drop = FALSE]
    zmat <- zmat[match(rownames(rclus), rownames(zmat)),
    ]
    out <- pheatmap::pheatmap(zmat, color = hmcol, show_colnames = TRUE,
                              cluster_cols = FALSE, cluster_rows = FALSE, main = title, border_color = border_color,
                              fontsize = fontsize, cellheight = cellheight,
                              annotation_row = rclus, annotation_colors = annoCols,
                              labels_row = labgenes, breaks = seq(from = lim[1],
                                                                  to = lim[2], length.out = 100))
  }
  if (isTRUE(retClus)) {
    return(out)
  }
}
