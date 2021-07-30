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
#' @param rclus Boolean indicating if rows (genes) should be clustered. Leave FALSE for genes to be ordered on decreasing expression in 'con'.
#' @param hmcol Color ramp palette of length 100 for custom heatmap colouring. Leave NULL for defaults.
#' @return Heatmap of z-score normalized gene expression with control condition samples appearing first.
#' @export

zheat<-function(genes=NULL, counts, conditions, con="WT", title="DEGs", labgenes=NULL, zscore=T, rclus=F, hmcol=NULL){
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
  tmp<-zmat[,which(conditions==levels(conditions)[1])]
  #print(head(tmp))
  for(i in 2:length(levels(conditions))){
    tmp<-cbind(tmp, zmat[,which(conditions==levels(conditions)[i])])
    #print(head(tmp))
  }
  zmat<-tmp
  #print(head(zmat))
  tmp=NULL
  #Get the absolute max value to centre the scale around 0:
  lim<-c(0, max(zmat[is.finite(zmat)]))
  if(isTRUE(zscore)){
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
  #Now make a heatmap
  if(is.logical(rclus)){
    pheatmap::pheatmap(zmat, color=hmcol, show_colnames=T, cluster_cols=F, cluster_rows=rclus, main=title, labels_row=labgenes,breaks=seq(from=lim[1], to=lim[2], length.out=100))
  } else {
    #rclus<-rclus[which(rownames(rclus) %in% rownames(zmat)),]
    rclus<-rclus[order(rclus[,1]),]
    zmat<-zmat[match(rownames(rclus), rownames(zmat)),]
    pheatmap::pheatmap(zmat, color=hmcol, show_colnames=T, cluster_cols=F, cluster_rows=F, main=title, annotation_row=rclus, labels_row=labgenes, breaks=seq(from=lim[1], to=lim[2], length.out=100))
  }
}
