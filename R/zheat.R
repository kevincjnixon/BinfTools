#' Z-score normalized heatmap of gene expression
#'
#' A function to generate a z-score normalized heatmap of gene expression from
#' a DESeq2 counts object.
#'
#' @param genes A character vector of genes to subset from rownames(counts) to create the heatmap
#' @param bnorm Boolean values indicating if z-score normalization should occur before subsetting genes. Leave "TRUE"
#' @param counts Normalized count matrix with rows as genes and columns as samples. Use 'counts(dds, normalized=T)'.
#' @param conditions Factor vector indicating conditions belonging to each sample (same order as colnames(counts))
#' @param con Character vector indicating the control condition level in 'conditions'. Default is "WT".
#' @param title Character vector indicating the title of the figure
#' @param labgenes Character vector corresponding to rownames(counts) of genes to be labelled on the side of the heatmap. Leave as NULL to label all genes. Use "" to label no genes.
#' @return Heatmap of z-score normalized gene expression with control condition samples appearing first.
#' @export

zheat<-function(genes=NULL,bnorm=T, counts, conditions, con="WT", title="DEGs", labgenes=NULL){
  hmcol<-colorRampPalette(c("blue","grey","red"))(100)
  zmat<-c()
  conditions<-as.factor(conditions)
  #If bnorm is true, zscore normalize the counts BEFORE pulling genes, otherwise, do it after
  if(!isFALSE(bnorm)){
	print("scaling to all genes...")
	#genes points to the rows of counts to use in the heatmap - if null, use all genes in counts
	normcounts<-t(scale(t(counts)))
	if(!is.null(genes)){
		print("pulling certain genes...")
		zmat<-normcounts[which(rownames(normcounts) %in% genes),]
	} else {
		zmat<-normcounts
	}
  } else {
	normcounts<-counts
	if(!is.null(genes)){
		print("pulling certain genes...")
		zmat<-normcounts[which(rownames(normcounts) %in% genes),]
	} else {
		zmat<-normcounts
	}
	#Now get the z-score for rows (genes) using the scale function
	print("scaling to select genes...")
	zmat<-t(scale(t(zmat)))
  }
  #Order the heatmap rows based on average control z-score (control condition is dicated by con
  #Start by releveling the conditions so that con is first:
  print(levels(conditions))
  ind<-which(levels(conditions)==con)
  if(ind != 1){
	x<-1:length(levels(conditions))
	x<-x[-ind]
	x<-c(ind,x)
	y<-levels(conditions)[x]
	levels(conditions)<-y
  }
  print(levels(conditions))
  zmat<-zmat[order(rowMeans(zmat[,which(conditions==con)]), decreasing=T),]
  #Now reorder the zmat to match the order of conditions (control then RNF)
  tmp<-zmat[,which(conditions==levels(conditions)[1])]
  for(i in 2:length(levels(conditions))){
	tmp<-cbind(tmp, zmat[,which(conditions==levels(conditions)[i])])
  }
  zmat<-tmp
  tmp=NULL
  #Get the absolute max value to centre the scale around 0:
  lim<-max(abs(zmat[is.finite(zmat)]))
  print(lim)
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
  pheatmap::pheatmap(zmat, color=hmcol, show_colnames=T, cluster_cols=F, cluster_rows=F, main=title, labels_row=labgenes,breaks=seq(from=-lim, to=lim, length.out=100))
}
