#' Calculate average gene expression between replicates
#'
#' @param counts Data.frame of counts (rows=genes, columns=samples)
#' @param cond character vector of length ncol(counts) assigning each column to a condition
#' @param method character of either "mean", "median", or "geoMean" of how to calculate average expression. Default is "mean"
#' @return Data frame with rows=genes and columns=condition with indicated average of gene expression beween replicates
#' @export

avgExp<-function(counts, cond, method=c("mean","median","geoMean")){
  levs<-unique(cond)
  x<-c()
  if(method[1]=="mean"){
    for(i in 1:length(levs)){
      x<-cbind(x, rowMeans(counts[,which(cond %in% levs[i])]))
    }
  }
  if(method[1]=="median"){
    for(i in 1:length(levs)){
      x<-cbind(x, Biobase::rowMedians(counts[,which(cond %in% levs[i])]))
    }
  }
  if(method[1]=="geoMean"){
    for(i in 1:length(levs)){
      x<-cbind(x, BinfTools:::rowGeoMean(counts[,which(cond %in% levs[i])]))
    }
  }
  rownames(x)<-rownames(counts)
  colnames(x)<-levs
  return(as.data.frame(x))
}

#' Sample correlation heatmap
#'
#' @param counts data.frame object with rows=genes, and columns=samples
#' @param method character indicating correlation method, default is "spearman"
#' @param title character indicating plot title
#' @param hmcol Color ramp palette to use for heatmap colours. Defaults to 'RdYlBu' colours
#' @param showNum Boolean indicating if correlation coefficients should be displayed on heatmap. Default=TRUE.
#' @param showTree Boolean indicating if hierarchical clustering trees should be shown. Default=F
#' @param annodf Data frame with rownames=sample names and columns of annotation data to be plotted. Leave NULL if none.
#' @param annoCols Colour palette for annotations
#' @return Heatmap with clustered rows/columns of correlation between samples
#' @export

corHeat<-function(counts, method="spearman", title="Spearman Correlation", hmcol=NULL, showNum=T, showTree=F, annodf=NULL, annoCols="Dark2"){
  M<-cor(counts, method=method)
  if (is.null(hmcol)) {
    hmcol <- colorRampPalette(rev(colPal("RdYlBu")))(100)
  }
  treeheight=50
  if(isFALSE(showTree)){
    treeheight=0
  }
  if(is.null(annodf)){
    annoCols=NULL
  } else {
    cols<-colPal(annoCols)
    annoCols<-list()
    for(i in 1:ncol(annodf)){
      annodf[,i]<-factor(annodf[,i])
      annoCols[[i]]<-cols[1:length(levels(annodf[,i]))]
      names(annoCols[[i]])<-levels(annodf[,i])
      names(annoCols)[i]<-colnames(annodf)[i]
    }
  }
  
  lim<-NULL
  if(min(M)<0){
    lim<-seq(from=-1, to=1, length.out=100)
  }
  pheatmap::pheatmap(M, color=hmcol, main=title, display_numbers=showNum, treeheight_col = treeheight, treeheight_row = treeheight,
                     annotation_row=annodf, annotation_col=annodf, annotation_colors = annoCols, breaks=lim)
}

## colInds is indices of columns in data.frame
custPairs <- function(colInds, data=iris, title="test", m="spearman") {
  n <- length(colInds)
  cols <- expand.grid(names(data)[colInds], names(data)[colInds])
  cInds <- unlist(mapply(function(a, b, c) a*n+b:c, 0:max(0,n-2), 2:n, rep(n, n-1)))
  cols <- cols[cInds,]  # indices will be in column major order
  #print(cols)
  ## These parameters are applied to each plot we create
  #pars <- list(geom_point(alpha=0.8, color="black"),
  #             geom_smooth(method=m, color="red", lwd=1))

  ## Create the plots (dont need the lower plots in the ggpairs call)
  #plots <- apply(cols, 1, function(cols)
  #  ggplot(data[,cols], aes_string(x=cols[2], y=cols[1])) + pars)
  plots<-apply(cols, 1, function(cols)
    ggpubr::ggscatter(data, x=cols[1], y=cols[2],
                      add="reg.line", conf.int=T,
                      cor.coef=F, cor.method=m,
                      add.params=list(color="red", fill="lightpink")))
  gg <- GGally::ggpairs(data[, colInds],
                diag=list(continuous="blankDiag"),
                upper=list(continuous=GGally::wrap("cor", method=m)),
                title=title)

  rowFromTop <- unlist(mapply(`:`, 2:n, rep(n, n-1)))
  colFromLeft <- rep(1:(n-1), times=(n-1):1)
  for (i in seq_along(plots))
    gg <- GGally::putPlot(gg, plots[[i]], rowFromTop[i], colFromLeft[i])
  return( gg )
}

corPlots<-function(avgmat, title="Correlation", method="spearman", transform="none", printInd=T){
  spPlot<-function(dat, x, y, xlab, ylab, title, method){
    g<-ggpubr::ggscatter(dat, x=x, y=y,
                         add="reg.line", conf.int=T,
                         cor.coef=T, cor.method=method,
                         add.params=list(color="red", fill="red"),
                         xlab=xlab, ylab=ylab, title=title)
    print(g)
  }
  pre<-"Mean Normalized Expression -"
  if(transform=="log10"){
    avgmat<-log10(1+avgmat)
    pre<-paste("Log10",pre)
  }
  if(transform=="log2"){
    avgmat<-log2(1+avgmat)
    pre<-paste("Log2",pre)
  }


  print(custPairs(colInds=c(1:ncol(avgmat)), data=avgmat, title=title, m=method))
  if(isTRUE(printInd)){
    for(i in 1:ncol(avgmat)){
      if(i<ncol(avgmat)){
        for(j in (i+1):ncol(avgmat)){
          x<-colnames(avgmat)[i]
          y<-colnames(avgmat)[j]
          message("Plotting",x,"vs",y,"...")
          spPlot(avgmat, x, y, xlab=paste(pre,x), ylab=paste(pre,y), title, method)
        }
      }
    }
  }
}

#'Expresion correlation plots of all gene expressions
#'
#'@param counts data.frame object with rows=genes and columns=samples of gene expression
#'@param cond character vector of length ncol(counts) assigning each sample to a condition
#'@param title character indicating the title of the plot(s)
#'@param method character indicating correlation method. Default is 'spearman'
#'@param transform character indicating if gene expression should be transformed. Default is "none", other options are "log2" and "log10"
#'@param printInd Boolean indicating if individual correlation plots should be printed in addition to paired plots. Default is FALSE.
#'@return correlation plots
#'@export

expCor<-function(counts, cond, title, method="spearman", transform="none", printInd=F){
  avgmat<-avgExp(counts, cond)
  corPlots(avgmat, title, method, transform, printInd)
}
