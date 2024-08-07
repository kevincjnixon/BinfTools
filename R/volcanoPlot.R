#'Custom volcano Plot
#'
#'This function creates a volcano plot from a DESeq2 results object.
#'
#'@param res A DESeq2 results object obtained from 'results(dds)' or a data.frame with the same column name values as a DESeq2 results object and rownames as genes
#'@param title A character vector indicating the title of the plot
#'@param p A number indicating the threshold for 'padj' where padj<p are significant genes. Should not be used if using 'pval'
#'@param pval A number indicating the threshold for 'pvalue' where pvalue<pval are significant genes. Should not be used if using 'p'
#'@param FC A number indicating the log2FoldChange threshold where abs(log2FC)>FC are significant genes. Default is 1 - can be 0 if not using fold-change threshold.
#'@param lab A character vector of genes to be labeled on the plot. Should correspond with rownames(res). Default is NULL. Enter 'labDEG' to label all DEGs given thresholds.
#'@param col A character vector of genes to be coloured orange on the plot. Should correspond with rownames(res). Default is NULL.
#'@param fclim A number indicating the maximum log2FoldChange value desired on the volcano plot (x-axis limits). Points outside this limit will appear as diamonds on the limits.
#'@param showNum A boolean indicating whether gene numbers should be displayed on the plot. Default is TRUE.
#'@param returnDEG A boolean indicating whether DEGs (using given thresholds) should be returned as a list (down, up)
#'@param ordBy Character indicating how DEGs should be ordered: "l2FC" order by descending absolute log2FoldChange, "sig" order by descending significance. Default is "sig".
#'@param expScale A boolean indicating if points should be scaled using average expression (baseMean) - higher expressed = larger point. Default=FALSE.
#'@param upcol A character vector indicating the colour (colour name or hexadecimal) of 'upregulated' genes. leave NULL for default red.
#'@param dncol A character vector indicating the colour (colour name or hexadecimal) of 'downregulated' genes. leave NULL for default blue.
#'@param labcol A character vector indicating the colour (colour name or hexadecimal) of genes to be specifically highlighted by the 'col' and 'lab' arguments. default is "orange".
#'@param labColumn A character vector indicating the column name (if left NULL, 'rownames(res)' will be used) to be used as labels. The gene names provided in the 'lab' argument must still correspond to rownames(res).
#'@return A volcano plot with x-axis indicating log2FoldChange expression and y-axis indicating -log10 pvalue. Blue dots are downregulated genes, red dots are upregulated genes.
#' @export
#'

volcanoPlot<-function(res, title, p=NULL, pval=NULL, FC=1, lab=NULL, col=NULL, fclim=NULL, showNum=TRUE, returnDEG=FALSE, ordBy="sig", expScale=FALSE, upcol=NULL, dncol=NULL, labcol="orange", labColumn=NULL){
  #Function to create a volcano plot from a results table (res)
  #p and pval are mutually exclusive and describe the adjusted pvalue or unadjusted pvalue thresholds of DEGs, respectively
  #FC is the absolute log2 fold-change threshold for a DEG
  #lab is a character vector of gene names (matching the rownames of res) that are to be labelled
  #col is a character vector of gene names (matching the rownames of res) that are to be coloured orange
  #If there are p-values = 0, remove these because they will cause an error with the -log p-value calculation, but save them in zeroes to include in the numbers of genes
  zeroes<-NULL
  #Set up variables to contain down/up regulated genes and genes with no change (nc)
  down<-NULL
  up<-NULL
  nc<-NULL
  pNA<-NULL
  #A variable to contain the coordinates of plots for labelling (if specified)
  labpoints<-NULL
  #A variable to store the significance threshold for later plotting:;
  sig<-NULL
  #print(head(res))
  extremes<-NULL
  #Figure out the largest abolute log fold-change so we can set the x-axis accordingly
  maxvalx<-fclim
  if(is.null(fclim)){
    if(!is.null(p)){
      #Ensure if using adjusted p-value, there are no fold-changes from padj=NA included
      maxvalx<-max(abs(na.omit(res$log2FoldChange[complete.cases(res$padj)])))
    }
    if(!is.null(pval)){
      #Ensure if using p-value, there are no fold-changes from pvalue=NA included
      maxvalx<-max(abs(na.omit(res$log2FoldChange[complete.cases(res$pvalue)])))
    }
    extremes<-0
  }
  if(!is.null(p)){
    #Set sig to the p-value threshold
    sig<--log(res[which(res$padj == max(subset(res, padj<p)$padj)),]$pvalue[1],10)
    #Separate out padj = 0
    zeroes<-subset(res, padj == 0)
    #Get the number of NA values for P values
    pNA<-length(which(is.na(res$padj) == TRUE))
    #Remove padj=0 from res
    res<-subset(res, padj > 0)
    #Subset out the down/up regulated genes and no change genes
    down<-subset(res,padj<p & log2FoldChange < -FC)
    up<-subset(res,padj<p & log2FoldChange > FC)
    nc<-subset(res,padj>p | abs(log2FoldChange) < FC)
    #Subset the extreme values (if necessary)
    if(is.null(extremes)){
      extremes<-subset(res, abs(log2FoldChange)>maxvalx)
      for(i in 1:nrow(extremes)){
        if(extremes$log2FoldChange[i] < 0){
          extremes$log2FoldChange[i]<-maxvalx*-1
        } else {
          extremes$log2FoldChange[i]<-maxvalx
        }
      }
    }
  }
  if(!is.null(pval)){
    #Set sig to the p-value threshold
    sig<--log(pval,10)
    #Separate out padj = 0
    zeroes<-subset(res, pvalue == 0)
    #Get the number of NA values for P values
    pNA<-length(which(is.na(res$pvalue) == TRUE))
    #Remove padj=0 from res
    res<-subset(res, pvalue > 0)
    #Subset out the down/up regulated genes and no change genes
    down<-subset(res,pvalue<pval & log2FoldChange < -FC)
    up<-subset(res,pvalue<pval & log2FoldChange > FC)
    nc<-subset(res,pvalue>pval | abs(log2FoldChange) < FC)
    #Subset the extreme values (if necessary)
    if(is.null(extremes)){
      extremes<-subset(res, abs(log2FoldChange)>maxvalx)
      for(i in 1:nrow(extremes)){
        if(extremes$log2FoldChange[i] < 0){
          extremes$log2FoldChange[i]<-maxvalx*-1
        } else {
          extremes$log2FoldChange[i]<-maxvalx
        }
      }
    }
  }
  #Set a column for point sizes:
  nc$cex<-rep(0.6, nrow(nc))
  up$cex<-rep(0.6, nrow(up))
  down$cex<-rep(0.6, nrow(down))
  cex_scale<-(diff(range(log10(res$baseMean[!is.na(res$baseMean)]))))
  if(isTRUE(expScale)){
    nc$cex<- log10(nc$baseMean)/cex_scale
    up$cex<- log10(up$baseMean)/cex_scale
    down$cex<- log10(down$baseMean)/cex_scale
  }

  if(is.null(upcol)){
    upcol <- rgb(1,0,0,0.75)
  }
  if(is.null(dncol)){
    dncol <- rgb(0,0,1,0.75)
  }

  #Plot the points, start with nc (black), then down (blue), and up (red)
  plot(nc$log2FoldChange,-log(nc$pvalue,10), pch=16, cex=nc$cex, main=title, xlab=expression(log[2](FoldChange)), ylab=expression(-log[10](pvalue)),
       ylim=c(0, max(c(10,max(na.omit(-log(res$pvalue,10)))))), xlim=c(-maxvalx, maxvalx), col=rgb(0,0,0,0.5))
  points(down$log2FoldChange,-log(down$pvalue,10), pch=16, cex=down$cex, col=dncol)
  points(up$log2FoldChange, -log(up$pvalue,10),pch=16, cex=up$cex, col=upcol)
  #Add in the extreme points (if necessary)
  if(length(extremes) >1){
    #Make sure if they are DEGs, they're coloured properly:
    tmp<-extremes[which(rownames(extremes)%in% rownames(nc)),]
    points(tmp$log2FoldChange, -log(tmp$pvalue,10), pch=18, cex=1, col=rgb(0,0,0,0.5))
    tmp<-extremes[which(rownames(extremes) %in% rownames(down)),]
    points(tmp$log2FoldChange, -log(tmp$pvalue,10), pch=18, cex=1, col=dncol)
    tmp<-extremes[which(rownames(extremes) %in% rownames(up)),]
    points(tmp$log2FoldChange, -log(tmp$pvalue,10), pch=18, cex=1, col=upcol)
    #Check to see if any extremems need to be labelled or coloured:
    if(!is.null(lab)){
      labpoints<-extremes[lab,]
      labs<-rownames(labpoints)
      if(!is.null(labColumn)){
        labs<-labpoints[,which(colnames(labpoints) %in% labColumn)]
      }
      text(labpoints$log2FoldChange, -log(labpoints$pvalue,10), labels=labs, cex=0.75, font=2,
           pos=4, col=labcol)
    }
    if(!is.null(col)){
      colpoints<-extremes[col,]
      points(colpoints$log2FoldChange, -log(colpoints$pvalue,10), pch=18, cex=1, col=labcol)
    }
  }
  #Plot dotted lines for the thresholds
  if(FC>0){
	abline(v=c(-(FC),FC),lty=c(2,2))
  }
  abline(h=c(sig),lty=c(2))
  if(!is.null(col)){
    colpoints<-res[col,]
    colpoints$cex<-rep(0.6, nrow(colpoints))
    if(isTRUE(expScale)){
      colpoints$cex<- log10(colpoints$baseMean)/cex_scale
    }
    points(colpoints$log2FoldChange, -log(colpoints$pvalue,10), pch=16, cex=colpoints$cex, col=labcol)
  }
  #Now calculate number of DEGs (including padj=0)
  numdown<-dim(down)[1]+dim(subset(zeroes, log2FoldChange < -FC))[1]
  numup<-dim(up)[1]+dim(subset(zeroes, log2FoldChange > FC))[1]
  numnc<-dim(nc)[1]+dim(subset(zeroes, abs(log2FoldChange)<FC))[1]+pNA
  #Get the DEGs
  DEdown <- rownames(down)
  DEup <- rownames(up)
  #Print the numbers of genes on the plot
  if(isTRUE(showNum)){
    text(-(maxvalx/2), max(c(10,max(na.omit(-log(res$pvalue,10))))), numdown, col = dncol)
    text((maxvalx/2), max(c(10,max(na.omit(-log(res$pvalue,10))))), numup, col= upcol)
    text(-maxvalx, 0, numnc, col="black", adj=c(0,0))
    text(maxvalx, 0, paste("Total:",(dim(res)[1])+dim(zeroes)[1]+pNA), col="purple", pos=2)
  }
  #Now add in the gene labels if they were specified
  if(!is.null(lab)){
    if(lab[1] == "labDEG"){
      lab=c(DEdown,DEup)
    }
    labpoints<-res[lab,]
    labs<-rownames(labpoints)
    if(!is.null(labColumn)){
      labs<-labpoints[,which(colnames(labpoints) %in% labColumn)]
    }
    text(labpoints$log2FoldChange, -log(labpoints$pvalue,10), labels=labs, cex=0.75, font=2,
         pos=4, col=labcol)
  }
  #Return a list of the numbers of genes (down, then up, then no change)
  if(isFALSE(returnDEG)){
    return(list(Down=numdown,Up=numup,No_Change=numnc))
  }
  if(isTRUE(returnDEG)){
    ord<-rownames(res)[order(abs(res$stat), decreasing=TRUE)]
    if(ordBy == "l2FC"){
      ord<-rownames(res)[order(abs(res$log2FoldChange), decreasing=TRUE)]
    }
    toRet<-list(Down=c(DEdown,rownames(subset(zeroes, log2FoldChange < -FC))), Up=c(DEup, rownames(subset(zeroes, log2FoldChange > FC))))
    toRet$Down<-BinfTools:::rmNA(toRet$Down[match(ord, toRet$Down)])
    toRet$Up<-BinfTools:::rmNA(toRet$Up[match(ord, toRet$Up)])
    return(toRet)
  }
}

#3D volcano plot - for combining two different contrasts - pvalues are combined using Fisher's combined probability test
#' Make a 3D volcano plot
#' combine 2 results objects into one volcano plot - note this is experimental and not usually done
#' @param x DESeq2 results object to be plotted for the x-axis
#' @param y DESeq2 results object to be plotted for the y-axis
#' @param xlab Contrast name for x. Default='res1'.
#' @param ylab Contrast name for y. Default='res2'.
#' @param pval significance threshold for differential genes. Default is 0.05
#' @param pcol Character of either 'pvalue' or 'padj' indicating which column to use to plot significance. Default is "pvalue".
#' @param usep numeric of either 1 or 2 indicating which p-values to plot (x or y, respectively). If left NULL (default). P-values from both x and y are combined using Fisher's combined probability test, and the resulting p-values are plotted.
#' @return 3D volcano plot with every combination of up/down-regulated genes labeled and a data frame of the combined data used to generate the plot.
#' @export

vol3d<-function(x, y, xlab="res1", ylab="res2", pval=0.05, pcol="pvalue", usep=NULL){
  x<-x[order(rownames(x)),]
  y<-y[order(rownames(y)),]
  x<-x[which(rownames(x) %in% rownames(y)),]
  y<-y[which(rownames(y) %in% rownames(x)),]
  zlab<-pcol
  if(!all.equal(rownames(x), rownames(y))){
    stop("Check again!")
  }
  pcol<-which(colnames(x) %in% pcol)
  pvals<- -log10(x[,pcol])
  if(is.null(usep)){
    pvals<-c()
    for(i in 1:nrow(x)){
      pvals<-c(pvals, survcomp::combine.test(c(x[,pcol][i], y[,pcol][i])))
    }
    pvals<- -log10(pvals)
  } else {
    if(usep[1] == 2){
      pvals<- -log10(y$pvalue)
    }
  }
  for3d<-data.frame(x=x$log2FoldChange, y=y$log2FoldChange, z=pvals, row.names=rownames(x))
  for3d$col<-ifelse(for3d$x > 0 & for3d$y >0 & for3d$z > -log10(pval), rgb(1,0,0,0.5), rgb(0,0,0,0.5))
  for3d$DE<-ifelse(for3d$x > 0 & for3d$y >0 & for3d$z > -log10(pval), "Up Both", "No Change")
  for3d$col[which(for3d$x <0 & for3d$y <0 & for3d$z > -log10(pval))] <- rgb(0,0,1,0.5)
  for3d$DE[which(for3d$x <0 & for3d$y <0 & for3d$z > -log10(pval))] <- "Down Both"
  for3d$col[which(for3d$x <0 & for3d$y >0 & for3d$z > -log10(pval))] <- rgb(0,1,0,0.5)
  for3d$DE[which(for3d$x <0 & for3d$y >0 & for3d$z > -log10(pval))] <- paste0("Down ", xlab,", Up ", ylab)
  for3d$col[which(for3d$x >0 & for3d$y <0 & for3d$z > -log10(pval))] <- rgb(0.5,0,1,0.5)
  for3d$DE[which(for3d$x >0 & for3d$y <0 & for3d$z > -log10(pval))] <- paste0("Up ",xlab,", Down ",ylab)
  scatterplot3d::scatterplot3d(for3d[,c(1:3)], pch=16, xlab=paste("log2FoldChange -", xlab), ylab=paste("log2FoldChange -", ylab), zlab=paste0("-log10(",zlab,")"), color=for3d$col,
                               grid=TRUE, box=FALSE)
  legend("topright", legend=c("Up Both", "Down Both", paste0("Up ",xlab,", Down ", ylab), paste0("Down ",xlab,", Up ",ylab), "No Change"),
         col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5), rgb(0.5,0,1,0.5), rgb(0,1,0,0.5), rgb(0,0,0,0.5)), pch=16)
  return(for3d)
}
