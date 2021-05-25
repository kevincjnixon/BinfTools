#'Interactive MA Plot
#'
#'This function creates an interactive MA plot from a DESeq2 results object. It is like
#''plotMA()' from DESeq2, but more customizable and interactive.
#'
#'@param res A DESeq2 results object obtained from 'results(dds)' or a data.frame with the same column name values as a DESeq2 results object and rownames as genes
#'@param title A character vector indicating the title of the plot
#'@param p A number indicating the threshold for 'padj' where padj<p are significant genes. Should not be used if using 'pval'
#'@param pval A number indicating the threshold for 'pvalue' where pvalue<pval are significant genes. Should not be used if using 'p'
#'@param FC A number indicating the log2FoldChange threshold where abs(log2FC)>FC are significant genes. Default is 1 - can be 0 if not using fold-change threshold.
#'@return An interactive MA plot with x-axis indicating gene expression 'baseMean' and y-axis indicating log2FoldChange. Blue dots are downregulated genes, red dots are upregulated genes.
#' @export

interMA<-function(res, title="MA Plot", p=NULL, pval=NULL, FC=1){
  res$colour<-rep("No Change", nrow(res))
  to_remove<-c()
  if(!is.null(p)){
    for(i in 1:nrow(res)){
      if(!is.na(res$padj[i])){
        if(res$log2FoldChange[i] >= FC && res$padj[i]<p){
          res$colour[i]<-"Up"
        }
        if(res$log2FoldChange[i] <= -FC && res$padj[i]<p){
          res$colour[i]<-"Down"
        }
      } else {
        to_remove<-c(to_remove, i)
      }
    }
  } else {
    for(i in 1:nrow(res)){
      if(!is.na(res$pvalue[i])){
        if(res$log2FoldChange[i] >= FC && res$pvalue[i]<pval){
          res$colour[i]<-"Up"
        }
        if(res$log2FoldChange[i] <= -FC && res$pvalue[i]<pval){
          res$colour[i]<-"Down"
        }
      } else{
        to_remove<-c(to_remove,i)
      }
    }
  }
  if(length(to_remove)>0){
    res<-res[-to_remove,]
  }
  vline<-function(val=0, color="green", col=column, mat=res){
    addon<-0.1*diff(range(mat[,`col`]))
    list(
      type="line",
      y0=min(mat[,`col`])-addon,
      y1=max(mat[,`col`])+addon,
      x0=val,
      x1=val,
      line=list(dash="dash",color=color)
    )
  }
  hline<-function(y=0, color="black"){
    list(
      type="line",
      x0=0,
      x1=1,
      xref="paper",
      y0=y,
      y1=y,
      line=list(color=color)
    )
  }
  q<-quantile(log10(res$baseMean))
  names(q)<-c()
  pal<-c("blue","black","red")
  fig<-plotly::plot_ly(res, type="scatter", x=~log10(baseMean), y=~log2FoldChange,
                       text=~rownames(res), color=~colour, colors=pal) %>%
    plotly::layout(shapes=list(hline(-FC), hline(FC))) %>% #, vline(q[1]),
                               #vline(q[2]), vline(q[3]), vline(q[4]), vline(q[5]))) %>%
    plotly::layout(title=title, xaxis=list(title="log10 normalized expression", showgrid=F,
                                           showline=T, zeroline=F, ticks="outside"),
                   yaxis=list(title="log2 Fold-Change Expression", showgrid=F,
                              showline=T, zeroline=F, ticks="outside"))

  print(fig)
  #return(res)
}


#'Interactive Volcano Plot
#'
#'This function creates an interactive volcano plot from a DESeq2 results object.
#'
#'@param res A DESeq2 results object obtained from 'results(dds)' or a data.frame with the same column name values as a DESeq2 results object and rownames as genes
#'@param title A character vector indicating the title of the plot
#'@param p A number indicating the threshold for 'padj' where padj<p are significant genes. Should not be used if using 'pval'
#'@param pval A number indicating the threshold for 'pvalue' where pvalue<pval are significant genes. Should not be used if using 'p'
#'@param FC A number indicating the log2FoldChange threshold where abs(log2FC)>FC are significant genes. Default is 1 - can be 0 if not using fold-change threshold.
#'@return An interactive volcano plot with x-axis indicating log2FoldChange and y-axis indicating significance Blue dots are downregulated genes, red dots are upregulated genes.
#' @export

interVP<-function(res, title="Volcano Plot", p=NULL, pval=NULL, FC=1){
  res$colour<-rep("No Change", nrow(res))
  to_remove<-c()
  sig<-NULL
  if(!is.null(p)){
    #Set sig to the p-value threshold
    sig<--log(res[which(res$padj == max(subset(res, padj<0.05)$padj)),]$pvalue[1],10)
    for(i in 1:nrow(res)){
      if(!is.na(res$padj[i])){
        if(res$log2FoldChange[i] >= FC && res$padj[i]<p){
          res$colour[i]<-"Up"
        }
        if(res$log2FoldChange[i] <= -FC && res$padj[i]<p){
          res$colour[i]<-"Down"
        }
      } else {
        to_remove<-c(to_remove, i)
      }
    }
  } else {
    sig<- -log10(pval)
    for(i in 1:nrow(res)){
      if(!is.na(res$pvalue[i])){
        if(res$log2FoldChange[i] >= FC && res$pvalue[i]<pval){
          res$colour[i]<-"Up"
        }
        if(res$log2FoldChange[i] <= -FC && res$pvalue[i]<pval){
          res$colour[i]<-"Down"
        }
      } else{
        to_remove<-c(to_remove,i)
      }
    }
  }
  if(length(to_remove)>0){
    res<-res[-to_remove,]
  }
  vline<-function(val=0, color="black", mat=res){
    addon<-0.1*diff(range(-log10(mat$pvalue)))
    list(
      type="line",
      y0=min(-log10(mat$pvalue))-addon,
      y1=max(-log10(mat$pvalue))+addon,
      x0=val,
      x1=val,
      line=list(dash="dash",color=color)
    )
  }
  hline<-function(y=0, color="black"){
    list(
      type="line",
      x0=0,
      x1=1,
      xref="paper",
      y0=y,
      y1=y,
      line=list(color=color)
    )
  }
  q<-quantile(log10(res$baseMean))
  names(q)<-c()
  pal<-c("blue","black","red")
  fig<-plotly::plot_ly(res, type="scatter", x=~log2FoldChange, y=~-log10(pvalue),
                       text=~rownames(res), color=~colour, colors=pal) %>%
    plotly::layout(shapes=list(vline(-FC), vline(FC), hline(sig))) %>% #, vline(q[1]),
    #vline(q[2]), vline(q[3]), vline(q[4]), vline(q[5]))) %>%
    plotly::layout(title=title, xaxis=list(title="log2 Fold-Change Expression", showgrid=F,
                                           showline=T, zeroline=F, ticks="outside"),
                   yaxis=list(title="-log 10 p-value", showgrid=F,
                              showline=T, zeroline=F, ticks="outside"))

  print(fig)
  #return(res)
}
