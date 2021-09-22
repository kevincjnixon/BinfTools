rowGeoMean<-function(a){
  x<-NULL
  for(i in 1:nrow(a)){
    x<-c(x, prod(a[i,]^(1/length(a[i,]))))
  }
  return(x)
}
rowMedian<-function(a){
  x<-NULL
  for(i in 1:nrow(a)){
    x<-c(x, median(as.numeric(a[i,])))
  }
  return(x)
}

perMean<-function(counts, condition){
  counts<-counts/rowMeans(counts)
  head(counts)
  res<-NULL
  for(i in 1:length(levels(factor(condition)))){
    #print(which(condition %in% levels(factor(condition))[i]))
    res<-cbind(res, rowMeans(counts[,which(condition %in% levels(factor(condition))[i])]))
  }
  colnames(res)<-levels(factor(condition))
  res<-as.data.frame(res)
  return(res)
}
#' A function to make a plot of normalized counts
#'
#' This function takes normalized counts of specific genes from a DESeq2 counts
#' object, scales them, and creates a plot with pairwise t-tests by condition
#'
#' @param counts Normalized counts from a DESeq2 object - use 'counts(dds, normalized=T)'
#' @param scaling Method used to scale counts per gene across samples. 'zscore', 'log10', or 'none'. Default is 'zscore'
#' @param genes Character vector of genes to subset from counts. Must correspond with rownames(counts).
#' @param condition Character vector of conditions in DESeq2 object. Must be in order of columns (counts).
#' @param con Character indicating the control condition. (To be displayed first in plot). Default NULL.
#' @param title Character vector indicating title of plot. Defaults to "expression"
#' @param compare List of character vectors (each of length 2) indicating pairwise comparisons. If NULL, all possible comparisons will be made. Default is NULL
#' @param col Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used. Default is "Dark2"
#' @param method Character indicating what to plot. One of "ind", "mean", "geoMean", or "median", or "perMean". Defaults to "ind" for individual data points (one point per sample).
#' @param pair Boolean indicating if t-test should be independent (F; default) or paired (T).
#' @param pc Numeric indicating the pseudocount to be added when scaling="log10". Default=1.
#' @param yax Character indicating the y-axis label. Leave NULL if going with default axis label.
#' @param showStat Boolean indicating if statistics should be plotted.
#' @param style Character indicating the style of plot ("violin" or "box"). Defaults to "violin".
#' @param textsize Numeric indicating text size for the plot. Leave NULL for default.
#' @return Generates a violin or box plot
#' @export

count_plot<-function(counts, scaling="zscore", genes, condition, con=NULL, title="expression", compare=NULL, col="Dark2", method="ind", pair=F, pc=1, yax=NULL, showStat=T, style="violin", textsize=NULL){
  #Pull the normalized counts of genes
  res<-counts[which(rownames(counts) %in% genes),]
  ylab="z-score Normalized Expression"
  if(scaling=="log10"){
    #Log10 of counts
    res<-as.data.frame(log(pc+res, 10))
    ylab=expression(log[10](NormalizedExpression))
  }
  if(scaling=="zscore"){
    #zscore normalized counts
    res<-as.data.frame(t(scale(t(res))))
  }
  if(scaling=="none"){
    message("No scaling method selected...")
    ylab="Normalized Expression"
  }
  #Check the method
  if(method=="mean"){
    message("Calculaing mean across each condition...")
    tmp<-NULL
    for(i in 1:length(levels(factor(condition)))){
      tmp<-cbind(tmp,rowMeans(res[,which(condition %in% levels(factor(condition))[i])]))
    }
    colnames(tmp)<-as.character(levels(factor(condition)))
    res<-as.data.frame(tmp)
    condition<-as.character(levels(factor(condition)))
  }
  if(method=="geoMean"){
    message("Calculating geometric mean across each condition...")
    tmp<-NULL
    for(i in 1:length(levels(factor(condition)))){
      tmp<-cbind(tmp,rowGeoMean(res[,which(condition %in% levels(factor(condition))[i])]))
    }
    colnames(tmp)<-levels(factor(condition))
    res<-as.data.frame(tmp)
    condition<-as.character(levels(factor(condition)))
  }
  if(method=="perMean"){
    message("Calculating percent mean for each condition...")
    res<-perMean(res, condition)
    condition<-as.character(levels(factor(condition)))
  }
  if(method=="median"){
    message("Calculating median across each condition...")
    tmp<-NULL
    for(i in 1:length(levels(factor(condition)))){
      tmp<-cbind(tmp,rowMedian(res[,which(condition %in% levels(factor(condition))[i])]))
    }
    colnames(tmp)<-levels(factor(condition))
    res<-as.data.frame(tmp)
    condition<-as.character(levels(factor(condition)))
  } else {
    message("Using individual data points...")
  }
  #now we need to convert the results to a format acceptable for ggplot
  times<-dim(res)[1]
  conditions<-as.factor(c(rep(condition, each=times)))
  x<-res %>% tidyr::gather(key="Sample", value="Expression") %>% dplyr::mutate(group=conditions) %>%
    dplyr::group_by(group)
  if(!is.null(con)){
    if(length(con)==length(levels(factor(x$group)))){
      x<- x %>% dplyr::mutate(group=forcats::fct_relevel(group, con))
    } else {
      newlev<-c(con, levels(factor(x$group))[!which(levels(factor(x$group)) %in% con)])
      x<- x %>% dplyr::mutate(group=forcats::fct_relevel(group, newlev))
    }
    #x$group<-relevel(as.factor(x$group), con)
  }
  x<-as.data.frame(x)
  #print(head(x))
  #Remove any infinite values
  to_remove<-c()
  for(i in 1:nrow(x)){
    if(any(is.infinite(x[i,2]))){
      #print("infinite value found...")
      to_remove<-c(to_remove, i)
    }
  }
  if(length(to_remove > 1)){
    print(paste("Removing", length(to_remove), "inifite values..."))
    x<-x[-to_remove,]
  }
  x<-as.data.frame(x)
  #Pairwise comparison
  pwc<- x%>% rstatix::pairwise_t_test(Expression ~ group, comparisons=compare, p.adjust.method="BH", paired=pair)
  pwc <- pwc %>% rstatix::add_xy_position(x="group")
  #Now for the plot:
  if(!is.null(yax)){
    ylab<-yax
  }
  p<-NULL
  if(style=="violin"){
    p<- ggpubr::ggviolin(x, x="group", y="Expression", fill="group") +
      ggplot2::geom_boxplot(width=0.1, fill="white") +
      ggplot2::labs(title=title, y=ylab, x="Condition") + ggplot2::theme_minimal() + ggplot2::scale_fill_manual(values=colPal(col))+ ggplot2::theme(text=ggplot2::element_test(size=textsize))
  } else {
    p<- ggpubr::ggboxplot(x, x="group", y="Expression", fill="group") +
      ggplot2::labs(title=title, y=ylab, x="Condition") + ggplot2::theme_minimal() + ggplot2::scale_fill_manual(values=colPal(col))+ ggplot2::theme(text=ggplot2::element_test(size=textsize))
  }
  if(isTRUE(showStat)){
    p <- p + ggpubr::stat_pvalue_manual(pwc, label="p.adj", tip.length=0, step.increase=0.1)
  }
  print(p)
}
