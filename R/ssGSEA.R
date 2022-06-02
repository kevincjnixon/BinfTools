#' Run ssGSEA and make a plot
#'
#' A function to run an ssGSEA analysis using normalized counts and create
#' a plot with pairwise comparisons
#'
#' @param counts Normalized count matrix - preferably z-score normalized by row. E.g. t(scale(t(counts(dds, normalized=T))))
#' @param geneset List with length > 1 of character vectors indicating gene sets of interest. Can be created using qusage::read.gmt("geneset.gmt")
#' @param method Method to employ in estimation of gene-set enrichment scores. One of "ssgsea","gsva","zscore","plage". Default is "ssgsea
#' @param stat.test Test to employ in comparison of enrichment scores across conditions. Can only be "t-test" or "wilcoxon". Default is "t-test".
#' @param condition Character vector of conditions in the same order they appear as columns in 'counts'.
#' @param con Character indicating the control condition (Condition to be plotted first). Default is NULL.
#' @param title Character vector indicating the title of the plot
#' @param compare List of character vectors (each of length 2) indicating the pairwise comparisons to be made between conditions. Leave NULL for all comparisons to be made
#' @param col Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used. Default is "Dark2".
#' @param style Character indicating the style of plot ("violin", "box", or "sina"). Defaults to "sina".
#' @param sinaPoint Character indicating colour of sina plot points. Defaults to "black".
#' @param showStat Boolean indicating if p-values should be plotted on figure
#' @param retRes Boolean indicating if GSVA results should be returned. Default=FALSE. If both retRes and retGP are TRUE, only GSVA results will be returned.
#' @param textsize Numeric indicating text size for plot. Leave NULL for default.
#' @param retGP Boolean indicating if ggplots2 object should be returned for further editing. Default=F.
#' @return Violin or box plot of normalized enrichment scores for genesets between conditions
#' @export

gsva_plot<-function(counts, geneset, method="ssgsea", stat.test = "t-test", condition, con=NULL, title="ssGSEA", compare=NULL, col="Dark2", style="sina", sinaPoint="black", showStat=T, retRes=F, textsize=NULL, retGP=F){
  if (!stat.test %in% c("t-test", "wilcoxon")){
    stop("Only t-test and wilcoxon tests are supported! Default is t-test.")
  }
  progressBar<-txtProgressBar()
  #run the gsva
	res<-as.data.frame(GSVA::gsva(as.matrix(counts),geneset, method=method))
	#now we need to convert the results to a format acceptible for ggplot
	times<-dim(res)[1]
	conditions<-as.factor(c(rep(condition, each=times)))
	x<-res %>% tidyr::gather(key="Sample", value="NES") %>% dplyr::mutate(group=conditions) %>%
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
	#Pairwise comparison
	if (stat.test == "t-test"){
	  pwc<- x%>% rstatix::pairwise_t_test(NES ~ group, comparisons=compare, p.adjust.method="BH")
	}
	if (stat.test == "wilcoxon"){
	  pwc<- x%>% rstatix::pairwise_wilcox_test(NES ~ group, comparisons=compare, p.adjust.method="BH")
	}
	pwc <- pwc %>% rstatix::add_xy_position(x="group")
	#Now for the plot:
	p<-NULL
	if(style=="violin"){
	p<- ggpubr::ggviolin(x, x="group", y="NES", fill="group") +
	  ggplot2::geom_boxplot(width=0.1, fill="white") +
	  ggplot2::labs(title=title, y="Normalized Enrichment Score", x="Condition") + ggplot2::theme_minimal() +
	  ggplot2::scale_fill_manual(values=colPal(col)) + ggplot2::theme(text=ggplot2::element_text(size=textsize))
	} 
	if(style=="sina"){
	  p<- ggpubr::ggviolin(x, x="group", y="NES", fill="group") +
	    ggforce::geom_sina(ggplot2::aes(colour=as.factor(x$group)), alpha=0.5) + ggplot2::scale_colour_manual(values=colPal(sinaPoint)) +
	    ggplot2::labs(title=title, y="Normalized Enrichment Score", x="Condition") + ggplot2::theme_minimal() +
	    ggplot2::scale_fill_manual(values=colPal(col)) + ggplot2::theme(text=ggplot2::element_text(size=textsize)) + ggplot2::guides(colour="none")
	}
	if(style=="box") {
	  p<- ggpubr::ggboxplot(x, x="group", y="NES", fill="group") +
	    ggplot2::labs(title=title, y="Normalized Enrichment Score", x="Condition") + ggplot2::theme_minimal() +
	    ggplot2::scale_fill_manual(values=colPal(col))+ ggplot2::theme(text=ggplot2::element_text(size=textsize))
	}
	if(isTRUE(showStat)){
	  p<-p+ggpubr::stat_pvalue_manual(pwc, label="p.adj", tip.length=0, step.increase=0.1)
	}
	if(isTRUE(retRes)){
	  print(p)
	  return(res)
	}
	if(isTRUE(retGP)){
	  return(p)
	} else {
	  print(p)
	}
}
