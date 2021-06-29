#' Run ssGSEA and make a plot
#'
#' A function to run an ssGSEA analysis using normalized counts and create
#' a plot with pairwise comparisons
#'
#' @param counts Normalized count matrix - preferably z-score normalized by row. E.g. t(scale(t(counts(dds, normalized=T))))
#' @param geneset List with length > 1 of character vectors indicating gene sets of interest. Can be created using qusage::read.gmt("geneset.gmt")
#' @param method Method to employ in estimation of gene-set enrichment scores. One of "ssgsea","gsva","zscore","plage". Default is "ssgsea
#' @param condition Character vector of conditions in the same order they appear as columns in 'counts'.
#' @param con Character indicating the control condition (Condition to be plotted first). Default is NULL.
#' @param title Character vector indicating the title of the plot
#' @param compare List of character vectors (each of length 2) indicating the pairwise comparisons to be made between conditions. Leave NULL for all comparisons to be made
#' @param col Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used. Default is "Dark2".
#' @param style Character indicating the style of plot ("violin" or "box"). Defaults to "violin".
#' @return Violin or box plot of normalized enrichment scores for genesets between conditions
#' @export

gsva_plot<-function(counts, geneset, method="ssgsea", condition, con=NULL, title="ssGSEA", compare=NULL, col="Dark2", style="violin"){
  progressBar<-txtProgressBar()
  #run the gsva
	res<-as.data.frame(GSVA::gsva(counts,geneset, method=method))
	#now we need to convert the results to a format acceptible for ggplot
	times<-dim(res)[1]
	conditions<-as.factor(c(rep(condition, each=times)))
	x<-res %>% tidyr::gather(key="Sample", value="NES") %>% dplyr::mutate(group=conditions) %>%
		dplyr::group_by(group)
	if(!is.null(con)){
	  x$group<-relevel(as.factor(x$group), con)
	}
	x<-as.data.frame(x)
	#Pairwise comparison
	pwc<- x%>% rstatix::pairwise_t_test(NES ~ group, comparisons=compare, p.adjust.method="BH")
	pwc <- pwc %>% rstatix::add_xy_position(x="group")
	#Now for the plot:
	p<-NULL
	if(style=="violin"){
	p<- ggpubr::ggviolin(x, x="group", y="NES", fill="group") +
	  ggplot2::geom_boxplot(width=0.1, fill="white") +
		ggpubr::stat_pvalue_manual(pwc, label="p.adj", tip.length=0, step.increase=0.1) +
	  ggplot2::labs(title=title, y="Normalized Enrichment Score", x="Condition") + ggplot2::theme_minimal() +
	  ggplot2::scale_fill_manual(values=colPal(col))
	} else {
	  p<- ggpubr::ggboxplot(x, x="group", y="NES", fill="group") +
	    ggpubr::stat_pvalue_manual(pwc, label="p.adj", tip.length=0, step.increase=0.1) +
	    ggplot2::labs(title=title, y="Normalized Enrichment Score", x="Condition") + ggplot2::theme_minimal() +
	    ggplot2::scale_fill_manual(values=colPal(col))
	}
	print(p)
}
