#' Run ssGSEA and make a plot
#'
#' A function to run an ssGSEA analysis using normalized counts and create
#' a violin plot with pairwise comparisons
#'
#' @param counts Normalized count matrix - preferably z-score normalized by row. E.g. t(scale(t(counts(dds, normalized=T))))
#' @param geneset List with length > 1 of character vectors indicating gene sets of interest. Can be created using qusage::read.gmt("geneset.gmt")
#' @param method Method to employ in estimation of gene-set enrichment scores. One of "ssgsea","gsva","zscore","plage". Default is "ssgsea
#' @param condition Character vector of conditions in the same order they appear as columns in 'counts'.
#' @param title Character vector indicating the title of the plot
#' @param compare List of character vectors (each of length 2) indicating the pairwise comparisons to be made between conditions. Leave NULL for all comparisons to be made
#' @param col Character indicating the RColorBrewer palette name to be used. Default is "Dark2".
#' @return Violin plot of normalized enrichment scores for genesets between conditions
#' @export

gsva_plot<-function(counts, geneset, method="ssgsea", condition, title="ssGSEA", compare=NULL, col="Dark2"){
  progressBar<-txtProgressBar()
  #run the gsva
	res<-as.data.frame(GSVA::gsva(counts,geneset, method=method))
	#now we need to convert the results to a format acceptible for ggplot
	times<-dim(res)[1]
	conditions<-as.factor(c(rep(condition, each=times)))
	x<-res %>% tidyr::gather(key="Sample", value="NES") %>% dplyr::mutate(group=conditions) %>%
		dplyr::group_by(group)
	x<-as.data.frame(x)
	#Pairwise comparison
	pwc<- x%>% rstatix::pairwise_t_test(NES ~ group, comparisons=compare, p.adjust.method="BH")
	pwc <- pwc %>% rstatix::add_xy_position(x="group")
	#Now for the plot:
	ggpubr::ggviolin(x, x="group", y="NES", fill="group") +
	  ggplot2::geom_boxplot(width=0.1, fill="white") +
		ggpubr::stat_pvalue_manual(pwc, label="p.adj", tip.length=0, step.increase=0.1) +
	  ggplot2::labs(title=title, y="Normalized Enrichment Score", x="Condition") + ggplot2::theme_minimal() +
	  ggplot2::scale_fill_brewer(palette=col)
}
