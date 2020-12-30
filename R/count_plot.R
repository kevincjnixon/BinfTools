#' A function to make a violin plot of normalized counts
#'
#' This function takes normalized counts of specific genes from a DESeq2 counts
#' object, scales them, and creates a violin plot with pairwise t-tests by condition
#'
#' @param counts Normalized counts from a DESeq2 object - use 'counts(dds, normalized=T)'
#' @param scaling Method used to scale counts per gene across samples. Either 'zscore' or 'log10'. Default is 'zscore'
#' @param genes Character vector of genes to subset from counts. Must correspond with rownames(counts).
#' @param condition Factor vector of conditions in DESeq2 object. Must be in order of columns (counts). - use 'dds$Condition'
#' @param title Character vector indicating title of plot. Defaults to "expression"
#' @param compare List of character vectors (each of length 2) indicating pairwise comparisons. If NULL, all possible comparisons will be made. Default is NULL
#' @return Generates a violin plot
#' @export

count_plot<-function(counts, scaling="zscore", genes, condition, title="expression", compare=NULL){
	#Pull the normalized counts of genes
	res<-counts[which(rownames(counts) %in% genes),]
	ylab="z-score Normalized Expression"
	if(scaling=="log10"){
		#Log10 of counts
		res<-as.data.frame(log(res, 10))
		ylab=expression(log[10](NormalizedExpression))
	}
	if(scaling=="zscore"){
		#zscore normalized counts
		res<-as.data.frame(t(scale(t(res))))
	}
	#now we need to convert the results to a format acceptible for ggplot
	times<-dim(res)[1]
	conditions<-as.factor(c(rep(condition, each=times)))
	x<-res %>% tidyr::gather(key="Sample", value="Expression") %>% dplyr::mutate(group=conditions) %>%
		dplyr::group_by(group)
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
	pwc<- x%>% rstatix::pairwise_t_test(Expression ~ group, comparisons=compare, p.adjust.method="BH")
	pwc <- pwc %>% rstatix::add_xy_position(x="group")
	#Now for the plot:
	ggplot2::ggviolin(x, x="group", y="Expression", fill="group") +
	  ggplot2::geom_boxplot(width=0.1, fill="white") +
		ggpubr::stat_pvalue_manual(pwc, label="p.adj", tip.length=0, step.increase=0.1) +
	  ggplot2::labs(title=title, y=ylab, x="Condition") + ggplot2::theme_minimal()

}
