
geoMean<-function(a){prod(a)^(1/length(a))}

#'Normalize gene expression to specific genes
#'
#'This function will normalize gene expression in a count matrix relative to a given gene
#'or vector of genes.
#'
#' @param counts count matrix (normalized or raw) with genes as rows and samples as columns
#' @param norm character vector of length >=1 of gene names (matching rownames(counts)) to use as reference genes for normalization. If length > 1, geometric mean of genes' expressions will be used.
#' @return data frame of dim(counts) of gene expression normalized relative to 'norm'
#' @export

normToGene<-function(counts, norm){
  #factor to normalize to
  fac<-c()
  if(length(norm)>1){
    #Calculate the geometric of the genes in each column (sample)
    tmp<-counts[which(rownames(counts) %in% norm),]
    for(i in 1:ncol(tmp)){
      fac<-c(fac, geoMean(tmp[,i]))
    }
  } else {
    fac<-counts[which(rownames(counts) %in% norm),]
  }
  #Normalize each column using each factor
  res<-counts[,1]/as.numeric(fac[1])
  for(i in 2:length(fac)){
    res<-cbind(res, counts[,i]/as.numeric(fac[i]))
  }
  rownames(res)<-rownames(counts)
  colnames(res)<-colnames(counts)
  return(as.data.frame(res))
}

#This function is for summarizing data and will be called in the next function
#This will calculate the mean and sd/se allelic reads for each type of cancer
#In some cases, there will be no sd/se (only one value/sample)
data_sum<-function(data, eb){
  #data is the data table and eb is what we want the error bars to be (sd or se)
  summary_func<-function(x, col, eb){
    sum<-NULL
    if(eb=="sd"){ #If we want sd, calcaulate sd
      #message("Calculating mean and standard deviation...")
      sum<-c(mean=mean(as.numeric(x[[col]]), na.rm=TRUE),
             eb=sd(as.numeric(x[[col]]), na.rm=T))
    }
    if(eb=="se"){ #If we want se, calculate se
      #message("Calculating mean and standard error...")
      sum<-c(mean=mean(as.numeric(x[[col]]), na.rm=T),
             eb=(sd(as.numeric(x[[col]]), na.rm=T)/sqrt(length(x[[col]]))))
    }
	if(eb==0){ #We don't want any error bars
		#message("eb=0. Error bars will not be showns...")
		sum<-c(mean=mean(as.numeric(x[[col]]), na.rm=T),
			  eb=0)
	}
    return(sum)
  }
  if(eb==0){ #We don't want any error bars
    message("eb=0. Error bars will not be shown...")
  }
  data_sum<-plyr::ddply(data, c("group","gene"), .fun=summary_func, "expression", eb)
  #data_sum<-rename(data_sum, c("mean"="reads"))
  #Set NA values for error bars (eb) to 0
  #data_sum$eb[which(is.na(data_sum$eb))]<-0
  return(data_sum)
}

#'Create a bar plot of gene expression
#'
#'This function will accept a list of genes and counts to produce a bar plot of gene expression in different conditions
#'
#'@param genes character vector of genes to plot. Must match rownames(counts)
#'@param counts count matrix of gene expression, where rows are genes and columns are samples
#'@param conditions character vector of length ncol(counts) describing the condition of each sampmle in counts
#'@param title character describing the title of the plot
#'@param norm character describing the condition to use as a reference for relative gene expression. Each gene's expression will be set relative to this condition. Leave 'NULL' if plotting raw values.
#'@param eb character describing the style of error bar. "sd" = standard deviation, "se" = standard error of the mean. Use '0' if no error bars to be plotted. Default is 'sd'.
#'@param returnDat Boolean indicating if list of raw and summarized data should be returned for further analysis. Default is FALSE.
#'@param col character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used. Default is "Dark2"
#'@param ord character indicating the order in which the samples should appear (overrides any ordering from using 'norm' argument). Default is NULL.
#'@return Bar plot of gene expression and list of length 2 containing 'rawData' and 'Summary' of gene expression data if 'returnDat' is TRUE.
#'@export

barGene<-function(genes, counts, conditions, title="Gene expression", norm=NULL, eb="sd", returnDat=F, col="Dark2", ord=NULL){
  ylab="Mean"
  #norm is the condition to normalize expression to for relative expression
  counts<-counts[which(rownames(counts) %in% genes),]
  #set up conditions
  condition<-as.factor(c(rep(conditions, each=nrow(counts))))
  genes2<-as.factor(rep(rownames(counts), length(conditions)))
  y<- counts %>% tidyr::gather(key="Sample", value="expression") %>% dplyr::mutate(group=condition) %>%
    dplyr::mutate(gene=genes2)

  #Summarize the data using the function above (get mean and sd/se for the bar plots)
  x<-data_sum(y, eb)
  x<- x %>% dplyr::mutate(gene=forcats::fct_relevel(gene, genes))
  if(!is.null(norm)){
    message("Normalizing values to ", norm, "...")
    y<-data.frame()
    for(i in 1:length(levels(x$gene))){
      tmp<-subset(x, gene==levels(x$gene)[i])
      tmp$mean<-tmp$mean/(tmp$mean[which(tmp$group %in% norm)])

      y<-rbind(y, tmp)
    }
    x<- y %>% dplyr::group_by(group)
    group_levels<-c(norm, levels(condition)[which(!levels(condition) %in% norm)])
    x<- x%>% dplyr::mutate(group=forcats::fct_relevel(group, group_levels))
    ylab="Relative"
  }
  if(!is.null(ord)){
    if(length(ord) == length(levels(factor(conditions)))){
      x<- x%>% dplyr::mutate(group=forcats::fct_relevel(group, ord))
    } else {
      message("length(ord) is not equal to length(levels(factor(conditions))). Ignoring...")
    }
  }
  #And now, we're ready for plotting
  p<-ggplot2::ggplot(x, ggplot2::aes(x=gene, y=mean, fill=group)) +
    ggplot2::geom_bar(stat="identity", color="black", position=ggplot2::position_dodge()) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=mean-eb, ymax=mean+eb), width=.2, position=ggplot2::position_dodge(.9)) +
    ggplot2::labs(title=title, y=paste0(ylab," Normalized Expression (+/-",eb,")"), x="Gene") +
    ggplot2::theme_minimal() + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60, hjust=1)) +
    ggplot2::scale_fill_manual(values=colPal(col))
  if(!is.null(norm)){
    p<- p + ggplot2::geom_hline(yintercept=1, linetype="dashed", color="black")
  }
  print(p) #Print the plot (was saved to 'p')
  if(isTRUE(returnDat)){
	  return(list(rawData=y, Summary=x)) #Return a list with the raw and summarized data - if you want to invesigate it later
  }
}


