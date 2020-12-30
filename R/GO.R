#' Run g:profiler enrichment analysis
#'
#' This function runs a g:profiler gost enrichment analysis on a query list of genes
#' and export a complete results table, a .gem file for use with EnrichmentMap Cytoscape app,
#' and a pdf of the top enriched and significant terms.
#'
#' @param geneList A character vector of genes to query
#' @param species A character vector of the organism name to use. Concatenate the first letter of the name and family name. Default is 'hsapiens'
#' @param bg A character vector describing the genes making up the background. Default NULL
#' @param source A vector of data sources to use. Currently: GO (GO:BP, GO:MF, GO:CC), KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP.
#' @param corr A character vector describing the correction method to use. One of 'gSCS', 'fdr', or 'bonferroni'. Default is 'fdr'
#' @param iea Boolean values indicating if electronic annotations should be excluded. Default FALSE.
#' @param prefix A character vector describing the path and prefix of the output files (should not include any file extensions)
#' @param ts Number indicating the minimum term size - the minimum number of genes per term when generating the plot of most enriched terms. Default is 10.
#' @return Exports a table of analysis results in 'prefix.GO.txt', a gem file in 'prefix.gem', and a pdf with two figures: top ten enriched terms followed by top ten significant terms in 'prefix.top10.pdf'
#' @export
#'

GO_GEM<-function(geneList,species="hsapiens",bg=NULL,source=NULL, corr="fdr", iea=FALSE, prefix="GO_analysis", ts=10){
  #ts is term size (for plotting, terms must have ts genes to make cutoff, default is 10)
  x<-gprofiler2::gost(geneList, organism=species, custom_bg=bg, sources=source, evcodes=TRUE, multi_query=FALSE, correction_method=corr, exclude_iea=iea)
  y<-x$result[,-14]
  y$enrichment <- (y$intersection_size/y$query_size)/(y$term_size/y$effective_domain_size)
  gem<-x$result[,c("term_id","term_name","p_value","intersection")]
  colnames(gem)<-c("GO.ID", "Description","p.Val","Genes")
  gem$FDR<-gem$p.Val
  gem$Phenotype="+1"
  gem<-gem[,c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")]
  write.table(gem, paste0(prefix,".gem.txt"), quote=FALSE, sep="\t", row.names = FALSE)
  y<-y[order(y$enrichment, decreasing=T),]
  GO_plot(y, prefix, ts)
  write.table(y, paste0(prefix,".GO.txt"), quote=FALSE, sep="\t", row.names=FALSE)
}

#' Generates GO figures
#'
#' Function that generates the top 10 term figures for GO_GEM()

GO_plot<-function(GOres, prefix, ts){
  tmp<-GOres[order(GOres$p_value),] #order by increasing p-value
  #Take the top ten terms (already sorted by enrihcment)
  GOres<-GOres[which(GOres$term_size >= ts),]
  GOres<-GOres[1:10,]
  ylim.prim<-c(0, max(GOres$enrichment)+1)
  ylim.sec<-c(0, max(-log10(GOres$p_value))+1)
  b<-diff(ylim.prim)/diff(ylim.sec)
  a<-b*(ylim.prim[1]=ylim.sec[1])
  p<- ggplot2::ggplot(GOres, ggplot2::aes(x=seq(1:length(term_name)), y=enrichment)) +
    ggplot2::geom_col(fill="blue", width=0.75) + ggplot2::geom_col(aes(x=seq(1:length(term_name)), y=a+(-log10(p_value))*b), fill="orange", width=0.375) +
    ggplot2::scale_x_continuous(name="GO Term", breaks=1:10, labels=GOres$term_name) +
    ggplot2::scale_y_continuous(name="Enrichment", sec.axis=sec_axis(~(. -a)/b, name="-Log10 P-value")) +
    ggplot2::theme_classic() + ggplot2::theme(axis.title.y=ggplot2::element_text(color="blue"), axis.title.y.right=ggplot2::element_text(color="orange"),
                            axis.text.x=ggplot2::element_text(angle=60, hjust=1))
  #Now get the top ten significant terms with less than 500 genes/term (in tmp)
  tmp<-tmp[which(tmp$term_size <=500),]
  tmp<-tmp[1:10,]
  q<- ggplot2::ggplot(tmp, ggplot2::aes(x=seq(1:length(term_name)), y=enrichment)) +
    ggplot2::geom_col(fill="blue", width=0.75) + ggplot2::geom_col(aes(x=seq(1:length(term_name)), y=a+(-log10(p_value))*b), fill="orange", width=0.375) +
    ggplot2::scale_x_continuous(name="GO Term", breaks=1:10, labels=tmp$term_name) +
    ggplot2::scale_y_continuous(name="Enrichment", sec.axis=sec_axis(~(. -a)/b, name="-Log10 P-value")) +
    ggplot2::theme_classic() + ggplot2::theme(axis.title.y=ggplot2::element_text(color="blue"), axis.title.y.right=ggplot2::element_text(color="orange"),
                            axis.text.x=ggplot2::element_text(angle=60, hjust=1))
  pdf(paste0(prefix,".Top10.pdf"))
  print(p)
  print(q)
  dev.off()
}

