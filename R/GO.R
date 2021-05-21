#' Run g:profiler enrichment analysis
#'
#' This function runs a g:profiler gost enrichment analysis on a query list of genes
#' and export a complete results table, a .gem file for use with EnrichmentMap Cytoscape app,
#' and a pdf of the top enriched and significant terms.
#'
#' @param geneList A character vector or named list of genes to query. If a named list is provided, the name of each list element is appended to the prefix for GO output.
#' @param species A character vector of the organism name to use. Concatenate the first letter of the name and family name. Default is 'hsapiens'
#' @param bg A character vector describing the genes making up the background. Default NULL
#' @param source A vector of data sources to use. Currently: GO (GO:BP, GO:MF, GO:CC), KEGG, REAC, TF, MIRNA, CORUM, HP, HPA, WP.
#' @param corr A character vector describing the correction method to use. One of 'gSCS', 'fdr', or 'bonferroni'. Default is 'fdr'
#' @param iea Boolean values indicating if electronic annotations should be excluded. Default FALSE.
#' @param prefix A character vector describing the path and prefix of the output files (should not include any file extensions)
#' @param ts Vector of length 2 indicating the minimum and maximum term size - the minimum/maximum number of genes per term when generating the plot of most enriched/significant terms. Default is c(10,500).
#' @param pdf Boolean indicating if bar plots should be exported to pdf. Default is TRUE.
#' @param fig Boolean indicating if bar plots should be printed to R output. Default is TRUE.
#' @param returnGost Boolean indicating if gost results should be returned for future use with gprofiler2 functions. Default is FALSE.
#' @return Exports a table of analysis results in 'prefix.GO.txt', a gem file in 'prefix.gem', and a pdf with two figures: top ten enriched terms followed by top ten significant terms in 'prefix.top10.pdf'
#' @export
#'

GO_GEM<-function(geneList,species="hsapiens",bg=NULL,source=NULL, corr="fdr", iea=FALSE, prefix="GO_analysis", ts=c(10,500),
                 pdf=T, fig=T, returnGost=F){
  GOfun<-function(genes, spec=species, cbg=bg, dsource=source, corrm=corr, exiea=iea, prefix=pre, termsz=ts, prpdf=pdf, prfig=fig, giveGost=returnGost){
    #ts is term size (for plotting, terms must have ts genes to make cutoff, default is 10)
    x<-gprofiler2::gost(genes, organism=spec, custom_bg=cbg, sources=dsource, evcodes=TRUE, multi_query=FALSE, correction_method=corrm, exclude_iea=exiea)
    y<-x$result[,-14]
    y$enrichment <- (y$intersection_size/y$query_size)/(y$term_size/y$effective_domain_size)
    gem<-x$result[,c("term_id","term_name","p_value","intersection")]
    colnames(gem)<-c("GO.ID", "Description","p.Val","Genes")
    gem$FDR<-gem$p.Val
    gem$Phenotype="+1"
    gem<-gem[,c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")]
    write.table(gem, paste0(prefix,".gem.txt"), quote=FALSE, sep="\t", row.names = FALSE)
    y<-y[order(y$enrichment, decreasing=T),]
    GO_plot(y, prefix, termsz, prpdf, prfig)
    write.table(y, paste0(prefix,".GO.txt"), quote=FALSE, sep="\t", row.names=FALSE)
    if(isTRUE(giveGost)){
      return(x)
    }
  }
  if(is.list(geneList)){
    message("geneList is list of character vectors...")
    if(isTRUE(returnGost)){
      gostList<-list()
      for(i in 1:length(geneList)){
        pre<-paste0(prefix,names(geneList)[i])
        message("Analyzing ",names(geneList)[i]," and saving as: ",pre,"...")
        gostList[[i]]<-tryCatch({GOfun(geneList[[i]], prefix=pre)},
                                error=function(e){
                                  message("No significant results for ",names(geneList)[i]," returning NA.")
                                  return(NA)
                                })
        names(gostList)[i]<-names(geneList)[i]
      }
      return(gostList)
    } else {
      for(i in 1:length(geneList)){
        pre<-paste0(prefix,names(geneList)[i])
        message("Analyzing ",names(geneList)[i]," and saving as: ",pre,"...")
        tryCatch({GOfun(geneList[[i]], prefix=pre)},
                 error=function(e){
                   message("No significant results for ",names(geneList)[i],"...")
                 })
      }
    }
  } else {
    message("Single character list of genes provided...")
    pre<-prefix
    if(isTRUE(returnGost)){
      x<-GOfun(geneList, prefix=pre)
      return(x)
    } else {
      GOfun(geneList, prefix=pre)
    }
  }
}

#' Generates GO figures
#'
#' Function that generates the top 10 term figures for GO_GEM()

GO_plot<-function(GOres, prefix, ts, pdf, fig){
  tmp<-GOres[order(GOres$p_value),] #order by increasing p-value
  #Take the top ten terms (already sorted by enrihcment)
  GOres<-GOres[which(GOres$term_size >= ts[1]),]
  GOres<-GOres[1:10,]
  ylim.prim<-c(0, max(GOres$enrichment)+1)
  ylim.sec<-c(0, max(-log10(GOres$p_value))+1)
  b<-diff(ylim.prim)/diff(ylim.sec)
  a<-b*(ylim.prim[1]=ylim.sec[1])
  p<- ggplot2::ggplot(GOres, ggplot2::aes(x=seq(1:length(term_name)), y=enrichment)) +
    ggplot2::geom_col(fill="blue", width=0.75) + ggplot2::geom_col(ggplot2::aes(x=seq(1:length(term_name)), y=a+(-log10(p_value))*b), fill="orange", width=0.375) +
    ggplot2::scale_x_continuous(name="GO Term", breaks=1:10, labels=GOres$term_name) +
    ggplot2::scale_y_continuous(name="Enrichment", sec.axis=ggplot2::sec_axis(~(. -a)/b, name="-Log10 P-value")) +
    ggplot2::theme_classic() + ggplot2::theme(axis.title.y=ggplot2::element_text(color="blue"), axis.title.y.right=ggplot2::element_text(color="orange"),
                            axis.text.x=ggplot2::element_text(angle=60, hjust=1)) +
    ggplot2::labs(title=paste0("Top Ten Enriched terms (>",ts[1],"genes/term)"),
                  subtitle=prefix)
  #Now get the top ten significant terms with less than 500 genes/term (in tmp)
  tmp<-tmp[which(tmp$term_size <=ts[2]),]
  tmp<-tmp[1:10,]
  q<- ggplot2::ggplot(tmp, ggplot2::aes(x=seq(1:length(term_name)), y=enrichment)) +
    ggplot2::geom_col(fill="blue", width=0.75) + ggplot2::geom_col(ggplot2::aes(x=seq(1:length(term_name)), y=a+(-log10(p_value))*b), fill="orange", width=0.375) +
    ggplot2::scale_x_continuous(name="GO Term", breaks=1:10, labels=tmp$term_name) +
    ggplot2::scale_y_continuous(name="Enrichment", sec.axis=ggplot2::sec_axis(~(. -a)/b, name="-Log10 P-value")) +
    ggplot2::theme_classic() + ggplot2::theme(axis.title.y=ggplot2::element_text(color="blue"), axis.title.y.right=ggplot2::element_text(color="orange"),
                            axis.text.x=ggplot2::element_text(angle=60, hjust=1))+
    ggplot2::labs(title=paste0("Top Ten Significant terms (<",ts[2],"genes/term)"),
                  subtitle=prefix)
  if(isTRUE(pdf)){
    pdf(paste0(prefix,".Top10.pdf"))
    print(p)
    print(q)
    dev.off()
  }
  if(isTRUE(fig)){
    print(p)
    print(q)
  }
}

