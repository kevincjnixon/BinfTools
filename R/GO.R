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
#' @param ord Boolean indicating if ranked GO analysis should be performed. Default=FALSE.
#' @param prefix A character vector describing the path and prefix of the output files (should not include any file extensions)
#' @param ts Vector of length 2 indicating the minimum and maximum term size - the minimum/maximum number of genes per term when generating the plot of most enriched/significant terms. Default is c(10,500).
#' @param pdf Boolean indicating if bar plots should be exported to pdf. Default is TRUE.
#' @param fig Boolean indicating if bar plots should be printed to R output. Default is TRUE.
#' @param figCols Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used for enrichment and significance, respectively in the output figures. Default is c("blue","orange").
#' @param returnGost Boolean indicating if gost results should be returned for future use with gprofiler2 functions. Default is FALSE.
#' @param writeRes Boolean indicating if GO.txt results should be written to file 'prefix.GO.txt'
#' @param writeGem Boolean indicating if gem.txt results should be written to file.
#' @param writeGene Boolean indicating if gene.txt results should be written to file.
#' @return Exports a table of analysis results in 'prefix.GO.txt', a gem file in 'prefix.gem', and a pdf with two figures: top ten enriched terms followed by top ten significant terms in 'prefix.top10.pdf'
#' @export
#'

GO_GEM<-function(geneList,species="hsapiens",bg=NULL,source=NULL, corr="fdr", iea=FALSE, ord=FALSE, prefix="GO_analysis", ts=c(10,500),
                 pdf=TRUE, fig=TRUE, figCols=c("blue","orange"), returnGost=FALSE, writeRes=TRUE, writeGem=FALSE, writeGene=FALSE, returnRes=FALSE){
  GOfun<-function(genes, spec=species, cbg=bg, dsource=source, corrm=corr, exiea=iea, rnk=ord, prefix=pre, termsz=ts, prpdf=pdf, prfig=fig, cols=figCols, giveGost=returnGost,
                  gemWrite=writeGem, resWrite=writeRes, genWrite=writeGene, giveRes=returnRes){
    #ts is term size (for plotting, terms must have ts genes to make cutoff, default is 10)
    if(isTRUE(genWrite)){
      write.table(genes, paste0(prefix,".genes.txt"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    }
    x<-gprofiler2::gost(genes, organism=spec, ordered_query=rnk, custom_bg=cbg, sources=dsource, evcodes=TRUE, multi_query=FALSE, correction_method=corrm, exclude_iea=exiea)
    y<-x$result[,-14]
    y$enrichment <- (y$intersection_size/y$query_size)/(y$term_size/y$effective_domain_size)
    gem<-x$result[,c("term_id","term_name","p_value","intersection")]
    colnames(gem)<-c("GO.ID", "Description","p.Val","Genes")
    gem$FDR<-gem$p.Val
    gem$Phenotype="+1"
    gem<-gem[,c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")]
    if(isTRUE(gemWrite)){
      write.table(gem, paste0(prefix,".gem.txt"), quote=FALSE, sep="\t", row.names = FALSE)
    }
    y<-y[order(y$enrichment, decreasing=TRUE),]
    GO_plot(y, prefix, termsz, prpdf, prfig, cols)
    if(isTRUE(resWrite)){
      write.table(y, paste0(prefix,".GO.txt"), quote=FALSE, sep="\t", row.names=FALSE)
    }
    if(isTRUE(giveGost)){
      return(x)
    } else {
      if(isTRUE(giveRes)){
        return(y)
      }
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
      if(isTRUE(returnRes)) {
        resList<-list()
        for(i in 1:length(geneList)){
          pre<-paste0(prefix,names(geneList)[i])
          message("Analyzing ",names(geneList)[i]," and saving as: ",pre,"...")
          resList[[i]]<-tryCatch({GOfun(geneList[[i]], prefix=pre)},
                                 error=function(e){
                                   message("No significant results for ",names(geneList)[i],"...")
                                   return(NA)
                                 })
          names(resList)[i]<-names(geneList)[i]
        }
        return(resList)
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
    }
  } else {
    message("Single character list of genes provided...")
    pre<-prefix
    if(isTRUE(returnGost)){
      x<-GOfun(geneList, prefix=pre)
      return(x)
    } else {
      if(isTRUE(returnRes)) {
        x<- GOfun(geneList, prefix=pre)
        return(x)
      } else {
        GOfun(geneList, prefix=pre)
      }
    }
  }
}

#' Generates GO figures
#'
#' Function that generates the top 10 term figures for GO_GEM()

GO_plot<-function(GOres, prefix, ts, pdf, fig, col, print=c("both","sig","enr")){
  tmp<-GOres[order(GOres$p_value),] #order by increasing p-value
  #Take the top ten terms (already sorted by enrihcment)
  GOres<-GOres[which(GOres$term_size >= ts[1]),]
  GOres<-GOres[1:10,]
  ylim.prim<-c(0, max(GOres$enrichment)+1)
  ylim.sec<-c(0, max(-log10(GOres$p_value))+1)
  b<-diff(ylim.prim)/diff(ylim.sec)
  a<-b*(ylim.prim[1]=ylim.sec[1])
  subt<-strsplit(prefix, split="/", fixed=TRUE)[[1]][length(strsplit(prefix, split="/", fixed=TRUE)[[1]])]
  p<- ggplot2::ggplot(GOres, ggplot2::aes(x=seq(1:length(term_name)), y=enrichment)) +
    ggplot2::geom_col(fill=colPal(col)[1], width=0.75) + ggplot2::geom_col(ggplot2::aes(x=seq(1:length(term_name)), y=a+(-log10(p_value))*b), fill=colPal(col)[2], width=0.375) +
    ggplot2::scale_x_continuous(name="GO Term", breaks=1:10, labels=GOres$term_name) +
    ggplot2::scale_y_continuous(name="Enrichment", sec.axis=ggplot2::sec_axis(~(. -a)/b, name="-Log10 P-value")) +
    ggplot2::theme_classic() + ggplot2::theme(axis.title.y=ggplot2::element_text(color=colPal(col)[1]), axis.title.y.right=ggplot2::element_text(color=colPal(col)[2]),
                                              axis.text.x=ggplot2::element_text(angle=60, hjust=1)) +
    ggplot2::labs(title=paste0("Top Ten Enriched terms (>",ts[1],"genes/term)"),
                  subtitle=subt)
  #Now get the top ten significant terms with less than 500 genes/term (in tmp)
  tmp<-tmp[which(tmp$term_size <=ts[2]),]
  tmp<-tmp[1:10,]
  q<- ggplot2::ggplot(tmp, ggplot2::aes(x=seq(1:length(term_name)), y=enrichment)) +
    ggplot2::geom_col(fill=colPal(col)[1], width=0.75) + ggplot2::geom_col(ggplot2::aes(x=seq(1:length(term_name)), y=a+(-log10(p_value))*b), fill=colPal(col)[2], width=0.375) +
    ggplot2::scale_x_continuous(name="GO Term", breaks=1:10, labels=tmp$term_name) +
    ggplot2::scale_y_continuous(name="Enrichment", sec.axis=ggplot2::sec_axis(~(. -a)/b, name="-Log10 P-value")) +
    ggplot2::theme_classic() + ggplot2::theme(axis.title.y=ggplot2::element_text(color=colPal(col)[1]), axis.title.y.right=ggplot2::element_text(color=colPal(col)[2]),
                                              axis.text.x=ggplot2::element_text(angle=60, hjust=1))+
    ggplot2::labs(title=paste0("Top Ten Significant terms (<",ts[2],"genes/term)"),
                  subtitle=subt)
  if(isTRUE(pdf)){
    pdf(paste0(prefix,".Top10.pdf"))
    if(print[1] == "both"){
      print(p)
      print(q)
    }
    if(print[1] == "sig"){
      print(q)
    }
    if(print[1] == "enr"){
      print(p)
    }
    dev.off()
  }
  if(isTRUE(fig)){
    if(print[1] == "both"){
      print(p)
      print(q)
    }
    if(print[1] == "sig"){
      print(q)
    }
    if(print[1] == "enr"){
      print(p)
    }
  }
}

#'Combine GO Results from DE analysis into one plot
#'
#'This function will combine GO results from GO_GEM when performed on DEGs into a single plot.
#'
#'@param GOresList list of GO_GEM results (set returnRes=TRUE) of length 2. Best for DEG GO results, where the first entry in the list is downregulated, and the second entry in the list is upregulated
#'@param title Character for the plot title
#'@param ts numeric vector of length 2 with the minimum and maximum term sizes for plotting. Default is c(10,500)
#'@param sig Boolean indicating if top significant results should be plotted. Set to FALSE to plot top enriched results
#'@param numTerm Numeric indicating the number of top terms to plot. Default is 10.
#'@param upcols character vector of length 2 indicating the colour for upregulated enrichment and significance bars
#'@param downcols character vector of length 2 indicating the colour for downregulated enrichment and significance bars
#'@param labs character vector of length 2 indicating the specific legend labels for the entries in GOresList. Default=c("Downregulated","Upgreulated")
#'@param textsize numeric indicating text size for plot. Leave NULL to keep default size.
#'@param retGP Booliean indicating if the ggplot2 object should be returned for further editing. Default=FALSE.
#'@return A plot of top GO results
#'@export

combGO_plot<-function(GOresList, title="GO results", ts=c(10,500), sig=TRUE, numTerm=10, upcols=c("lightpink","red"),
                      downcols=c("lightblue","blue"), labs=c("Downregulated","Upregulated"), textsize=NULL, retGP=FALSE){
  if(length(GOresList) != 2){
    stop("length(GOresList) must be 2!")
  }
  #print(names(GOresList))
  #message("Down should be first. Up second.")
  GOresList[[1]]$dir<-factor(rep(labs[1], nrow(GOresList[[1]])))
  GOresList[[1]]$enrichment<-GOresList[[1]]$enrichment*-1
  GOresList[[1]]$sig<-(-log10(GOresList[[1]]$p_value))*-1
  GOresList[[2]]$dir<-factor(rep(labs[2], nrow(GOresList[[2]])))
  GOresList[[2]]$sig<-(-log10(GOresList[[2]]$p_value))
  filtRes<-NULL
  if(isTRUE(sig)){
    GOresList[[1]]<-GOresList[[1]][order(GOresList[[1]]$p_value),]
    GOresList[[2]]<-GOresList[[2]][order(GOresList[[2]]$p_value),]
    GOresList[[1]]<-subset(GOresList[[1]], term_size <= ts[2])
    GOresList[[2]]<-subset(GOresList[[2]], term_size <= ts[2])
    filtRes<-rbind(head(GOresList[[1]], n=numTerm),
                   head(GOresList[[2]], n=numTerm)[order(head(GOresList[[2]], n=numTerm)$p_value, decreasing=TRUE),])
  } else {
    GOresList[[1]]<-GOresList[[1]][order(GOresList[[1]]$enrichment),]
    GOresList[[2]]<-GOresList[[2]][order(GOresList[[2]]$enrichment, decreasing=TRUE),]
    GOresList[[1]]<-subset(GOresList[[1]], term_size >= ts[1])
    GOresList[[2]]<-subset(GOresList[[2]], term_size >= ts[1])
    filtRes<-rbind(head(GOresList[[1]], n=numTerm),
                   head(GOresList[[2]], n=numTerm)[order(head(GOresList[[2]], n=numTerm)$enrichment),])
  }
  #print(head(GOresList[[1]]$term_name))


  ylim.prim<-c(0, max(filtRes$enrichment)+1)
  ylim.sec<-c(0, max(-log10(filtRes$p_value))+1)
  b<-diff(ylim.prim)/diff(ylim.sec)
  a<-b*(ylim.prim[1]=ylim.sec[1])

  # upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$dir == "Up"))
  # downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$dir == "Down"))
  #upcols<-c("red","orange")
  #downcols<-c("darkblue", "purple")
  colos = c(upcols, downcols)
  #names(colos) = 1:length(colos)
  names(colos)=c("Up_Enrichment","Up_Significance","Down_Enrichment", "Down_Significance")
  #filtRes$Index = as.factor(1:nrow(filtRes))
  filtRes$fill_E=c(rep(downcols[1], length(filtRes$dir[which(filtRes$dir %in% labs[1])])),
                   rep(upcols[1], length(filtRes$dir[which(filtRes$dir %in% labs[2])])))
  filtRes$Enrichment=filtRes$dir
  filtRes$fill_S=c(rep(downcols[2], length(filtRes$dir[which(filtRes$dir %in% labs[1])])),
                   rep(upcols[2], length(filtRes$dir[which(filtRes$dir %in% labs[2])])))
  filtRes$Significance=filtRes$dir
  #print(head(filtRes))
  #print(tail(filtRes))
  g = ggplot2::ggplot(filtRes, ggplot2::aes(x=seq(1:length(term_name)), y=enrichment)) +
    ggplot2::geom_col( ggplot2::aes(fill = dir) , width=0.9) +
    ggplot2::scale_fill_manual(name="Enrichment", values = c(downcols[1], upcols[1])) +
    ggplot2::geom_col(ggplot2::aes(x=seq(1:length(term_name)), y=a+sig*b, color=Significance), fill=filtRes$fill_S, width=0.375) +
    ggplot2::scale_colour_manual(values=c(downcols[2], upcols[2])) +
    ggplot2::guides(colour=ggplot2::guide_legend(override.aes = list(fill=c(downcols[2],upcols[2]), colour=c(downcols[2], upcols[2]), name="Significance"))) +
    #ggplot2::scale_fill_manual(name="Significance", values=c(downcols[2], upcols[2])) +
    ggplot2::scale_x_continuous(name="GO Term", breaks=1:length(filtRes$term_name), labels=filtRes$term_name) +
    ggplot2::scale_y_continuous(name="Enrichment", sec.axis=ggplot2::sec_axis(~(. -a)/b, name="-Log10 P-value")) +
    ggplot2::labs(title=title) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position="right", text=ggplot2::element_text(size=textsize))
  if(isTRUE(retGP)){
    return(g)
  } else {
    print(g)
  }
  #return(filtRes)
}

#'Make a heatmap showing significance of groups of GO terms from multiple sets of results
#'
#'@param GOresList list of GO_GEM results (set returnRes=TRUE) of length >= 2.
#'@param termList named list of character vectors where the names represent the group names for GO terms and vectors contain GO term names corresponding to the 'term_name' column in the GO_GEM results data frams.
#'@param hmcol colorRampPalette of length 100 that will direct the colour palette of the heatmap. Default is colorRampPalette(c("white","darkblue"))(100)
#'@param width numeric indicating the cell width (Default=NA will automatically direct the cell width)
#'@param height numeric indicating the cell height (Default=NA will automatically direct the cell height)
#'@param maxVal numeric indicating the maximum -log10 p-value for the heatmap. Default NA sets it automatically.
#'@param minVal numeric indicating the minimum -log10 p-value (1.3 should be absolute mimum - pvalue 0.05). Can only be set if maxVal is not NA. Default NA sets it automatically.
#'@param NAcol character indicating the colour of NA values in heatmap. Default is "darkgrey"
#'@param ret Boolean indicating if table used for the heatmap (p-values not in -log10) should be returned. Default=FALSE.
#'@return A grouped heatmap showing significance of GO terms across analyses as -log10(p-value)
#'@export

GOHeat<-function(GOresList, termList, hmcol=colorRampPalette(c("white","darkblue"))(100), width=NA, height=NA, maxVal=NA, minVal=NA, NAcol="darkgrey", ret=FALSE){
  if(any(duplicated(unlist(termList)))){
    message("Duplicate term names found in termList. Removing all but the first instance of the term...")
    #creates a named vector where the names give the location of the duplicated instance of the term
    x<-data.frame(dups=names(unlist(termList)[which(duplicated(unlist(termList)))]),
                  term=unlist(termList)[which(duplicated(unlist(termList)))])
    #separate the column 'dups' into the names and indices for the duplicated instances
    loc<-tidyr::separate(x, col=dups, into=c("name","index"), sep="(?<=[A-Za-z])(?=[0-9])")
    #Now we remove them
    for(i in unique(loc$name)){
      tmp<-subset(loc[which(loc$name %in% i),])
      termList[[which(names(termList) %in% i)]]<-termList[[which(names(termList) %in% i)]][-as.numeric(tmp$index)]
    }
  }
  retP<-function(GOres, term){
    p_val<-GOres$p_value[which(GOres$term_name %in% term)]
    if(length(p_val)<1){
      p_val<-NA
    }
    return(p_val)
  }
  forHeat<-suppressWarnings(unique(tidyr::gather(as.data.frame(do.call("cbind", termList)), key="Group",value="Term")))
  tmp<-c()
  for(i in forHeat$Term){
    tmp<-rbind(tmp, unlist(lapply(GOresList, retP, term=i)))
  }
  colnames(tmp)<-names(GOresList)
  forHeat<-cbind(forHeat, tmp)

  gaps<-c()
  if(length(termList)<2){
    gaps<-NULL
  } else {
    for(i in 1:(length(termList)-1)){
      gaps<-c(gaps, sum(gaps[length(gaps)],length(termList[[i]])))
    }
  }
  #print(head(forHeat))
  tmp <- -log10(forHeat[,-c(1:2)])
  rownames(tmp)<-forHeat$Term
  rowAnno<-data.frame(row.names=forHeat$Term, Group=forHeat$Group)
  rowAnno$Group<-suppressWarnings(as.factor(rowAnno$Group))
  rowAnno$Group<-suppressWarnings(forcats::fct_relevel(rowAnno$Group, names(termList)))
  if(!is.na(maxVal)){
    if(!is.na(minVal)){
      maxVal<-seq(from=minVal, to=maxVal, length.out=100)
    } else {
      maxVal<-seq(from=min(as.matrix(tmp)[is.finite(as.matrix(tmp))]), to=maxVal, length.out=100)
    }
  }
  pheatmap::pheatmap(tmp, scale="none", cluster_rows=FALSE, cluster_cols=FALSE, legend=TRUE, annotation_row=rowAnno,
                     gaps_row=gaps, color=hmcol, border_color="black",
                     cellwidth = width, cellheight = height, breaks=maxVal, na_col=NAcol)
  if(isTRUE(ret)){
    return(forHeat)
  }
}

#'Custom GO Function
#'
#'Perform a hypergeometric enrichment test for a custom gene set
#'
#' @param genes character vector of genes to be queried
#' @param gmt Named list of gene set(s) to be analyzed for enrichment
#' @param gsName Character specifying geneSet source - to be incorportated into output and help indicate source of gene sets
#' @param bg character vector of genes indicating the background of genes. If left NULL, argument 'sp' will indicate species and assume all genes are in background.
#' @param sp character of either "human" or "mouse" indicating the species. This will provide the number of background genes if bg=NULL.
#' @param FDR Boolean indicating if p-values should be FDR corrected (Benjamini-Hochberg). Default is TRUE.
#' @param byRegion Boolean. If query genes has duplicate gene symbols, set to TRUE to remove duplicates (relevant only if query genes are from annotated peaks from genomic data). Default=FALSE.
#' @param enr character of either "pos" or "neg" indicating if significance should be calculated for positive or negative hypergeometric enrichment, respectively. default="pos".
#' @param significant Boolean indicating if only significant (p<0.05) results should be returned. defaulte=TRUE.
#' @param minp numeric indicating the minimum p-value possible reported (to avoid zeros)
#' @return Data.frame of same structure of results table from GO_GEM() when returnRes=TRUE.
#' @export

customGO<-function(genes, gmt, gsName="custom GeneSet", bg=NULL, sp="human", FDR=TRUE, byRegion=FALSE, enr="pos", significant=TRUE, minp=1e-300){
  if(!is.list(gmt)){
    gmt<-BinfTools::read.gmt(gmt)
  }
  filtBG<-function(x, bg){
    return(x[which(x %in% bg)])
  }
  if(!is.null(bg)){
    message("Filtering gene sets to correspond with provided background")
    gmt<-lapply(gmt, filtBG, bg)
  } else {
    if(is.numeric(sp)){
      message("No custom background, using provided background size in 'sp' argument: ", sp)
      bg<-sp
    } else {
      message("No custom background, using genome size for ",sp,"...")
      if(sp == "human"){
        bg<-c(1:18123) #Number used from g:profiler
      } else {
        if(sp == "mouse"){
          bg<-c(1:18172) #Number from g:profiler
        } else {
          message("Using gmt as background...")
          bg<-length(unique(unlist(gmt)))
        }
      }
    }
  }
  findOL<-function(gs, x, retVal=FALSE, unique){
    if(isTRUE(unique)){
      x<-unique(x)
    }
    if(isFALSE(retVal)){
      return(length(x[which(x%in% gs)]))
    } else{
      return(x[which(x %in% gs)])
    }
  }
  OL<-unlist(lapply(gmt, findOL, genes, unique=byRegion))
  OL2<-lapply(gmt, findOL, genes, retVal=TRUE, unique=byRegion)
  res<-data.frame(query=rep("query_1", length(gmt)),
                  significant=rep("FALSE", length(gmt)),
                  p_value=rep(1, length(gmt)),
                  term_size=lengths(gmt),
                  query_size=rep(length(genes), length(gmt)),
                  intersection_size=OL,
                  precision=rep(0, length(gmt)),
                  recall=rep(0, length(gmt)),
                  #term_id=names(gmt),
                  source=rep(gsName, length(gmt)),
                  term_name=names(gmt),
                  effective_domain_size=length(bg),
                  intersection=unlist(lapply(OL2, toString)),
                  enrichment=rep(0, length(gmt)))
  vectorPhyper<-function(intersection_size, term_size, domain_size,query_size, adj=TRUE, lt=FALSE){
    pvals<-c()
    #pb<-txtProgressBar(0,1,style=3)
    for(i in 1:length(intersection_size)){
      pvals<-c(pvals, phyper(intersection_size[i], term_size[i], domain_size[i]-term_size[i],query_size[i], lower.tail = lt))
      #setTxtProgressBar(pb, i)
    }
    if(isTRUE(adj)){
      pvals<-p.adjust(pvals, "BH")
    }
    return(pvals)
  }
  if(enr == "pos"){
    enr<-FALSE
  } else {
    enr<-TRUE
  }
  print(dim(res))
  res$p_value<-vectorPhyper(res$intersection_size, res$term_size, res$effective_domain_size, res$query_size, adj=FDR, lt=enr)
  res$significant<-ifelse(res$p_value<0.05, TRUE, FALSE)
  res$precision<-res$intersection_size/res$query_size
  res$recall<-res$intersection_size/res$term_size
  res$enrichment<-(res$intersection_size/res$query_size)/(res$term_size/res$effective_domain_size)
  res<-res[order(res$enrichment, decreasing=TRUE),]
  if(!is.null(minp)){
    res$p_value[res$p_value<minp]<-minp
  }
  if(isTRUE(significant)){
    if(nrow(subset(res, significant == "TRUE"))<1){
      message("No significant results. Returning NA.")
      return(NA)
    } else {
      return(subset(res, significant == "TRUE"))
    }
  } else {
    return(res)
  }

}

#Grep a certain key from a GO results table's "term_name" column, and return provided columns
#Great for making a termList for GOHeat (return "term_name" only) or getting intersection genes
#' Grep a key word from a GO results table and return specific columns
#'
#' @param GOtable Data frame returned from GO_GEM or clusGO with returnRes=TRUE
#' @param key Character to be matched in the column 'term_name'. Note that the key is not case sensitive and will return parital matches. Be sure to review the results.
#' @param cols Character vector of any of colnames(GOtable) indicating which columns should be returned. Default is c("term_name","intersection").
#' @return Data frame of subset GOtable where key is found in term_name and columns match cols. If no matches are found, NA is returned.
#' @export

GOgrep<-function(GOtable, key, cols=c("term_name","intersection")){
  x<-GOtable[grep(key, GOtable$term_name, ignore.case=TRUE),which(colnames(GOtable) %in% cols), drop=FALSE]
  if(nrow(x)<1){
    x<-NA
  }
  return(x)
}

#' Remove NAs
#'
#' @param x vector of any type
#' @return x with NAs removed
#' @export
rmNA<-function(x){
  return(x[!is.na(x)])
}

#Filter a list of GO results to either remove tables with no significant results (replace=FALSE), or replace them with a dummy table (replace=TRUE).
#Replacing with a dummy table is good when you want to show all GO analysis names in a GOHeat analysis.
#' Filter a list of GO results to remove/replace tables with no enriched terms
#'
#' @param x list of GO results tables returned when returnRes=TRUE in GO_GEM() or clusGO()
#' @param replace Boolean indicating if slots with no enriched GO terms should be replaced with a 'dummy table' - this is good if you want to run a GOHeat analysis. Default is FALSE.
#' @param sources Character vector indicating which values should be placed in the 'sources' column if replace = TRUE. Defaults to c("GO:BP","GO:MF","GO:CC","CORUM","REAC","WP","HP","TF","KEGG")
#' @return list of GO results tables with tables containing no enriched terms being removed or replaced.
#' @export

filtGO<-function (x, replace = FALSE, sources=c("GO:BP","GO:MF","GO:CC","CORUM","REAC","WP","HP","TF","KEGG"))
{
  to_rm <- c()
  for (i in 1:length(x)) {
    if (is.null(nrow(x[[i]]))) {
      to_rm <- c(to_rm, i)
    }
  }
  if (length(to_rm) > 0) {
    if (isFALSE(replace)) {
      message("Removing ", length(to_rm), " GO tables with no results.")
      x <- x[-to_rm]
    }
    else {
      message("Replacing ", length(to_rm), " GO tables with no results with a dummy table.")
      tmp <- data.frame(p_value = rep(1, length(sources)), term_name = rep("No Enriched Terms",
                                                                           length(sources)), enrichment = rep(0, length(sources)), precision = rep(0, length(sources), source=sources))
      for (i in 1:length(to_rm)) {
        x[[to_rm[i]]] <- tmp
      }
    }
  }
  return(x)
}

#' Easily make a termlist for GOHeat()
#'
#' @param GOres List of GO results tables from GO_GEM() with returnRes=TRUE.
#' @param keylist named list of character vectors of keys to match in the term_name column of the GO results
#' @param retCol character or character vectors of column names from GOres to be returned. Default is "term_name".
#' @return named list of term names enriched in at least one of the GO tables in GOres for use with GOHeat(). Be sure to double check the termlist before making the heatmap as unwanted terms could be returned.
#' @export
makeTermList<-function(GOres, keylist, retCol="term_name"){
  termlist<-list()
  for(i in 1:length(keylist)){
    tmp<-c()
    for(k in 1:length(keylist[[i]])){
      tmp<-c(tmp, unlist(lapply(GOres, GOgrep, key=keylist[[i]][k], cols=retCol)))
    }
    termlist[[i]]<-tmp
    termlist[[i]]<-rmNA(unique(termlist[[i]]))
    names(termlist)[i]<-names(keylist)[i]
  }
  return(termlist)
}

#' Create a dot plot of GO results
#'
#' @param GOres Single results table from GO_GEM() output with returnRes=T
#' @param title Character indicating the title of the plot
#' @return Dot plot representing enrichment, precision (as 'Gene Ratio'), and significance of enriched GO terms.
#' @export

GO_Dot<-function(GOres, title=""){
  require(ggplot2, quietly=T)
  x<-data.frame(Enrichment=log2(GOres$enrichment), Gene_Ratio=GOres$precision, FDR=GOres$p_value, term_name=GOres$term_name)
  p<-ggplot(data=x, aes(x=Enrichment, y=term_name, color=FDR, size=Gene_Ratio)) +
    geom_point() + scale_color_gradient(low="black", high="lightblue") + theme_minimal() +
    ylab("") + xlab("log2(Enrichment)") + ggtitle(title) + labs(size="Gene Ratio")
  print(p)
}

#' Dot heatmap for GO results
#'
#' @param GOresList Named list containing GO_GEM() results tables when returnRes=T
#' @param termList Named list of relevant GO term names - made with makeTermList()
#' @param ret Boolean indicating if values used for heatmap should be returned. Default is FALSE.
#' @return Heatmap figure representing enrichment and precision (as 'Gene Ratio') of selected enriched pathways.
#' @export

GOHeat2<-function (GOresList, termList, ret = FALSE)
{
  if (any(duplicated(unlist(termList)))) {
    message("Duplicate term names found in termList. Removing all but the first instance of the term...")
    x <- data.frame(dups = names(unlist(termList)[which(duplicated(unlist(termList)))]),
                    term = unlist(termList)[which(duplicated(unlist(termList)))])
    loc <- tidyr::separate(x, col = dups, into = c("name",
                                                   "index"), sep = "(?<=[A-Za-z])(?=[0-9])")
    for (i in unique(loc$name)) {
      tmp <- subset(loc[which(loc$name %in% i), ])
      termList[[which(names(termList) %in% i)]] <- termList[[which(names(termList) %in%
                                                                     i)]][-as.numeric(tmp$index)]
    }
  }
  retEnr <- function(GOres, term) {
    p_val <- log2(GOres$enrichment[which(GOres$term_name %in% term)])
    if (length(p_val) < 1) {
      p_val <- NA
    }
    return(p_val)
  }
  retGR <- function(GOres, term){
    GR <- GOres$precision[which(GOres$term_name %in% term)]
    if(length(GR) < 1) {
      GR <- NA
    }
    return(GR)
  }

  forHeat <- suppressWarnings(unique(tidyr::gather(as.data.frame(do.call("cbind",
                                                                         termList)), key = "Group", value = "Term")))
  tmp <- c()
  for (i in forHeat$Term) {
    tmp <- rbind(tmp, unlist(lapply(GOresList, retEnr, term = i)))
  }
  colnames(tmp) <- names(GOresList)

  forHeat <- cbind(forHeat, tmp)

  tmp <- c()
  for (i in forHeat$Term) {
    tmp <- rbind(tmp, unlist(lapply(GOresList, retGR, term = i)))
  }
  colnames(tmp) <- names(GOresList)

  forHeat <- cbind(forHeat, tmp)

  #return(forHeat)

  gaps <- c()
  if (length(termList) < 2) {
    gaps <- NULL
  }
  else {
    for (i in 1:(length(termList) - 1)) {
      gaps <- c(gaps, sum(gaps[length(gaps)], length(termList[[i]])))
    }
  }

  d<-data.frame(
    x=rep(names(GOresList), each= nrow(forHeat)),
    y=rep(forHeat$Term, length(GOresList)),
    z=rep(forHeat$Group, length(GOresList)),
    Enrichment=unlist(forHeat[,3:(2+length(GOresList))]),
    Gene_Ratio=unlist(forHeat[,(3+length(GOresList)):(2+2*length(GOresList))]))
  d$y<-factor(d$y, levels=do.call("c",termList))
  d$z<-factor(d$z, levels=names(termList))
  #return(d)

  p<-ggplot(d, aes(x, forcats::fct_rev(y), fill=Enrichment, size=Gene_Ratio)) +
    geom_point(shape=21, stroke=0) + geom_hline(yintercept=seq(0.5,(nrow(forHeat)+0.5),1), size=.2) +
    geom_vline(xintercept=seq(0.5, length(GOresList)+0.5,1), size=.2) +
    scale_x_discrete(position = "top") +
    scale_radius(range=c(1,15)) +
    scale_fill_gradient(low="black", high="lightblue") +
    theme_minimal() +
    theme(legend.position = "bottom",
          panel.grid.major = element_blank(),
          legend.text = element_text(size=8),
          legend.title=element_text(size=8)) +
    guides(size = guide_legend(override.aes = list(fill=NA, color="black", stroke=0.25),
                               label.position="bottom", title.position="right", order=1),
           fill =  guide_colorbar(ticks.colour = NA, title.position = "top", order=2)) +
    labs(size = "Gene Ratio", fill="log2(Enrichment)", x=NULL, y=NULL)
  plot(p)

  p <- ggplot(d, aes(x, forcats::fct_rev(y), fill = Enrichment)) +
    geom_tile(color = "white") +  # Use geom_tile for heatmap squares
    geom_point(aes(size = Gene_Ratio), shape = 16) +  # Add dots for gene ratio
    geom_hline(yintercept = seq(0.5, (nrow(forHeat) + 0.5), 1), size = 0.2) +
    geom_vline(xintercept = seq(0.5, length(GOresList) + 0.5, 1), size = 0.2) +
    scale_x_discrete(position = "top") +
    scale_size(range = c(5, 20)) +
    scale_fill_gradient(low = "white", high = "red", na.value = "white") +  # Adjust color range to red
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 8)
    ) +
    guides(
      size = guide_legend(
        override.aes = list(fill = NA, color = "black", stroke = 0.25),
        label.position = "bottom", title.position = "right", order = 1
      ),
      fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)
    ) +
    labs(size = "Gene Ratio", fill = "log2(Enrichment)", x = NULL, y = NULL)

  #plot(p)
  plot(p + facet_grid(rows=vars(z), scales="free",space="free"))

  if (isTRUE(ret)) {
    return(forHeat)
  }
}
