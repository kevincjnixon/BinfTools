#' Gene Set Enrichment Analysis
#'
#' Perform a GSEA using fgsea and gage analyses
#'
#' @param rnk named vector of gene rankings (names are genes with matching nomenclature of geneset file)
#' @param gmt named list of gene sets or character indicating gmt file path for gene sets to be used in GSEA analysis
#' @param pval Adjusted p-value threshold to filter results. Set to 1 to return all results. Default=1
#' @param ts numeric vector of length two indicating the minimum and maximum gene set sizes
#' @param nperm integer indicating number of permutations to run. Default is 10000
#' @param parseBader Boolean indicating if gene set names should be parsed for the output figure following Bader Lab nomenclature (i.e. '%' as delimiter)
#' @param plot.title Title for the plot
#' @return A data frame of GSEA results and a figure showing normalized enrichment scores (NES) for top positive and negative enriched gene sets.
#' @export

GSEA = function(rnk, gmt, pval=1, ts=c(10,600), nperm=10000, parseBader=T, plot.title = "Gene Set Enrichment Analysis") {
  set.seed(54321)
  require(dplyr, quietly=T)

  if ( any( duplicated(names(rnk)) )  ) {
    warning("Duplicates in gene names")
    rnk = rnk[!duplicated(names(rnk))]
  }
  if  ( !all( order(rnk, decreasing = TRUE) == 1:length(rnk)) ){
    warning("Gene list not sorted")
    rnk = sort(rnk, decreasing = TRUE)
  }
  if(is.character(gmt)){
    if(grep("http", gmt)>0){
      myGO = fgsea::gmtPathways(url(gmt))
   } else {
      myGO = fgsea::gmtPathways(gmt)
   }
  } else {
    myGO = gmt
  }

  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = rnk,
                        minSize=ts[1],
                        maxSize=ts[2],
                        nperm=nperm) %>%
    as.data.frame() %>%
    dplyr::filter(padj < !!pval)
  #print(dim(fgRes))

  ## Filter FGSEA by using gage results. Must be significant and in same direction to keep
  gaRes = gage::gage(rnk, gsets=myGO, same.dir=TRUE, set.size =ts)

  ups = as.data.frame(gaRes$greater) %>%
    tibble::rownames_to_column("Pathway") %>%
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")

  downs = as.data.frame(gaRes$less) %>%
    tibble::rownames_to_column("Pathway") %>%
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")

  #print(dim(rbind(ups,downs)))
  ## Define up / down pathways which are significant in both tests
  keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
  keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]

  fgRes = fgRes[ !is.na(match(fgRes$pathway,
                              c( keepups$pathway, keepdowns$pathway))), ] %>%
    arrange(desc(NES))
  fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")

  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))


  upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
  downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
  colos = c(upcols, downcols)
  names(colos) = 1:length(colos)
  filtRes$Index = as.factor(1:nrow(filtRes))
  if(isTRUE(parseBader)){
    filtRes$pathway<-sapply(strsplit(filtRes$pathway, "%", T),'[[',1)
  }
  g = ggplot2::ggplot(filtRes, ggplot2::aes(reorder(pathway, NES), NES)) +
    ggplot2::geom_col( ggplot2::aes(fill = Index )) +
    ggplot2::scale_fill_manual(values = colos ) +
    ggplot2::coord_flip() +
    ggplot2::labs(x="Pathway", y="Normalized Enrichment Score",
         title=plot.title) +
    ggplot2::theme_minimal() + ggplot2::theme(legend.position="none")#, title=element_text(size=1))

  print(g)
  #output = list("Results" = fgRes, "Plot" = g)
  output<-fgRes
  return(output)
}

plotEnrichment<-function (pathway, stats, gseaParam = 1, ticksSize = 0.2, NES=NULL, title)
{
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  ES<-max(tops)
  #hm<- -diff/1.25
  #axis<-0
  if(NES<0){
    ES<-min(bottoms)
  }
    hm<-min(bottoms)-(diff*1.25)
    axis<-min(bottoms)-diff/2.25
  #}
  hmcol<-colorRampPalette(c("red","grey","blue"))(length(rnk))
  sz<-12
  if(nchar(title)>60){
    sz<-10
  }
  if(nchar(title)>80){
    sz<-8
  }
  #print(paste(title, ":",nchar(title)))
  g <- ggplot2::ggplot(toPlot, ggplot2::aes(x = x, y = y, colour=x)) + ggplot2::geom_point(color = "green", size = 0.1) +
    ggplot2::geom_segment(mapping=ggplot2::aes(x=0, xend=length(rnk), y=ES, yend=ES), colour = "red", linetype = "dashed", size=1) +
    #geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed", size=1) +
    ggplot2::geom_segment(colour = "black", size=1, mapping=ggplot2::aes(x=0, xend=length(rnk), y=axis, yend=axis)) +
    ggplot2::geom_segment(data = data.frame(x=0:length(rnk)), mapping = ggplot2::aes(x=x, y=hm-diff/5, xend=x, yend=hm+diff/5), size=ticksSize) +
    ggplot2::geom_line(color = "green", size=2) + ggplot2::theme_classic() +
    ggplot2::geom_segment(data = data.frame(x = pathway), mapping = ggplot2::aes(x = x, y = axis-diff/2, xend = x, yend = axis+diff/2), size = ticksSize,
                 colour = "black") +
    ggplot2::theme(axis.line.x=ggplot2::element_blank()) +
    ggplot2::labs(x = "Rank", y = "Enrichment Score", title=title) +
    ggplot2::scale_colour_gradientn(colours=hmcol) + ggplot2::theme(legend.position="none", title=ggplot2::element_text(size=sz))
  g
}

#' Generate enrichment plot from GSEA results
#'
#' Generate an enrichment plot for specific GSEA results
#'
#' @param gseaRes Data frame output of GSEA function containing specific row(s) to make enrichment plots
#' @param rnk named vector of gene rankings (names are genes with matching nomenclature of geneset file). Must be same as one used to create gseaRes input
#' @param gmt named list of gene sets or character indicating path to gmt file used to generate gseaRes input
#' @param title character vector of length nrow(gseaRes) with a title for each plot. Leave NULL for automated titles (parsed by '%' or ':")
#' @return Enrichment plot with rank of genes on x-axis and running enrichment score on y-axis
#' @export

enPlot<-function(gseaRes, rnk, gmt, title=NULL){
  if(!is.null(title) && length(title)!=nrow(gseaRes)){
    stop("title must have length of nrow(gseaRes)")
  }
  myGO<-gmt
  if(is.character(gmt)){
     if(grep("http", gmt)>0){
      myGO = fgsea::gmtPathways(url(gmt))
   } else {
    myGO<-fgsea::gmtPathways(gmt)
   }
  }
  #parse through names in gseaRes table, pull the unique identifiers from the names and use it with myGO:
  for(i in 1:nrow(gseaRes)){
    id<-c()
    main<-c()
    if(length(grep(":", gseaRes$pathway[i]))>0){
      #message("Standard names detected...")
      id<-paste0(":",sapply(strsplit(gseaRes$pathway[i],":",T),'[[',2),"$")
      if(is.null(title)){
        main<-paste(gseaRes$pathway[i], "NES:", round(gseaRes$NES[i], digits=3), "padj:", signif(gseaRes$padj[i], digits=3))
      }
    }
    if(length(grep("%", gseaRes$pathway[i]))>0){
      #message("BaderLab names detected....")
      id<-sapply(strsplit(gseaRes$pathway[i],"%",T),'[[',3)
      if(is.null(title)){
        main<-paste(sapply(strsplit(gseaRes$pathway[i],"%",T),'[[',1), "NES:", round(gseaRes$NES[i], digits=3), "padj:", signif(gseaRes$padj[i], digits=3))
      }
    }

    if(!is.null(title)){
      main<-title[i]
    }
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout=grid::grid.layout(2,1, heights=grid::unit(c(0.75,0.25),"npc"))))
    print(plotEnrichment(myGO[[grep(id, names(myGO))[1]]], rnk, NES=gseaRes$NES[i], title=main), vp=grid::viewport(layout.pos.row = 1))

    rnk2<-as.data.frame(rnk, row.names=names(rnk))
    rnk2<-as.data.frame(rnk2[order(rnk2$rnk, decreasing = T),,drop=F])
    g<-ggplot2::ggplot(rnk2, ggplot2::aes(x=seq(1:nrow(rnk2)), y=rnk))+
      ggplot2::geom_bar(stat="identity", fill="lightgrey")+ggplot2::theme_classic() +
      ggplot2::labs(x="Rank", y="Score") +ggplot2::theme(axis.line.x=ggplot2::element_blank())
    print(g, vp=grid::viewport(layout.pos.row = 2))
  }
}

#' Make a custom gmt from GSEA terms
#'
#' @param terms Character vector of GSEA terms (pathway column from results of GSEA() function)
#' @param gmt GMT file name or R object used to generate the gsea results using GSEA(). Can be GSEA results object if setting leadingEdge=T.
#' @param leadingEdge Boolean indicating if the leading edge genes only should be extracted. Default=F
#' @return List gmt object of genesets with terms from GSEA analysis
#' @export
gsea_gmt<-function(terms, gmt, leadingEdge=F){
  if(is.character(gmt)){
     if(grep("http", gmt)>0){
      myGO = fgsea::gmtPathways(url(gmt))
   } else {
    gmt<-fgsea::gmtPathways(gmt)
   }
  }
  pb<-txtProgressBar(min=0, max=length(terms), style=3)
  res<-list()
  index<-1
  for(i in 1:length(terms)){
    if(isFALSE(leadingEdge)){
      tryCatch({res[[index]]<-gmt[[grep(terms[i], names(gmt), ignore.case=T)]]
                names(res)[index]<-terms[i]
      index<-index+1}, error=function(e) NULL)
    } else {
        res[[index]]<-gmt[which(gmt$pathway == terms[i]),]$leadingEdge
      names(res)[index]<-terms[i]
      index<-index+1
    }

    setTxtProgressBar(pb, i)
  }
  return(res)
}


GSEA_EM<-function(gsea, rnk, prefix=NULL, gmt=NULL){
  N=length(rnk)
  df<-data.frame(NAME=gsea$pathway,
                 "GS <br> follow link to MSigDB"=gsea$pathway,
                 "GS DETAILS"=rep("Details...", nrow(gsea)),
                 SIZE=gsea$size,
                 ES=gsea$ES,
                 NES=gsea$NES,
                 "NOM p-val"=gsea$pval,
                 "FDR p-val"=gsea$padj,
                 "FWER p-val"=gsea$padj,
                 "RANK AT MAX"=rep(0, nrow(gsea)),
                 "LEADING EDGE"=rep(NA, nrow(gsea)),
                 check.names=F)
  for(i in 1:nrow(gsea)){
    LE<-lengths(gsea$leadingEdge[i])
    tags<-LE/gsea$size[i]
    genes<-LE/N
    signal<-((tags*(1-genes))*(N/(N-gsea$size[i])))
    df$`LEADING EDGE`[i]<-paste0("tags=",round(tags*100),"%, list=",round(genes*100),"%, signal=",round(signal*100),"%")
  }
  neg<-subset(df, NES < 0)
  pos<-subset(df, NES > 0)
  if(is.null(prefix)){
    write.table(neg, file="GSEA_neg_report.tsv", quote=F, row.names=F, sep="\t")
    write.table(pos, file="GSEA_pos_report.tsv", quote=F, row.names=F, sep="\t")
    if(!is.null(gmt)){
      if(is.character(gmt)){
        if(grep("http",gmt)>0){
          gmt<-fgsea::gmtPathways(url(gmt))
        } else {
          gmt<-fgsea::gmtPathways(gmt)
        }
      }
      write.gmt(gmt, filename="GSEA_gmt.gmt")
    }
  } else {
    write.table(neg, file=paste0(prefix, "_neg_report.tsv"), quote=F, row.names=F, sep="\t")
    write.table(pos, file=paste0(prefix, "_pos_report.tsv"), quote=F, row.names=F, sep="\t")
    if(!is.null(gmt)){
      if(is.character(gmt)){
        if(grep("http",gmt)>0){
          gmt<-fgsea::gmtPathways(url(gmt))
        } else {
          gmt<-fgsea::gmtPathways(gmt)
        }
      }
      write.gmt(gmt, filename=paste0(prefix, "_gmt.gmt"))
    }
  }
}
