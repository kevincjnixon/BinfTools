#' Plot top n pathways according to given stat
#'
#' @param fgsea GSEA results from GSEA() function
#' @param n number of pathways to plot, default is 30.
#' @param by The stat for selecting top pathways, default is padj
#' @param decreasing Should the sort order be increasing or decreasing? Default is False.
#' @param parseBader Boolean indicating if gene set names should be parsed for the output figure following Bader Lab nomenclature (i.e. '%' as delimiter)
#' @export
barplotTopGSEA <- function(fgRes, n = 30, by = "padj", decreasing = FALSE, parseBader = T){

  if (!by %in% colnames(fgRes)){
    stop("The input dataframe does not contain the provided stat!")
  }
  fgRes <- fgRes[order(fgRes[[by]], decreasing = decreasing), ]
  filtRes = head(fgRes, n = n)

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

}
