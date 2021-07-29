#' Make a Venn diagram
#'
#' Make a Venn diagram from a named list of genes. Uses functions from the VennDiagram package.
#'
#' @param x A named list (max length of 5) of characters (genes) to plot in a Venn Diagram
#' @param cols Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used. Default is "Dark2"
#' @param lty Line type (1=solid line, 2=dashed line, default="blank")
#' @return An image of a Venn diagram showing overlaps between groups
#' @export
