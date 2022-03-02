#'Make a heatmap showing significance of groups of GO terms from multiple sets of results
#'
#'@param GOresList list of GO_GEM results (set returnRes=T) of length >= 2.
#'@param termList named list of character vectors where the names represent the group names for GO terms and vectors contain GO term names corresponding to the 'term_name' column in the GO_GEM results data frams.
#'@param hmcol colorRampPalette of length 100 that will direct the colour palette of the heatmap. Default is colorRampPalette(c("white","darkblue"))(100)
#'@param width numeric indicating the cell width (Default=NA will automatically direct the cell width)
#'@param height numeric indicating the cell height (Default=NA will automatically direct the cell height)
#'@return A grouped heatmap showing significance of GO terms across analyses as -log10(p-value)
#'@export
