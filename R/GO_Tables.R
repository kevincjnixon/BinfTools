#' GO Table
#'
#' Make a publication-quality table from your GO results
#'
#' @param GOres GO results table returned from GO_GEM when returnRes=TRUE
#' @param title Character for Table title, leave NULL to avoid a title
#' @param subtitle Character for subtitle of table. Leave NULL to not have subtitle
#' @param sig Boolean indeicating if terms should be ordered by significance. FALSE lets terms be ordered by enrichment. Default=TRUE.
#' @param ts numeric vector of length two indicating minimum and maximum terms size cutoffs for enrichment and significance, respectively. Default=c(10,500)
#' @param colnames Character vector of length 3 indicating the output column names for the Term name, enrichment, and adjusted p-value. Defaults to c("Term", "Enrichment", "FDR")
#' @param numTerm Number indicating the number of terms to print in the table. Default=10.
#' @param retGT Boolean indicating if the gt table object should be returned - useful for custom formatting. Default=FALSE.
#' @param printGT Boolean indicating if the gt table object should be printed to Viewer. Default =TRUE.
#' @return A table of class gt of the top GO results
#' @export

GO_Table<-function(GOres, title=NULL, subtitle=NULL, sig=TRUE, ts=c(10,500), colnames=c("Term", "Enrichment", "FDR"), numTerm=10, retGT=FALSE, printGT=TRUE){
  require(gt, quietly=TRUE)
  require(dplyr, quietly = TRUE)
  foot<-" genes per term"
  if(isTRUE(sig)){
    GOres<-GOres[order(GOres$p_value),]
    GOres<-subset(GOres, term_size < ts[2])
    foot<-paste0("*<",ts[2],foot)
  } else {
    GOres<-GOres[order(GOres$enrichment, decreasing = TRUE),]
    GOres<-subset(GOres, term_size > ts[1])
    foot<-paste0("*>",ts[1],foot)
  }
  GOres<-GOres[1:numTerm,c("term_name","enrichment","p_value")]
  GOres<-GOres[complete.cases(GOres),]
  GOres$enrichment<-round(GOres$enrichment, digits=2)
  GOres$p_value<-formatC(GOres$p_value, format="e", digits=2)

  colnames(GOres)<-colnames

  gt_tbl<- gt(GOres) %>%
    tab_style(style=list(cell_text(weight="bold")),
              locations=cells_column_labels()) %>%
    tab_source_note(source_note=foot) %>%
    tab_style(style=list(cell_text(style="italic")),
              locations=cells_source_notes())

  if(!is.null(title)){
    gt_tbl<-
      gt_tbl %>%
      tab_header(title=title, subtitle=subtitle)
  }

  if(isTRUE(printGT)){
    print(gt_tbl)
  }
  if(isTRUE(retGT)){
    return(gt_tbl)
  }
}

#' GO Table from list of GO results
#'
#' Make a publication-quality table from your GO results
#'
#' @param GOres List of GO results table returned from GO_GEM or clusGO when returnRes=TRUE.
#' @param title Character for Table title, leave NULL to avoid a title
#' @param subtitle Character for subtitle of table. Leave NULL to not have subtitle
#' @param sig Boolean indeicating if terms should be ordered by significance. FALSE lets terms be ordered by enrichment. Default=TRUE.
#' @param ts numeric vector of length two indicating minimum and maximum terms size cutoffs for enrichment and significance, respectively. Default=c(10,500)
#' @param colnames Character vector of length 3 indicating the output column names for the Term name, enrichment, and adjusted p-value. Defaults to c("Term", "Enrichment", "FDR")
#' @param numTerm Number indicating the number of terms to print in the table. Default=10.
#' @param retGT Boolean indicating if the gt table object should be returned - useful for custom formatting. Default=FALSE.
#' @param printGT Boolean indicating if the gt table object should be printed to Viewer. Default =TRUE.
#' @return A table of class gt of the top GO results for each entry in the list of GO results.
#' @export

listGO_Table<-function(GOres, title=NULL, subtitle=NULL, sig=TRUE, ts=c(10,500), colnames=c("Term", "Enrichment", "FDR"), numTerm=10, retGT=FALSE, printGT=TRUE){
  require(gt, quietly=TRUE)
  require(dplyr, quietly = TRUE)
  if(length(which(unlist(lapply(GOres, class))  != "data.frame"))>0){
    clusters<-names(which(unlist(lapply(GOres, class))  != "data.frame"))
    message("No results for ",toString(clusters),"... Removing them." )
  }
  GOres<-GOres[which(unlist(lapply(GOres, class)) == "data.frame")]
  foot<-" genes per term"
  if(isTRUE(sig)){
    GOres<-lapply(GOres, function(x){
      x<-x[order(x$p_value),]
      x<-subset(x, term_size < ts[2])})
    foot<-paste0("*<",ts[2],foot)
  } else {
    GOres<-lapply(GOres, function(x){
      x<-x[order(x$enrichment, decreasing=TRUE),]
      x<-subset(x, term_size > ts[1])})
    foot<-paste0("*>",ts[1],foot)
  }
  GOres<-lapply(GOres, function(x) x[1:numTerm,c("term_name","enrichment","p_value")])
  GOres<-lapply(GOres, function(x) x[complete.cases(x),])
  GOres<-lapply(GOres, function(x){
    x$enrichment<-round(x$enrichment, digits=2)
    x$p_value<-formatC(x$p_value, format="e", digits=2)
    return(x)})

  GOres<-lapply(GOres, function(x) {colnames(x)<-colnames
  return(x)})
  clus<-rep(names(GOres), unlist(lapply(GOres, nrow)))
  GOres<-do.call("rbind", GOres)
  GOres$cluster<-clus

  gt_tbl<- gt(GOres, groupname_col="cluster") %>%
    tab_style(style=list(cell_text(style="italic", weight="bold")),
              locations=cells_row_groups())%>%
    tab_style(style=list(cell_text(weight="bold", size="large")),
              locations=cells_column_labels()) %>%
    tab_source_note(source_note=foot) %>%
    tab_style(style=list(cell_text(style="italic")),
              locations=cells_source_notes())

  if(!is.null(title)){
    gt_tbl<-
      gt_tbl %>%
      tab_header(title=title, subtitle=subtitle)
  }

  if(isTRUE(printGT)){
    print(gt_tbl)
  }
  if(isTRUE(retGT)){
    return(gt_tbl)
  }
}

#' Pander table of top n GO results
#'
#' Make a Pandoc table of the top n results from a GO_GEM() results object. This is useful for making markdown reports.
#'
#' @param GOres A single results table from GO_GEM() function with returnRes = T.
#' @param ts Numeric vector of length 2 indicating the minimum and maximum term sizes to filter GO results. Defaults to c(10,500).
#' @param sig Boolean indicating if the top significant terms should be shown. If FALSE, top enriched terms are shown. Default is TRUE.
#' @param title Character for the table legend description
#' @param n number of top terms to show. Default is 10. If there are fewer than n terms after filtering, only those terms will be shown.
#' @return Prints a pandoc.table of results to stdout.
#' @export

GOtab<-function(GOres, ts=c(10,500), sig=TRUE, title="", n=10){
  GOres<-subset(GOres, term_size >= ts[1] & term_size <= ts[2])
  if(isTRUE(sig)){
    GOres<-GOres[order(GOres$p_value, decreasing=F),]
  } else {
    GOres<-GOres[order(GOres$enrichment, decreasing=T),]
  }
  tmp<-data.frame(Term=GOres$term_name[1:n], FDR=GOres$p_value[1:n], Enrichment=GOres$enrichment[1:n])
  tmp<-tmp[complete.cases(tmp),]
  if(nrow(tmp)>0){
    pander::pandoc.table(tmp, emphasize.rownames=FALSE, caption=title)
  }
}
