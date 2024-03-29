#'Convert to gene symbols
#'
#'This function makes use of gprofiler2's 'gconvert' function to convert gene
#' ids in the rownames of the results object or the normalized counts object to
#' a gene symbol for easy labeling in plots.
#'
#' This will be integrated into the volcanoPlot(), MA_Plot(), and zheat() functions
#' in the future.
#'
#' @param object a data.frame object (results or normalized counts) with the rownames as the gene IDs to convert or a character vector of gene IDs to convert
#' @param obType Character vector (one of 'res', 'counts', or 'vector') indicating whether *object* is a results object, normalized counts object, or character vector of gene IDs, respectively.
#' @param species Character vector indicating the species (e.g. "hsapiens" (*default*), "mmusculus").
#' @param target Character vector indicating which format of the gene IDs to return (e.g. "HGNC", ENSG", "REFSEQ_MRNA","ENTREZGENE", see https://biit.cs.ut.ee/gprofiler/convert for all options).
#' @param addCol Boolean indicating if an additional column named "SYMBOL" should be added, leaving the rownames as the original format. If FALSE (*default*), rownames of *object* will be replaced with gene symbols - duplicate gene symbols will be handled by keeping the symbol with the highest gene expression values.
#' @return A data.frame object in the format indicated by *obType* with gene symbols as the rownames (if *addCol*=FALSE) or gene symbols added as their own column named "SYMBOL".
#' @export

getSym<-function(object, obType=c("res","counts","vector"), species="hsapiens", target="HGNC", addCol=FALSE){
  x<-as.data.frame(object)
  if(obType == "vector"){
    genes<-object
  } else {
    genes<-rownames(object)
  }
  if(target == "ENSGV"){
    genes<-sapply(strsplit(genes, split=".", fixed=TRUE), '[[',1)
    target<-"ENSG"
  }
  if(obType == "vector"){
    y<-gprofiler2::gconvert(genes, organism=species, target=target, mthreshold=1, filter_na=FALSE)
    sym<-y$target
    return(sym)
  }
  if(obType == "res"){
    y<-gprofiler2::gconvert(genes, organism=species, target=target, mthreshold=1, filter_na=FALSE)
    sym<-y$target
    if(isTRUE(addCol)){
      x$SYMBOL<-sym
    }else{
      x$tmp<-sym
      #Order based on expression (decreasing)
      x<-x[order(x$baseMean, decreasing=TRUE),]
      #Remove NAs
      x<-x[complete.cases(x$tmp),]
      #Check to see if there are any duplicated symbols in x$tmp
      dups<-unique(x$tmp[duplicated(x$tmp)])
      if(length(dups)>=1){
        print("Duplicate gene symbols detected. Keeping symbols with highest overall expression...")
        to_remove<-c()
        for(i in 1:length(dups)){
          #Get the array row numbers of the duplicates (first will be highest expression)
          ind<-which(x$tmp%in%dups[i], arr.ind=TRUE)
          to_remove<-c(to_remove, ind[-1])
        }
        #Remove the rows
        x<-x[-to_remove,]
      }
      #Set rownames to symbol and remove tmp column
      rownames(x)<-x$tmp
      x<-x[,-(which(colnames(x) %in% "tmp"))]
    }
  }
  if(obType == "counts"){
    y<-gprofiler2::gconvert(genes, organism=species, target=target, mthreshold=1, filter_na=FALSE)
    sym<-y$target
    if(isTRUE(addCol)){
      x$SYMBOL<-sym
    }else{
      x$tmp<-sym
      #Order based on expression (decreasing)
      x<-x[order(rowMeans(object), decreasing = TRUE),]
      #Remove NAs
      x<-x[complete.cases(x$tmp),]
      #Check to see if there are any duplicated symbols in x$tmp
      dups<-unique(x$tmp[duplicated(x$tmp)])
      if(length(dups)>=1){
        print("Duplicate gene symbols detected. Keeping symbols with highest overall expression...")
        to_remove<-c()
        for(i in 1:length(dups)){
          #Get the array row numbers of the duplicates (first will be highest expression)
          ind<-which(x$tmp%in%dups[i], arr.ind=TRUE)
          to_remove<-c(to_remove, ind[-1])
        }
        #Remove the rows
        x<-x[-to_remove,]
      }
      #Set rownames to symbol and remove tmp column
      rownames(x)<-x$tmp
      x<-x[,-(which(colnames(x) %in% "tmp"))]
    }
  }
  return(x)
}
