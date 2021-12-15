#' Load HTSeq-count data from a sampleTable into a countMatrix and colData object
#'
#' Function to use a sampleTable (same format for DESeqDataSetFromHTSeqCount) to generate a count matrix and colData object
#'
#' @param sampleTable Character indicating the path to a tab-delimited (.txt) or comma-separated (.csv) sampleTable, or a data frame.
#' @return A list of length 2. Counts = count matrix and colData = metadata for samples
#' @export

loadData<-function(sampleTable){
  if(is.character(sampleTable)){
    if(length(grep(".txt", sampleTable))>0 | length(grep(".tsv", sampleTable))>0){
      sampleTable<-read.delim(sampleTable)
    }
    if(length(grep(".csv", sampleTable))>0){
      sampleTable<-read.csv(sampleTable)
    }
  } else {
    if(!is.data.frame(sampleTable)){
      stop("sampleTable should be character or data.frame.")
    }
  }
  x<-list()
  for(i in 1:nrow(sampleTable)){
    x[[i]]<-read.delim(sampleTable[i,2], header=F, row.names=1)
    colnames(x[[i]])<-sampleTable[i,1]
    names(x)[i]<-sampleTable[i,1]
  }
  x<-do.call("cbind", x)
  colData<-sampleTable[,c(1,3:ncol(sampleTable))]
  return(list(Counts=x, colData=colData))
}

#' Function to load featureCounts data
#'
#' @param filename Character indicating the path to the featureCounts table to import
#' @return A list of length 2: Counts=featureCount count matrix and annotation= Gene Annotation information from featureCounts
#' @export

loadFC<-function(filename){
  x<-read.delim(filename, skip=1)
  anno<-x[,c(1:6)]
  counts<-x[,c(7:ncol(x))]
  return(list(Counts=counts, annotation=anno))
}
