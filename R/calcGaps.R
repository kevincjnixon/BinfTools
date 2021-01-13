#'Calculate optimal number of clusters
#'
#'Calculate the gap statistic for a series of k-means clusterings using the SAGx
#'package. Use the number of clusters with the highest gap statistic for your
#'clustering analysis.
#'
#'@param resList List of results data frames
#'@param maxClus Maximum number of clusters to test, minimum 3. Default 10. NOTE: The more, the longer it takes.
#'@param iter Number of bootstrapping iterations for calculating gap statistic. Default 10. NOTE: the more, the longer it takes.
#'@return A vector of length 2 containing the maximum gap statistic and the optimal number of k-means clusters, respectively.
#'@export

calcGaps<-function(resList, maxClus=10, iter=10){
  #Make sure there is more than one entry in the list
  if(length(resList)<2){
    stop("resList must have at least 2 entries...")
  }
  #Make sure that maxClus is greater than 2
  if(maxClus<3){
    stop("maxClus must be at least 3...")
  }
  #Get the names of the comparisons
  compNames<-names(resList)
  #Create a matrix for clustering based on the log2FoldChanges of each comparison
  mat<-data.frame(genes=rownames(resList[[1]]), comp1=resList[[1]]$log2FoldChange)
  for(i in 2:length(resList)){
    tmp<-data.frame(genes=rownames(resList[[i]]), tmp=resList[[i]]$log2FoldChange)
    mat<-merge(mat, tmp, by.x="genes", by.y="genes")
  }
  #Set genes to rownames
  rownames(mat)<-mat$genes
  #remove 'genes' column
  mat<-mat[,-(which(colnames(mat) %in% "genes"))]
  colnames(mat)<-compNames
  #Remove any NAs
  mat<-as.matrix(mat[complete.cases(mat),])
  #Calculate optimal number of clusters:
  message(paste("Calculating gap statistics for",maxClus,"k-means clusters..."))
  gaps<-list()
  pb<-txtProgressBar(min=2, max=maxClus, initial=2, style=3)
  for(i in 2:maxClus){
    #Set the seed
    set.seed(1234)
    tryCatch({ kset<-kmeans(mat, centers=i)
    gaps[[i]]<-SAGx::gap(mat, kset$cluster, B=iter)} ,
    warning=function(e){
      message(paste0("Warning occurred clustering with ",i," centres..."))
      #Put 0 as gaps in place:
      gaps[[i]]<-c(0,0)
    }
    )
    gc()
    setTxtProgressBar(pb,i)
  }
  #Put all the gap statistics in a single vector
  x<-unlist(sapply(gaps, '[[',1))
  #Max Gap
  res<-c(x[which(x==max(x))], which(x==max(x), arr.ind=T)+1)
  names(res)<-c("Max Gap Statistic","Optimal Cluster Number")
  return(res)
}
