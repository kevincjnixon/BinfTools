##### Functions for upset plot-like overlaps #####
fromList2<-function (input) 
{
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  #Now add original value names as rownames
  row.names(data)<-elements
  return(data)
}

overlapGroups <- function (listInput, sort = TRUE) {
  listInputmat    <- fromList2(listInput) == 1
  # condensing matrix to unique combinations elements
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  # going through all unique combinations and collect elements for each in a list
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- names(myelements)
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  return(grouplist)
}

filtOL<-function(OLs, delim, length=2){
  #Start by splitting the names of the list using the delimiter
  nameList<-strsplit(names(OLs),delim, T)
  #Now filter down to only entries matching the desired length
  OLs<-OLs[which(lengths(nameList) == length)]
  nameList<-nameList[which(lengths(nameList)==length)]
  #Now we want to filter to the entries in which the names from all parties match each other
  #Best way to do this would be to check the length of unique names and pick those with just '1'
  OLs<-OLs[which(unlist(lapply(nameList, function(x){return(length(unique(x)))})) == 1)]
  nameList<-nameList[which(unlist(lapply(nameList, function(x){return(length(unique(x)))})) == 1)]
  return(OLs)
}