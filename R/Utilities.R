#' grep for multiple patterns
#'
#' Search for multiple regex patterns in a given vector
#'
#' @param p Character vector of regex patterns to search
#' @param x Target vector to search for patterns
#' @param value Boolean indicating if the value of matching patterns should be returned. If FALSE, only the index of matching patterns will be returned. Default is FALSE.
#' @param ic Boolean indicating if case should be ignored. Default is FALSE.
#' @param v Boolean indicating if the returned values should be inverted. Default is FALSE.
#' @param unique Boolean indicating if duplicate values should be removed before returning. Default is FALSE.
#' @return Vector of indices of x that contain matching patterns of p. If value = TRUE, values of x containing p will be returned.
#' @export

mgrep<-function(p,x, value=F, ic=F, v=F, unique=F){
  res<-c()
  for(i in p){
    res<-c(res, grep(i,x,ic,value=value, invert=v))
  }
  if(isTRUE(unique)){
    res<-unique(res)
  }
  return(res)
}
