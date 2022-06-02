#Take a numeric vector (i.e. a vector of gene expression) and assign it into named quantiles of a desired length - great for making a 'condition' vector based off of gene expression
quart_group<-function(x, groups=c("1st","2nd","3rd","4th")){
  q<-quantile(x, probs=seq(0,1, (1/length(groups))))
  g<-c(rep(groups[length(groups)], length(x)))
  for(i in 1:length(g)){
    for(k in (length(q)-1):2){
      if(x[i]<q[k]){
        g[i]<-groups[k-1]
      }
    }
  }
  return(g)
}

#Grep a certain key from a GO results table's "term_name" column, and return provided columns
#Great for making a termList for GOHeat (return "term_name" only) or getting intersection genes
GOgrep<-function(GOtable, key, cols=c("term_name","intersection")){
  x<-GOtable[grep(key, GOtable$term_name, ignore.case=T),which(colnames(GOtable) %in% cols), drop=F]
  if(nrow(x)<1){
    x<-NA
  }
  return(x)
}

na.rm<-function(x){
  return(x[!is.na(x)])
}
