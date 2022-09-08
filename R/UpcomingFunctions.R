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

#Plot a correlation between two columns from two different or the same data frame or two named vectors
#Great for comparing gene expression changes or values, or two gene rankings for GSEA
plotCor<-function(x,y, xCol=F, yCol=F, xlab="", ylab="", title=""){
  if(isFALSE(xCol)){
    if(isFALSE(yCol)){
      x<-x[which(names(x) %in% names(y))]
      y<-y[which(names(y) %in% names(x))]
    } else {
      x<-x[which(names(x) %in% rownames(y))]
      tmpNames<-rownames(y)[which(rownames(y) %in% names(x))]
      y<-y[which(rownames(y) %in% names(x)),which(colnames(y) %in% yCol)]
      names(y)<-tmpNames
    }
  } else {
    if(isFALSE(yCol)){
      tmpNames<-rownames(x)[which(rownames(x) %in% names(y))]
      x<-x[which(rownames(x) %in% names(y)),which(colnames(x) %in% xCol)]
      names(x)<-tmpNames
      y<-y[which(names(y) %in% names(x))]
    } else {
      tmpNames<-rownames(x)[which(rownames(x) %in% rownames(y))]
      x<-x[which(rownames(x) %in% rownames(y)), which(colnames(x) %in% xCol)]
      names(x)<-tmpNames
      tmpNames<-rownames(y)[which(rownames(y) %in% names(x))]
      y<-y[which(rownames(y) %in% names(x)), which(colnames(y) %in% yCol)]
      names(y)<-tmpNames
    }
  }
  x<-x[order(names(x))]
  y<-y[order(names(y))]
  print(length(x))
  print(length(y))
  fit<-lm(y~x)
  plot(x,y, main=title, pch=16, xlab=xlab, ylab=ylab)
  abline(fit, lty=2, col="red")
  legend("topleft", legend=paste0("R=",round(sqrt(summary(fit)$adj.r.squared),digits=4)))
}

#If you have n=1 for one sample, you can still make a 'pseudo MA plot'
pseudoMA<-function(counts, treat, con, FoldChange=1, title="", retDA=F, retRes=F){
  x<-data.frame(row.names=rownames(counts),
                baseMean=rowMeans(counts[,which(colnames(counts) %in% c(treat, con)), drop=F]),
                log2FoldChange=log2((1+rowMeans(counts[,which(colnames(counts) %in% treat),drop=F]))/(1+rowMeans(counts[,which(colnames(counts) %in% con),drop=F]))),
                pvalue=rep(0.5,nrow(counts)),
                padj=rep(0.5, nrow(counts)))
  #Filter genes with baseMean==0
  x<-subset(x, baseMean>0)
  print(range(x$baseMean))
  DA<-BinfTools::MA_Plot(x, title, 1, NULL, FC=FoldChange, returnDEG=T)
  if(isTRUE(retDA)){
    return(DA)
  }
  if(isTRUE(retRes)){
    return(x)
  }
}
