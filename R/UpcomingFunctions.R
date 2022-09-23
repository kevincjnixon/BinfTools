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
plotCor<-function(x,y, xCol=F, yCol=F, xlab="", ylab="", title="", scale=F){
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
  #print(length(x))
  #print(length(y))
  if(isTRUE(scale)){
    messgae("z-score...")
    x<-scale(x)
    y<-scale(y)
  }
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

#3D volcano plot - for combining two different contrasts - pvalues are combined using Fisher's combined probability test
vol3d<-function(x, y, xlab="res1", ylab="res2", pval=0.05, pcol="pvalue", usep=NULL){
  x<-x[order(rownames(x)),]
  y<-y[order(rownames(y)),]
  x<-x[which(rownames(x) %in% rownames(y)),]
  y<-y[which(rownames(y) %in% rownames(x)),]
  zlab<-pcol
  if(!all.equal(rownames(x), rownames(y))){
    stop("Check again!")
  }
  pcol<-which(colnames(x) %in% pcol)
  pvals<- -log10(x[,pcol])
  if(is.null(usep)){
    pvals<-c()
    for(i in 1:nrow(x)){
      pvals<-c(pvals, survcomp::combine.test(c(x[,pcol][i], y[,pcol][i])))
    }
    pvals<- -log10(pvals)
  } else {
    if(usep[1]==2){
      pvals<- -log10(y$pvalue)
    }
  }
  for3d<-data.frame(x=x$log2FoldChange, y=y$log2FoldChange, z=pvals, row.names=rownames(x))
  for3d$col<-ifelse(for3d$x > 0 & for3d$y >0 & for3d$z > -log10(pval), rgb(1,0,0,0.5), rgb(0,0,0,0.5))
  for3d$DE<-ifelse(for3d$x > 0 & for3d$y >0 & for3d$z > -log10(pval), "Up Both", "No Change")
  for3d$col[which(for3d$x <0 & for3d$y <0 & for3d$z > -log10(pval))] <- rgb(0,0,1,0.5)
  for3d$DE[which(for3d$x <0 & for3d$y <0 & for3d$z > -log10(pval))] <- "Down Both"
  for3d$col[which(for3d$x <0 & for3d$y >0 & for3d$z > -log10(pval))] <- rgb(0,1,0,0.5)
  for3d$DE[which(for3d$x <0 & for3d$y >0 & for3d$z > -log10(pval))] <- paste0("Down ", xlab,", Up ", ylab)
  for3d$col[which(for3d$x >0 & for3d$y <0 & for3d$z > -log10(pval))] <- rgb(0.5,0,1,0.5)
  for3d$DE[which(for3d$x >0 & for3d$y <0 & for3d$z > -log10(pval))] <- paste0("Up ",xlab,", Down ",ylab)
  scatterplot3d::scatterplot3d(for3d[,c(1:3)], pch=16, xlab=paste("log2FoldChange -", xlab), ylab=paste("log2FoldChange -", ylab), zlab=paste0("-log10(",zlab,")"), color=for3d$col, 
                               grid=T, box=F)
  legend("topright", legend=c("Up Both", "Down Both", paste0("Up ",xlab,", Down ", ylab), paste0("Down ",xlab,", Up ",ylab), "No Change"),
         col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5), rgb(0.5,0,1,0.5), rgb(0,1,0,0.5), rgb(0,0,0,0.5)), pch=16)
  return(for3d)
}

#Make a heatmap of gene expression grouped by gene sets
gmtHeat<-function(counts, cond, gmt, con=NULL, labgenes=NULL, avgExp=T, zscore=T, retGroups=F){
  #Filter gmt to have only genes found in rownames(counts)
  gmt<-lapply(gmt, function(x){return(x[which(x %in% rownames(counts))])})
  #Make an annotation data frame:
  annodf<-suppressWarnings(unique(tidyr::gather(as.data.frame(do.call("cbind", gmt)), key="Term", value="Genes")))
  rownames(annodf)<-make.names(annodf$Genes, unique=T)
  if(isTRUE(zscore)){
    message("Z-scoring counts")
    counts<-t(scale(t(counts)))
  }
  #Custom Order columns
  if(!is.null(con)){
    if (length(con) == length(levels(factor(cond)))) {
      cond <- forcats::fct_relevel(cond, con)
    }
    else {
      ind <- which(levels(cond) == con)
      if (ind != 1) {
        x <- 1:length(levels(cond))
        x <- x[-ind]
        x <- c(ind, x)
        y <- levels(cond)[x]
        cond <- forcats::fct_relevel(cond, y)
      }
    }
    tmp.counts <- counts[, which(cond == levels(cond)[1])]
    tmp.cond <- as.character(cond[which(cond == 
                                                      levels(cond)[1])])
    for (i in 2:length(levels(cond))) {
      tmp.counts <- cbind(tmp.counts, counts[, which(cond == 
                                                 levels(cond)[i])])
      tmp.cond <- c(tmp.cond, as.character(cond[which(cond == 
                                                                          levels(cond)[i])]))
    }
    counts <- tmp.counts
    cond <- tmp.cond
  }
  if(isTRUE(avgExp)){
    message("Averaging expression values accross replicates")
    counts<-avgExp(counts, cond)
    cond<-colnames(avgExp)
  }
  #Make the heatmap data frames
  forHeat<-list()
  for(i in 1:length(gmt)){
    forHeat[[i]]<-counts[which(rownames(counts) %in% gmt[[i]]),]
    forHeat[[i]]<-forHeat[[i]][match(gmt[[i]], rownames(forHeat[[i]])),]
    names(forHeat)[i]<-names(gmt)[i]
  }
  if(isTRUE(retGroups)){
    message("Returning list of groups")
    return(forHeat)
  }
  forHeat<-suppressWarnings(as.data.frame(do.call("rbind", forHeat)))
  rownames(forHeat)<-rownames(annodf) 
  gaps <- c()
  if (length(gmt) < 2) {
    gaps <- NULL
  } else {
    for (i in 1:(length(gmt) - 1)) {
      gaps <- c(gaps, sum(gaps[length(gaps)], length(gmt[[i]])))
    }
  }
  if (!is.null(labgenes)) {
    tmp <- rep(" ", nrow(forHeat))
    for (i in 1:length(labgenes)) {
      tmp[which(annodf$Genes %in% labgenes[i], arr.ind = T)] <- labgenes[i]
    }
    labgenes <- tmp
  } else {
    labgenes<-annodf$Genes
  }
  pheatmap::pheatmap(forHeat, annotation_row=annodf[,-2,drop=F], gaps_row=gaps, cluster_rows = F,
                     labels_row=labgenes, cluster_cols=F)
}

#Filter a list of GO results to either remove tables with no significant results (replace=F), or replace them with a dummy table (replace=T). 
#Replacing with a dummy table is good when you want to show all GO analysis names in a GOHeat analysis.
filtGO<-function(x, replace=F){
  to_rm<-c()
  for(i in 1:length(x)){
    if(is.null(nrow(x[[i]]))){
      to_rm<-c(to_rm,i)
    }
  }
  if(length(to_rm)>0){
    if(isFALSE(replace)){
      message("Removing ",length(to_rm)," GO tables with no results.")
      x<-x[-to_rm]
    } else {
      message("Replacing ", length(to_rm)," GO tables with no results with a dummy table.")
      tmp<-data.frame(p_value=c(1,1), term_name=rep("No Enriched Terms",2))
      for(i in 1:length(to_rm)){
        x[[to_rm[i]]]<-tmp
      }
    }
  }
  return(x)
}
