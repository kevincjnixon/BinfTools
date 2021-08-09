getOL<-function(x, retVals=F){
    if(length(x)>5){
      stop("max length of x is 5...")
    }
    OL=list()
    index<-1
    #start with 1 on one comparisons:
    for(i in 1:length(x)){
      if(i < length(x)){
        for(j in (i+1):length(x)){
          OL[[index]]<-x[[i]][which(x[[i]] %in% x[[j]])]
          names(OL)[index]<-paste0("n",i,j)
          index<-index+1
        }
      }
    }
    if(length(x)>2){
      #Now three way comparisons
      for(i in 1:length(x)){
        if(i<length(x)){
          for(j in (i+1):length(x)){
            if(j<length(x)){
              tmp<-x[[i]][which(x[[i]] %in% x[[j]])]
              for(k in (j+1):length(x)){
                OL[[index]]<-tmp[which(tmp %in% x[[k]])]
                names(OL)[index]<-paste0("n",i,j,k)
                index<-index+1
              }
            }
          }
        }
      }
    }
    if(length(x)>3){
      #Now for the four way comparisons
      for(i in 1:length(x)){
        if(i<length(x)){
          for(j in (i+1):length(x)){
            if(j<length(x)){
              tmp<-x[[i]][which(x[[i]] %in% x[[j]])]
              for(k in (j+1):length(x)){
                if(k<length(x)){
                  tmp<-tmp[which(tmp %in% x[[k]])]
                  for(d in (k+1):length(x)){
                    OL[[index]]<-tmp[which(tmp %in% x[[d]])]
                    names(OL)[index]<-paste0("n",i,j,k,d)
                    index<-index+1
                  }
                }
              }
            }
          }
        }
      }
    }
    if(length(x)==5){
      #Now for the final comparison:
      tmp<-x[[1]]
      for(i in 2:length(x)){
        tmp<-tmp[which(tmp %in% x[[i]])]
      }
      OL[[index]]<-tmp
      names(OL)[index]<-"n12345"
    }
    #print(names(OL))
    if(isFALSE(retVals)){
      OL<-as.list(lengths(OL))
    }
    return(OL)
 }
  
nonOL<-function(x, OL){
  OLgenes<-unique(unlist(OL))
  y<-list()
  for(i in 1:length(x)){
    y[[i]]<-x[[i]][which(!x[[i]] %in% OLgenes)]
    names(y)[i]<-paste("Only",i, sep="_")
  }
  return(y)
}

#' Make a Venn diagram
#'
#' Make a Venn diagram from a named list of genes. Uses functions from the VennDiagram package.
#'
#' @param x A named list (max length of 5) of characters (genes) to plot in a Venn Diagram
#' @param title Character indicating the title of the plot. Default="Venn Diagram"
#' @param cols Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used. Default is "Dark2"
#' @param lty Line type (1=solid line, 2=dashed line, default="blank")
#' @param scale Boolean indicating if circles should be scaled to size (works only for 2 or 3-way Venn diagrams). Default=F.
#' @param retVals Boolean indicating if the overlaps should be returned (default=F)
#' @return An image of a Venn diagram showing overlaps between groups
#' @export

plotVenn<-function(x, title="Venn Diagram", cols="Dark2", lty="blank", scale=F, retVals=F){
  require(VennDiagram, quietly=T)
  y<-getOL(x)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(2,1, heights=unit(c(0.25, 10),"null"))))
  venn<-NULL
  if(length(x)==1){
    venn<-draw.single.venn(area=length(x),
                     fill=BinfTools::colPal(cols)[1:length(x)],
                     lty=lty,
                     category=names(x))
  }
  if(length(x)==2){
    venn<-draw.pairwise.venn(area1=length(x[[1]]),
                       area2=length(x[[2]]),
                       cross.area=y$n12,
                       fill=BinfTools::colPal(cols)[1:length(x)],
                       lty=lty,
                       category = names(x),
                       scaled=scale,
                       euler.d=scale)
  }
  if(length(x)==3){
    if(isTRUE(scale)){
        assign("overrideTriple", TRUE, envir=.GlovalEnv)
    }
    venn<-draw.triple.venn(area1=length(x[[1]]),
                     area2=length(x[[2]]),
                     area3=length(x[[3]]),
                     n12=y$n12,
                     n23=y$n23,
                     n13=y$n13,
                     n123=y$n123,
                     fill=BinfTools::colPal(cols)[1:length(x)],
                     lty=lty,
                     category = names(x),
                     scaled=scale,
                     euler.d=scale)
    if(isTRUE(scale)){
       rm(overrideTriple, pos=".GlobalEnv")
    }
  }
  if(length(x)==4){
    venn<-draw.quad.venn(area1=length(x[[1]]),
                   area2=length(x[[2]]),
                   area3=length(x[[3]]),
                   area4=length(x[[4]]),
                   n12=y$n12, n13=y$n13, n14=y$n14,
                   n23=y$n23, n24=y$n24, n34=y$n34,
                   n123=y$n123, n124=y$n124, n134=y$n134, n234=y$n234,
                   n1234=y$n1234,
                   fill=BinfTools::colPal(cols)[1:length(x)],
                   lty=lty,
                   category=names(x))
  }
  if(length(x)==5){
    venn<-draw.quintuple.venn(area1=length(x[[1]]), area2=length(x[[2]]),
                        area3=length(x[[3]]), area4=length(x[[4]]), area5=length(x[[5]]),
                        n12=y$n12, n13=y$n13, n14=y$n14, n15=y$n15, n23=y$n23, n24=y$n24,
                        n25=y$n25, n34=y$n34, n35=y$n35, n45=y$n45, n123=y$n123, n124=y$n124,
                        n125=y$n125, n134=y$n134, n135=y$n135, n145=y$n145, n234=y$n234, n235=y$n235,
                        n245=y$n245, n345=y$n345, n1234=y$n1234, n1235=y$n1235, n1245=y$n1245,
                        n1345=y$n1345, n2345=y$n2345, n12345=y$n12345,
                        fill=BinfTools::colPal(cols)[1:length(x)],
                        lty=lty,
                        category=names(x))
  }
  #print(venn, vp=viewport(layout.pos.row=2))
  grid.text(title, vp=viewport(layout.pos.row=1))
  
  if(isTRUE(retVals)){
      res<-getOL(x, retVals=retVals)
      res<-c(nonOL(x, res), res)
      return(res)
  }
}
