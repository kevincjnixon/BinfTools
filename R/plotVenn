getOL<-function(x){
    if(length(x)>5){
      stop("max length of x is 5...")
    }
    OL<-list(n12=NA, n13=NA, n14=NA, n15=NA,
             n23=NA, n24=NA, n25=NA,
             n34=NA, n35=NA, n45=NA,
             n123=NA, n124=NA, n125=NA, n134=NA, n135=NA, n145=NA)
    OL=list()
    index<-1
    #start with 1 on one comparisons:
    for(i in 1:length(x)){
      if(i < length(x)){
        for(j in (i+1):length(x)){
          OL[[index]]<-length(x[[i]][which(x[[i]] %in% x[[j]])])
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
                OL[[index]]<-length(tmp[which(tmp %in% x[[k]])])
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
                    OL[[index]]<-length(tmp[which(tmp %in% x[[d]])])
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
        tmp<-tmp[which(tmp) %in% x[[i]]]
      }
      OL[[index]]<-length(tmp)
      names(OL)[index]<-"n12345"
    }
    print(names(OL))
    return(OL)
 }
  
#' Make a Venn diagram
#'
#' Make a Venn diagram from a named list of genes. Uses functions from the VennDiagram package.
#'
#' @param x A named list (max length of 5) of characters (genes) to plot in a Venn Diagram
#' @param cols Character indicating the RColorBrewer palette name or list of colours (hex, name, rgb()) to be used. Default is "Dark2"
#' @param lty Line type (1=solid line, 2=dashed line, default="blank")
#' @return An image of a Venn diagram showing overlaps between groups
#' @export

plotVenn<-function(x, cols="Dark2", lty="blank"){
  require(VennDiagram)
  y<-getOL(x)
  grid.newpage()
  if(length(x)==1){
    draw.single.venn(area=length(x),
                     fill=BinfTools::colPal(cols)[1:length(x)],
                     lty=lty,
                     category=names(x))
  }
  if(length(x)==2){
    draw.pairwise.venn(area1=length(x[[1]]),
                       area2=length(x[[2]]),
                       cross.area=y$n12,
                       fill=BinfTools::colPal(cols)[1:length(x)],
                       lty=lty,
                       category = names(x))
  }
  if(length(x)==3){
    #print(lengths(y))
    draw.triple.venn(area1=length(x[[1]]),
                     area2=length(x[[2]]),
                     area3=length(x[[3]]),
                     n12=y$n12,
                     n23=y$n23,
                     n13=y$n13,
                     n123=y$n123,
                     fill=BinfTools::colPal(cols)[1:length(x)],
                     lty=lty,
                     category = names(x))
  }
  if(length(x)==4){
    draw.quad.venn(area1=length(x[[1]]),
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
    draw.quintuple.venn(area1=length(x[[1]]), area2=length(x[[2]]),
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
}
