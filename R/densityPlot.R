#' Density Plot
#'
#' Create a plot showing the densities of values in columns of a data frame
#'
#' @param x data.frame object where the density of columnar values will be plotted
#' @param title Character indicating plot title
#' @param xlab Character indicating the x-axis label. Default is 'score'
#' @param ylab Character indicating the y-axis label. Default is 'Frequency'
#' @param ymax Numeric value indicating the maximum y-value. If values are potentially greater than ymax, ymax will be adjusted to relfect the new maximum (rounded to 6 decimal points).
#' @param xmin Numeric value indicating the minimum x-value. Default is NULL, meaning the x-axis limits will range from the negative maximum absolute x-value to the positive maximum absolue x-value.
#' @return A plot with filled densities for each column.
#' @export

dPlot<-function(x, title="", xlab="score", ylab="Frequency", ymax=0, xmin=NULL){
  for(i in 1:ncol(x)){
    ymax<-max(c(ymax, density(x[,i])$y))
  }
  message("setting ymax to ", ymax)
  xmin <- xmin
  if(is.null(xmin)){
    xmin<- -max(abs(x), na.rm=T)
  }
  plot(density(x[,1]), col=1, lty=1, xlab=xlab, ylab=ylab, xlim=c(xmin, max(abs(x), na.rm=T)), ylim=c(0,round(ymax, digits=6)),
       main=title)
  t<-smooth.spline(density(x[,1]))
  xy <- predict(t, seq(min(x[,1]), max(x[,1]), by=0.05)) # Some vertices on the curve
  m <- length(xy$x)
  x.poly <- c(xy$x, xy$x[m], xy$x[1])
  y.poly <- c(xy$y, 0, 0)
  polygon(x.poly, y.poly, col=ggplot2::alpha(1,0.5), border=1)
  for(i in 2:ncol(x)){
    lines(density(x[,i]), col=i, lty=1)
    t<-smooth.spline(density(x[,i]))
    xy <- predict(t, seq(min(x[,i]), max(x[,i]), by=0.05)) # Some vertices on the curve
    m <- length(xy$x)
    x.poly <- c(xy$x, xy$x[m], xy$x[1])
    y.poly <- c(xy$y, 0, 0)
    polygon(x.poly, y.poly, col=ggplot2::alpha(i,0.5), border=i)
  }
  legend("topright", legend=colnames(x), col=c(1:ncol(x)), lty=1)
}
