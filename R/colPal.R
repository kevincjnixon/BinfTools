#' Show colour palettes from RColorBrewer
#' @return Image of colour palettes from RColorBrewer with names
#' @export

showPals<-function(){
  par(mar=c(3,4,2,2))
  RColorBrewer::display.brewer.all()
}

#'Colour palettes
#'
#'This function creates either a custom colour palette or can return the colours
#'from an RColorBrewer palette
#'
#'@param col A character vector indicating the custom colours ("rgb", "colour name", "Hexadecimal") or RColourBrewer palette to be used. If list of colour values, the must all be of the same format ("color name", "hexadecimal", or "rgb()") If RColourBrewer palette name, length must be 1.
#'@return A character vector of hexadecimal values representing the colour palette
#' @export

colPal<-function(col){
  if(length(col)>1){
    if(length(grep("#", col))>0){
      #If rgb or hex is provided, they're already in hex format, so we can return that
      return(col)
    } else {
      #must be colour names (e.g."red", "blue"...)
      #Need to convert to hex
      return(gplots::col2hex(col))
    }
  } else {
    tryCatch({
      col<-suppressWarnings(RColorBrewer::brewer.pal(Inf, col))
      return(col)},
      error=function(e){
        #message("Not RColorBrewer palette, trying colour name...")
        if(length(grep("#", col))>0){
          return(col)
        } else {
          tryCatch(return(gplots::col2hex(col)),
                   error=function(e){
                     message("Not a valid colour.")
                     stop()
                   })
        }
      }
      )
  }
}
