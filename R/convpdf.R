#'
#' \code{convpdf} This function ...
#'
#' @export

convpdf<-function (origpdf, origx, convx){
  # input
  #origpdf
  #origx
  #convx
  #
  olen <- length(origx)
  clen <- length(convx)
  convpdf <- vector('numeric', length=clen)
  convpdf[1]<- 0
  for(cpos in 2:(clen-1)){
    opos <- max(sum(ifelse(origx < convx[cpos], 1, 0)),1) 
    cslope<- (origpdf[opos+1]-origpdf[opos])/(origx[opos+1]-origx[opos])
    convpdf[cpos] <-origpdf[opos]+(convx[cpos]-origx[opos])*cslope
  }
  convpdf[clen]=1
  return(convpdf)
}