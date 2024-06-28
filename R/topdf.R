#'
#' \code{topdf} This function ...
#'
#' @export

topdf <- function(inset,scol,ecol){
  mattemp = apply(t(inset[,scol:ecol]), 2, function(x){x/sum(x)})
  sumtemp = apply(inset[,scol:ecol], 1, sum)
  outset = inset
  outset[,scol:ecol]=t(mattemp)#normalized function
  outset[,ecol+1] = sumtemp#sum of original function
  #      browser()
  return(outset)
}