#'
#' \code{tocdf} This function ...
#'
#' @export

tocdf <- function(inset, scol, ecol){
  mattemp =apply(t(inset[,scol:ecol]), 2, cumsum)
  outset = inset
  outset[,scol:ecol]=t(mattemp)
  return(outset)
}