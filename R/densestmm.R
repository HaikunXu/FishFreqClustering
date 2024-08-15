#'
#' \code{densestmm} This function ...
#'
#' @export

densestmm <- function(densx, samplvec, bw=0){ 
  mm = length(samplvec)
  if(bw <= 0){     
    bw = 1.06*sd(samplvec)/mm^0.2
  }
  outdens = rep(0, length(densx))
  for(j in 1:mm){
    xxs = (densx-samplvec[j])/bw
    outdens = outdens+dnorm(xxs)
  }
  outdens = outdens/(mm*bw)
  return(outdens)
}