#'
#' \code{adjinf} This function ...
#'
#' @export

adjinf <- function(lats, lons, mindist = 5){ 
  nn <- length(lats)
  adjinf <- matrix(0, nn, nn)
  mindistsq <- mindist*mindist
  for(i in 1:(nn-1)){
    for(j in (i+1):nn) { 
      rr <- (lats[i]-lats[j])^2+(lons[i]-lons[j])^2
      if(rr < mindistsq*1.001) {
        adjinf[i,j] = 1
        adjinf[j,i] = 1
      }
    }
  }
  return(adjinf)
}