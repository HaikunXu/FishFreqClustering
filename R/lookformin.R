#'
#' \code{lookformin} This function ...
#'
#' @export

lookformin <- function(distantmat, arukaind,adj, adjmat){
  nn <- ncol(distantmat)
  ijpair <-  c(0,0)
  mindistance <- max(abs(distantmat))*10
  for(i in 1:(nn-1)){
    if(arukaind[i] != 0){
      for(j in (i+1):nn){
        if((adj == FALSE) || (adjmat[i,j] == 1)){
          if((arukaind[j]!=0) &&(distantmat[i,j] < mindistance)){
            ijpair <-  c(i, j)
            mindistance <- distantmat[i,j]
          } 
        }
      }
    }
  }
  #    browser()
  if(ijpair[2] == 0) browser()
  return(list(ijpair = ijpair, mindistance = mindistance))
}