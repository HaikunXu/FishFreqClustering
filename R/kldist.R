#'
#' \code{kldist} This function ...
#'
#' @export

kldist <- function(dist1, rr1, dist2, rr2){
  cmbdist <- (dist1*rr1+dist2*rr2)/(rr1+rr2)
  temp1 <- ifelse(dist1 !=0, dist1*log(dist1), 0)
  temp2 <- ifelse(dist2 !=0, dist2*log(dist2), 0)
  temp3 <- ifelse(cmbdist != 0, cmbdist*log(cmbdist), 0)
  #    browser()
  outvalue <- sum(rr1*temp1+rr2*temp2 -(rr1+rr2)*temp3) 
  return(outvalue)
}