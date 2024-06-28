#'
#' \code{drawcells2} This function ...
#'
#' @export

drawcells2 <- function(lats, lons, colseq,putnum = TRUE){
  plot(lons-2.5, lats+2.5 ,col=colseq,pch=15,cex=4,
       ylim=c(-30,30),xlim=c(-150,-70), 
       ylab='latitude', xlab='longitude')#,main=paste('クラスター数 = ',kk))
  if(putnum) text(lons-2.5,lats+2.5,seq(1,length(lats), 1))
  for(ss in seq(-150, -70, 10)) abline(v=ss, lty=2)
  for(ss in seq(-30, 30, 10)) abline(h=ss, lty=2)
}