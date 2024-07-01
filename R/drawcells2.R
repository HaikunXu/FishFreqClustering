#'
#' \code{drawcells2} This function ...
#'
#' @export

drawcells2 <- function(lats, lons, dlat, dlon, colseq, putnum = TRUE){
  plot(lons, lats ,col=colseq,pch=15,cex=4,
       xlim=c(min(lons)-dlat,max(lons)+dlon),ylim=c(min(lats)-dlat,max(lats)+dlat), 
       ylab='latitude', xlab='longitude')#,main=paste('クラスター数 = ',kk))
  if(putnum) text(lons,lats,seq(1,length(lats), 1))
  for(ss in seq(min(lons)-dlon/2, max(lons)+dlon/2, dlon)) abline(v=ss, lty=2)
  for(ss in seq(min(lats)-dlat/2, max(lats)+dlat/2, dlat)) abline(h=ss, lty=2)
}