#'
#' \code{clusthistd3} This function ...
#'
#' @export

clusthistd3 <- function(kk, colseq, colcol=rep(c(2,3,5,4,6,1,7,8),3), plot = FALSE) {
  tcol = colseq
  if(plot == TRUE) {
    plot(c(0,1,2), c(0,1,4), pch=' ', ylab='Density functions', xlab = "Meter") 
  }

  tempmat = diag(rrs) %*% as.matrix(densmaty)
  
  for(ii in 1:kk){
    if(length(tempmat[tcol==(ii+1),]) < 2*ncol(tempmat)){ 
      avrdensy = tempmat[tcol==ii+1,]/rrs[tcol==(ii+1)]  # corrected 2023.8.17
      #               browser()
    }else{
      avrdensy = apply(tempmat[tcol==(ii+1),],2,sum)/sum(rrs[tcol==(ii+1)])
    }
    
    if(ii == 1) densy_df <- data.frame("Length" = densmatx[1,], "Density" = avrdensy, "Cell" = ii)
    else densy_df <- rbind(densy_df,
                           data.frame("Length" = densmatx[1,], "Density" = avrdensy, "Cell" = ii))
    
    #           browser()
    if(plot == TRUE) {
      lines(densmatx[1,], avrdensy, col=colcol[ii], lwd=3, lty=1)
      text(densmatx[1,which(avrdensy ==max(avrdensy))], max(avrdensy)+0.2, 
         ii , col=colcol[ii],cex=2)
    }
  }

  densy_df$Cell <- as.factor(densy_df$Cell)
  return(densy_df)
}