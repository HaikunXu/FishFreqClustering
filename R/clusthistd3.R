#'
#' \code{clusthistd3} This function ...
#'
#' @export

clusthistd3 <- function(kk, colseq, textchar, colcol=rep(c(2,3,5 , 4, 6,1,7,8),3),ylims=c(0,3.6)) {
  tcol = colseq
  plot(c(0,1,2),c(0,1,4),xlim=c(0,2), ylim=ylims,pch='  ', 
       ylab='density functions',xlab=' ') 
  tempmat =diag(rrs)%*%as.matrix(densmaty)
  for(ii in 1:kk){
    if(length(tempmat[tcol==(ii+1),]) < 2*ncol(tempmat)){ 
      avrdensy = tempmat[tcol==ii+1,]/rrs[tcol==(ii+1)]  #corrected 2023.8.17
      #               browser()
    }else{
      avrdensy = apply(tempmat[tcol==(ii+1),],2,sum)/sum(rrs[tcol==(ii+1)])
    }
    #           browser()
    lines(densmatx[1,], avrdensy, col=colcol[ii], lwd=3, lty=1)
    text(densmatx[1,which(avrdensy ==max(avrdensy))], max(avrdensy)+0.8, 
         paste(sum(rrs[tcol==(ii+1)])) , col=colcol[ii],cex=1.4)
    text(densmatx[1,which(avrdensy ==max(avrdensy))], max(avrdensy)+0.3, 
         paste('(',sum(rep(1,length(tcol))[tcol==(ii+1)]),')') , col=colcol[ii],cex=1.4)
  }
  title(textchar)
}