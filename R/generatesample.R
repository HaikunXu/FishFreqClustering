#'
#' \code{generatesample} This function ...
#'
#' @export

generatesample <- function(distfunc, densx, trrs,nsize=200,bw,Shift=0){
  mm<- length(trrs)
  fremat<- vector('list', mm)
  ctable<- matrix(0, mm, length(densx))
  dtable<- matrix(0, mm, 420)
  for(i in 1:mm){
    temp <- runif(trrs[i]*nsize)  
    if(Shift > 0){
      distfunct = distfunc[(1+Shift):(length(distfunc)-Shift)]
      fremat[[i]]<-cut(temp, breaks=distfunct,labels=FALSE)+Shift
      #         browser()
      ctable[i,] <- table(c(1:length(densx), fremat[[i]]))-1
      #         browser()
    }else{
      fremat[[i]]<-cut(temp, breaks=distfunc,labels=FALSE)
      ctable[i,] <- table(c(1:length(densx), cut(temp, breaks=distfunc,labels=FALSE)))-1
    }
    temp <- density(densx, weights =ctable[i,]/sum(ctable[i,]),bw)
    
    ii <- which(temp$x>=min(densx) & temp$x<=max(densx))

    dtablex <- temp$x[ii]
    dtable[i,]<-temp$y[ii]
  }

  rlist<- list(fremat=fremat, ctable=ctable, dtable=dtable,dtablex=dtablex)
  #   browser()
  return(rlist)
}
