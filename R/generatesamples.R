#'
#' \code{generatesamples} This function ...
#'
#' @export

generatesamples <- function(distfunc, densx, trrs,bw=0,Shift=0){
  mm<- length(trrs)
  fremat<- vector('list', mm)
  ctable<- matrix(0, mm, length(densx))
  dtable<- matrix(0, mm,length(densx))
  for(i in 1:mm){
    temp <- runif(trrs[i])  
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
    temp <- densestmm(densx, densx[fremat[[i]]], bw)
    dtable[i,]<-temp
  }
  ########################################
  dtablex <- densx
  rlist<- list(fremat=fremat, ctable=ctable, dtable=dtable,dtablex=dtablex)
  #   browser()
  return(rlist)
}
