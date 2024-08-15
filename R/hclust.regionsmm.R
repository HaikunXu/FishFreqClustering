#'
#' \code{hclust.regionsmm} This function ...
#'
#' @export

hclust.regionsmm <- function(distribmatin,rr,adj=FALSE, adjmat=NULL,KL="KL"){
  distribmat <- distribmatin
  nn <- nrow(distribmat)
  distantmat <- matrix(0, nn,nn) 
  if(is.null(adjmat)){adjmat <- matrix(1, nn,nn)}
  arukaind <- seq(-1, -nn,-1)# 
  # netative => leaf,  0 ,  positive => node
  #     nodetocol <- rep(0, nn)    # pointer for node id and column
  rrn <- rr
  merges <- matrix(0, nn-1,2)
  distseq <- rep(0, nn-1)
  #     cdf and pmf for all population when methods are for cdf
  if(KL != "KL"){
    allcdf = rep(0, ncol(distribmatin))
    for(i in 1:nrow(distribmatin)) allcdf = allcdf+distribmatin[i,]*rr[i]
    allcdf  = allcdf/sum(rr)
    allpmf = c(allcdf[2:length(allcdf)],1)-allcdf
    #      browser()
  }
  # initial value = distance between distributions
  for(i in 1:(nn-1)){
    for(j in (i+1):nn){
      if(KL=="KL"){
        distantmat[i,j] <- kldist(distribmat[i,],rr[i], distribmat[j,],rr[j])
      }else{
        distantmat[i,j] <- vmdist(distribmat[i,],rr[i], distribmat[j,],rr[j],allpmf, KL)
      }
    }
  }
  #   browser()
  for(i in 1:(nn-1)){ 
    minparout <- lookformin(distantmat, arukaind,adj,adjmat)
    distseq[i] <-minparout$mindistance          
    ii <- minparout$ijpair[1]
    jj <- minparout$ijpair[2]
    merges[i, ] <- c(arukaind[ii], arukaind[jj])
    arukaind[ii] <-  i#custer  number
    arukaind[jj] <- 0#
    rr[ii] <- rr[ii]+rr[jj]
    distribmat[ii] <- (rr[ii]*distribmat[ii]+rr[jj]*distribmat[jj])/(rr[ii]+rr[jj])
    for(it in 1:nn){
      if (adjmat[jj, it] > 0.5)  {
        adjmat[ii, it] <- 1
        adjmat[it, ii] <- 1
      }
    } 
    #         browser()
    for(i1 in 1:(ii-1)) {
      if(KL=="KL"){
        distantmat[i1, ii]<- kldist(distribmat[i1,],rr[i1], distribmat[ii,],rr[ii])  
      }else{
        distantmat[i1, ii]<- vmdist(distribmat[i1,],rr[i1], distribmat[ii,],rr[ii],allpmf, KL)  
      }
    }
    for(i1 in (ii+1):nn) {
      if(KL=="KL"){
        distantmat[ii, i1]<- kldist(distribmat[i1,],rr[i1], distribmat[ii,],rr[ii])  
      }else{
        distantmat[ii, i1]<- vmdist(distribmat[i1,],rr[i1], distribmat[ii,],rr[ii],allpmf, KL)  
      }
    }
  }
  #         browser()
  return(list(merges=merges,distseq=distseq,distantmat=distantmat) )
}