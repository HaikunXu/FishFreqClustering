#'
#' \code{heterodist} This function ...
#'
#' @export

heterodist <- function(merges,distseq,densx,densmaty,rrs, KL='KL',doko=c(1,1),BB=5, bw=0.05, bins) {
  #
  #     Computing standardized distance from  homogenity 
  #                                              2020.3
  #     Input 
  #         merges : merges component of hclust.regionsmm output(m-1 by 2)
  #         distseq : distseq component of hclust.regionsmm output(m-1)
  #         densx : values correspond to densmaty (distributionsize)
  #         densmaty : matrix of individual distributions   (mm by distributionsize)
  #         rrs : sample sizes of individual distributions  (mm)
  #         KL : 'KL'
  #         doko : doko[1] node id from top  doko[2] the number to compute distance
  #         BB : the number to repeat clustering at each node
  # 
  #     memo:
  #        seq(0.00, 2.01, 0.01) is the breakpoints of original data
  #
  mm = length(rrs)
  #########################################
  #     Find member distributions for each node
  #########################################
  childnodes <- vector("list", length=mm-1)
  ichildnodes <- vector("list", length=2)
  for(i in 1:(mm-1)){
    for(j in 1:2){
      if(merges[i,j] < 0){  
        ichildnodes[[j]] <- abs(merges[i,j])
      }else{
        ichildnodes[[j]] <- childnodes[[merges[i,j]]]
      }
    }
    childnodes[[i]] <- c(ichildnodes[[1]],ichildnodes[[2]])     
  }
  
  #########################################
  #     Compute distance criteria for nodes
  #########################################
  distmat <- matrix(0, mm-1, BB)
  statmat<- matrix(0, mm-1, 6)
  dimnames(statmat)[[2]]<- c('data dist.','sim dist. mean', 'sim dist. var', 'stdized v','sim dist. max','data dist./max')
  #   browser()
  statmat[,1]<- distseq
  startnode <- mm - doko[1]
  endnode <- mm-doko[1]-doko[2]+1
  for (i in seq(startnode, endnode, -1)){ 
    print(paste0('node = ',i))   #20200905
    
    # 
    #   computation for parent node
    #
    trrs = rrs[childnodes[[i]]]
    tdistmat = densmaty[childnodes[[i]],]
    overall_pmf <- as.vector((trrs%*%tdistmat)/sum(trrs))
    overall_pdf  <- cumsum(overall_pmf)
    overall_pdf <- overall_pdf/overall_pdf[length(overall_pdf)]

    converted_pdf<- convpdf(overall_pdf, densx, bins)
    
    #
    #    computation of distance measures
    #
    for(j in 1:BB){
      shiftj = 0
      if( sum(ifelse(diff(converted_pdf)==0, 1,0)) > 0 ) {
        shiftj = max(seq(1, length(converted_pdf)-1)*ifelse(diff(converted_pdf)==0, 1,0))
      }
      #        if(i == 8) browser()
      tsamples <- generatesample(converted_pdf, bins, trrs, 200,bw=bw,Shift = shiftj)  
      tclust <- hclust.regionsmm(tsamples$dtable, trrs, KL=KL)
      distmat[i,j] <- tclust$distseq[length(tclust$distseq)]
      #         browser()
    }
    statmat[i,2] <- mean(distmat[i,])
    statmat[i,3] <- var(distmat[i,])
    statmat[i,4] <- (statmat[i,1]-statmat[i,2])/sqrt(statmat[i,3])
    statmat[i,5] <- max(distmat[i,])
    statmat[i,6] <- statmat[i,1]/statmat[i,2]
    #      browser()
  }
  return(list(statmat=statmat, distmat=distmat, childnodes=childnodes))
}
