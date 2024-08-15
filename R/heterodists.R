#'
#' \code{heterodists} This function ...
#'
#' @export

heterodists <- function(merges,distseq,densx,densmaty,rrs, KL='KL',BB=c(100, 1000),
                        bounds = c(0.005, 0.1)) {
  #
  #      Modifyied from heterodist                                       2024.2
  #
  #     Performing hierarchical test
  #                                              2020.9
  #     Input 
  #         merges : merges component of hclust.regionsmm output(m-1 by 2)
  #         distseq : distseq component of hclust.regionsmm output(m-1)
  #         densx : values correspond to densmaty (distributionsize)
  #         densmaty : matrix of individual distributions   (mm by distributionsize)
  #         rrs : sample sizes of individual distributions  (mm)
  #         KL : 'KL'
  #         BB : the number to repeat clustering at each node BB[1] first test BB[2] second test
  #         bounds : If the standardized value < bounds[1] or the value > bounds[2], stop after the first test
  #stop the whole test if the value > bounds[3]
  # 
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
  distmat <- matrix(0, mm-1, BB[1])
  statmat<- matrix(0, mm-1, 7)
  dimnames(statmat)[[2]]<- c('data dist.','sim dist. mean', 'sim dist. var', 'stdized v','sim dist. max',
                             'out of 100', 'out of 1000')
  #   browser()
  statmat[,1]<- distseq
  ifterm <- rep(0, mm-1)
  
  for(ii in seq(mm-1, 1, -1) ){
    # 
    #   computation for parent node
    #
    if(ifterm[ii] == 0){
      trrs = rrs[childnodes[[ii]]]
      tdistmat = densmaty[childnodes[[ii]],]
      overall_pmf <- as.vector((trrs%*%tdistmat)/sum(trrs))
      overall_pdf  <- cumsum(overall_pmf)
      overall_pdf <- overall_pdf/overall_pdf[length(overall_pdf)]
      converted_pdf<- convpdf(overall_pdf, densx, densx)
      #
      #       computation of distance measures
      #
      for(j in 1:BB[1]){
        shiftj = 0
        if( sum(ifelse(diff(converted_pdf)==0, 1,0)) > 0 ) {
          shiftj = max(seq(1, length(converted_pdf)-1)*ifelse(diff(converted_pdf)==0, 1,0))
        }
        tsamples <- generatesamples(converted_pdf, densx, trrs, bw=0,Shift = shiftj)  
        tclust <- hclust.regionsmm(tsamples$dtable, trrs, KL=KL)
        distmat[ii,j] <- tclust$distseq[length(tclust$distseq)]
      }
      #- c('data dist.','sim dist. mean', 'sim dist. var', 'stdized v','sim dist. max',  'out of BB[1]', 'out of BB[2]')
      statmat[ii,2] <- mean(distmat[ii,])
      statmat[ii,3] <- var(distmat[ii,])
      statmat[ii,4] <- (statmat[ii,1]-statmat[ii,2])/sqrt(statmat[ii,3])
      statmat[ii,5] <- max(distmat[ii,])
      statmat[ii,6] <- sum(ifelse(distmat[ii,]>statmat[ii,1], 1,0))
      #      browser()
      
      cheb = 1/statmat[ii,4]^2
      pest = statmat[ii,6]/BB[1]
      if(cheb < bounds[1]){ 
        #                  ii = ii-1
        print(' cheb <  bounds[1] --> reject H0 \n')
      }else if( pest > bounds[2]){
        #                  ii = 0
        if(merges[ii,1]>0) ifterm[merges[ii,1]] <- 1
        if(merges[ii,2]>0) ifterm[merges[ii,2]] <- 1
        print('pvalue  > bounds[2] --> accept H0 \n')
      }else{
        print('between bounds\n')
        temp = rep(0, BB[2]) 
        for(j in 1:BB[2]){
          #                     print(paste('ii=',ii,'j=',j))   #20200905
          shiftj = 0
          if( sum(ifelse(diff(converted_pdf)==0, 1,0)) > 0 ) {
            shiftj = max(seq(1, length(converted_pdf)-1)*ifelse(diff(converted_pdf)==0, 1,0))
          }
          tsamples <- generatesamples(converted_pdf, densx, trrs, bw=0,Shift = shiftj)  
          tclust <- hclust.regionsmm(tsamples$dtable, trrs, KL=KL)
          temp[j] <- tclust$distseq[length(tclust$distseq)]
        }
        statmat[ii,7] <- sum(ifelse(temp>statmat[ii,1], 1,0))
      }
    }else{      # node is a homogeneous cluster.
      if(merges[ii,1]>0) ifterm[merges[ii,1]] <- 1
      if(merges[ii,2]>0) ifterm[merges[ii,2]] <- 1
    }
    #      browser()
  }
  return(list(statmat=statmat, distmat=distmat, childnodes=childnodes))
}
