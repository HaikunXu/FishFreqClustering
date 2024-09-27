#'
#' \code{cplotu} This function ...
#'
#' @export

cplotu <- function(merges, distseq,hopt = 'unit',plotnum=10000,...){
  #                                                   plotnum   # 2020.08.29
  nnum = nrow(merges)+1
  #
  # step 1)  find positions
  #
  xpos = rep(0, nnum)        #x position for leaf
  xgroup = rep(0,nnum)      # group for leaf
  nheight = rep(0, nnum)    # height for node
  nxpos  = rep(0, nnum)               #  x position for node
  ngroup = rep(0, nnum)              # group for node 
  minx = 0    #leaf id from left
  if(hopt == 'dev') cdev = 0
  #  
  maxdist = max(distseq)
  for(k in 1:(nnum-1)){ 
    if(merges[k,1] < 0 && merges[k,2] < 0){
      xpos[abs(merges[k,1])] = 1
      xpos[abs(merges[k,2])] = 2
      xgroup[abs(merges[k,1])] = k
      xgroup[abs(merges[k,2])] = k
      nxpos[k] = 1.5
      if(hopt=='unit'){ 
        nheight[k] = 1
      } else if (hopt=='dev'){
        cdev = cdev+distseq[k]; nheight[k] = cdev
      } else {
        nheight[k]=distseq[k]
      }
      ngroup[k] = k
      #print('case1'); browser()
    }else if(merges[k,1] < 0 && merges[k,2] > 0){
      thisgroup = which(xgroup == ngroup[merges[k,2]])
      xpos[abs(merges[k,1])] = max(xpos[thisgroup])+1
      xgroup[abs(merges[k,1])] = k
      xgroup[thisgroup] = k
      thisngroup = which(ngroup == ngroup[merges[k,2]])
      nxpos[k] = nxpos[merges[k,2]]/2+xpos[abs(merges[k,1])]/2
      if(hopt == 'unit'){
        nheight[k] = nheight[merges[k,2]]+1
      }else if(hopt=='dev'){
        cdev = cdev+distseq[k]; nheight[k] = cdev
      }else{
        nheight[k] = distseq[k]
      }
      ngroup[k] = k
      ngroup[thisngroup] = k  
      #print('case2'); browser()
    }else if(merges[k,1] > 0 && merges[k,2] < 0){
      #       browser()
      thisgroup = which(xgroup == ngroup[merges[k,1]])
      xpos[abs(merges[k,2])] = max(xpos[thisgroup])+1
      xgroup[abs(merges[k,2])] = k
      xgroup[thisgroup] = k
      thisngroup = which(ngroup == ngroup[merges[k,1]])
      nxpos[k]=nxpos[merges[k,1]]/2+xpos[abs(merges[k,2])]/2
      if(hopt =='unit'){
        nheight[k] = nheight[merges[k,1]]+1
      }else if(hopt=='dev'){
        cdev = cdev+distseq[k]; nheight[k] = cdev
      }else{
        nheight[k] = distseq[k]
      }
      ngroup[k]=k
      ngroup[thisngroup] = k
      #print('caes3');browser()
    }else{
      thisgroup1 = which(xgroup == ngroup[merges[k,1]])
      thisgroup2 = which(xgroup == ngroup[merges[k,2]])
      xpos[thisgroup2] = xpos[thisgroup2]+max(xpos[thisgroup1])
      xgroup[thisgroup1] = k
      xgroup[thisgroup2] = k
      ngroup[k] = k
      thisngroup1 = which(ngroup == ngroup[merges[k,1]])
      thisngroup2 = which(ngroup == ngroup[merges[k,2]])
      nxpos[thisngroup2] = nxpos[thisngroup2]+max(xpos[thisgroup1])
      ngroup[thisngroup1] = k
      ngroup[thisngroup2] = k
      nxpos[k] = nxpos[merges[k,1]]/2+nxpos[merges[k,2]]/2
      if(hopt=='unit'){
        nheight[k] = max(nheight[merges[k,1]],nheight[merges[k,2]])+1
      }else if(hopt=='dev'){
        cdev = cdev+distseq[k]; nheight[k] = cdev
      } else{
        nheight[k] = distseq[k]
      } 
      #print('case4');browser()    
    }
  }
  #   browser()
  #
  # plot dendrogram
  #
  plot(c(0, nnum),c(0,max(nheight)),pch=' ',axes=FALSE,xlab=' ',ylab=' ',...)
  for(k in 1:(nnum-1)){
    if(merges[k,1] < 0 && merges[k,2] < 0){
      x1 = xpos[abs(merges[k,1])]; x2 = xpos[abs(merges[k,2])]
      y1 = 0; y2 = nheight[k]; y3 = 0
      #         lines(c(x1,x1,x2,x2), c(0,y2,y2,0)) 
    }else if(merges[k,1] < 0 && merges[k,2] > 0){
      x2 = xpos[abs(merges[k,1])]; x1 = nxpos[merges[k,2]]
      y1=0; y2 = nheight[k];y3=nheight[merges[k,2]]
      #         lines(c(x1,x1,x2,x2), c(0,y2,y2,y3)) 
    }else if(merges[k,1] > 0 && merges[k,2] < 0){
      x2 = xpos[abs(merges[k,2])]; x1 = nxpos[merges[k,1]]
      y1 = nheight[merges[k,1]]; y2 = nheight[k];y3=0
      #         lines(c(x1,x1,x2,x2), c(y1,y2,y2,0)) 
    }else{
      x1 = nxpos[merges[k,1]]; x2 = nxpos[merges[k,2]]
      y1 = nheight[merges[k,1]]; y2 = nheight[k]; y3 = nheight[merges[k,2]]
      #         lines(c(x1,x1,x2,x2), c(y1,max(y1,y2)+1,max(y1,y2)+1,y2)) 
    } 
    lines(c(x1,x1,x2,x2), c(y1, y2, y2, y3))
    if(sum(plotnum==k)==1) text((x1+x2)/2, y2, k,cex=2)       # 2020.08.29
  }
  #   return(list(xpos=xpos, ngroup=ngroup, nxpos=nxpos, nheight=nheight))
  return()
}