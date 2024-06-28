#'
#' \code{putcolor} This function ...
#'
#' @export

putcolor <- function(mergep,kk){
  nn <- nrow(mergep)+1
  temp <- rep(0, nn)    #nodeへの色
  temp2 <- rep(0, nn)  #leaf への色
  #
  coln <- 2
  for(i in (nn-1):(nn-kk+1)){
    for(j in 1:2){
      if(mergep[i,j] <0){#leaf がクラスタートップ
        temp2[abs(mergep[i,j])] <- coln
        coln <- coln + 1
      } else if(mergep[i,j]< (nn-kk+1) ){  #nodeがクラスタートップ
        temp[mergep[i,j]] <- coln
        coln <- coln + 1
      }
    }
  }
  for(i in (nn-kk):1){
    for(j in 1:2){
      if(mergep[i,j] <0){#leaf に色を与える
        temp2[abs(mergep[i,j])] <- temp[i]
        #           browser()
      }else {  #nodeに色を与える
        temp[mergep[i,j]] <- temp[i]
      }
    }  
  }
  #  browser()
  return(temp2)
}
