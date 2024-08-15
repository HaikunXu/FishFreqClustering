#'
#' \code{find_children} This function ...
#'
#' @export

find_children <- function(nodes, node_id){
  
  node_cells <- teststat$childnodes[node_id][[1]]
  node_children <- rep(NA, 2)
  
  # left child
  for (i in seq(node_id - 1, 1, -1)) {
    if(head(teststat$childnodes[i][[1]],1) == head(node_cells,1)) {
      node_children[1] <- i
      break
    }
  }
  
  # right child
  for (i in seq(node_id - 1, 1, -1)) {
    if(tail(teststat$childnodes[i][[1]],1) == tail(node_cells,1)) {
      node_children[2] <- i
      break
    }
  }
  
  return(node_children)
}