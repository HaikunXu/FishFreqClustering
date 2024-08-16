#'
#' \code{find_clusters} This function ...
#'
#' @export

find_clusters <- function(MJS_statistics, distance_threshold, catch_threshold){
  
  # put the first node in
  node_record <- max(MJS_statistics$Node_Number)
  
  for (i in 1:100) {
    flag <- 0
    for (j in 1:length(node_record)) {
      node_id <- node_record[j]
      if(MJS_statistics$STD_distance[node_id] >= distance_threshold & MJS_statistics$Catch_Proportion[node_id] >= catch_threshold) {
        # find children
        node_children <- find_children(node_id)
        if(sum(is.na(node_children))==0) {
          # add children nodes
          node_record <- c(node_record, node_children)
          # remove parent node
          node_record <- node_record[-j]
          
          flag <- 1
          break
        }
      }
    }
    
    if(j == length(node_record) & flag == 0) break
  }
  
  return(node_record)
  
}