join_nodes <- function(V, current){
  # pick two an merge
  if (nrow(current) > 1){
    to_join <- sample(1:nrow(current), 2)
    joined <- current[to_join[1],] + current[to_join[2],]
    current <- current[-to_join,]
    V <- join_nodes(rbind(joined, V), rbind(current, joined))
    rownames(V) <- NULL
  }
  return(V)
}

get_V_random_tree <- function(n){
  V <- diag(rep(1, n))
  return(join_nodes(V, V))
}
