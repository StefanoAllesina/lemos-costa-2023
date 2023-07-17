# Tree distances vs likelihoods
#devtools::install_github("ms609/TreeDist")
library(TreeDist)

fn <- "../organized_results/bioII_2014-8_1.RData"
#fn <- "../organized_results/cadotte_2013_1.RData"
get_tree <- function(V){
  # Order the columns of V so that they have consecutive ones
  cst <- apply(V, 2, paste, collapse = "")
  V <- V[,order(cst,decreasing = TRUE)]
  A <- t(V[-1,]) %*% V[-1,]
  tree <- ape::vcv2phylo(A)
  tree$edge.length <- NULL
  return(reorder(tree))
}

distance_mat <- function(S, V){
  S2 <- t(V) %*% V
  return(mean((S - S2)^2))
}

load(fn)
likelihoods <- results$lik_orig
distances_matrix <- c(0)
trees <- list()
V <- results$V
S <- t(V) %*% V
trees[[1]] <- get_tree(V)
original_tree <- trees[[1]]
# add random trees
for (i in 1:dim(results$random_V)[3]){
  VV <- results$random_V[,,i]
  colnames(VV) <- colnames(V)
  tree <- get_tree(VV)
  lik <- results$lik_rnd[i]
  trees[[length(trees) + 1]] <- tree
  likelihoods <- c(likelihoods, lik)
  distances_matrix <- c(distances_matrix, distance_mat(S, VV))
}
# add shuffled trees
for (i in 1:dim(results$scrambles)[2]){
  VV <- results$V[,results$scrambles[i,]]
  colnames(VV) <- colnames(V)
  tree <- get_tree(VV)
  lik <- results$lik_scr[i]
  trees[[length(trees) + 1]] <- tree
  likelihoods <- c(likelihoods, lik)
  distances_matrix <- c(distances_matrix, distance_mat(S, VV))
}

# compute distances
distances <- numeric(0)
for (tt in trees) distances <- c(distances, TreeDistance(original_tree, tt))

toplot <- tibble(distance = distances, distance_mat = distances_matrix, logLik = likelihoods - likelihoods[1])
toplot <- toplot[-1,]
ggplot(toplot) + aes(x = distance, y = logLik) + geom_point(alpha = 0.3) + geom_smooth()
ggplot(toplot) + aes(x = distance_mat, y = logLik) + geom_point(alpha = 0.3) + geom_smooth() + scale_x_sqrt()
