library(abind) # to combine arrays
library(tidyverse)

# Combine all randomizations and remove duplicates

# get stems of file names
files_results <- list.files("../results/", pattern = "RData")
files_and_data <- tibble(file = files_results, 
                         data_model = str_replace_all(string = files_results, replacement = "", pattern = "_\\d+\\.RData"))

unique_data_sets <- sort(unique(files_and_data$data_model))

# load broom trees
unique_broom_trees <- unique_data_sets[grepl(pattern = "*._1.$", x = unique_data_sets)]
unique_data_sets <- setdiff(unique_data_sets, unique_broom_trees)
broom_trees <- tibble()
for (i in unique_broom_trees){
  tmp <- files_and_data$file[files_and_data$data_model == i]  
  load(paste0("../results/", tmp))
  broom_trees <- rbind(broom_trees, 
                       tibble(label = results$label, 
                              modelnum = results$modelnum - 10,
                              likelibroom = results$lik_orig))
}

for (i in 1:length(unique_data_sets)){
#i <- 1
current_data_model <- unique_data_sets[i]
randomizations <- files_and_data$file[files_and_data$data_model == current_data_model]
# load the first randomization
load(paste0("../results/", randomizations[1]))
combined_results <- results
# get likelihood of broom tree
likelibroom <- broom_trees$likelibroom[(broom_trees$label == results$label) & (broom_trees$modelnum == results$modelnum)]
combined_results$broom <- likelibroom
# add the other randomizations
for (i in 2:length(randomizations)){
  load(paste0("../results/", randomizations[i]))  
  # check whether to update original fit
  if (results$lik_orig > combined_results$lik_orig){
    combined_results$lik_orig <- results$lik_orig
    combined_results$Pred <- results$Pred
    combined_results$pars_orig <- results$pars_orig
    print(c(combined_results$label, combined_results$modelnum, combined_results$lik_orig))
    
  }
  # add scrambled trees
  combined_results$lik_scr <- c(combined_results$lik_scr, results$lik_scr)
  combined_results$scrambles <- rbind(combined_results$scrambles, results$scrambles)
  # add random trees
  combined_results$lik_rnd <- c(combined_results$lik_rnd, results$lik_rnd)
  combined_results$random_V <- abind(combined_results$random_V, results$random_V)
}

print(combined_results$lik_orig)
# save the combined results
# now that we have all results remove duplicate topologies
getSigma <- function(V){
  return(paste(as.vector(t(V) %*% V), collapse = ""))
}

# isomorphic trees for scrambles
V <- combined_results$V
labels <- getSigma(V)
likelihoods_scramble <- combined_results$lik_orig
unique_scrambles <- t(as.matrix(1:ncol(V), nrows = 1, ncols = ncol(V)))

for (i in 1:nrow(combined_results$scrambles)){
  sc <- combined_results$scrambles[i,]
  V2 <- V[, sc]
  Sigma <- getSigma(V2)
  my_lik <- combined_results$lik_scr[i]
  for (j in 1:length(labels)){
    if (labels[j] == Sigma){
      if (likelihoods_scramble[j] < my_lik) likelihoods_scramble[j] <- my_lik
      break
    }
  }
  if (j == length(labels)){
    # this is a new tree
    labels <- c(labels, Sigma)
    unique_scrambles <- rbind(unique_scrambles, sc)
    likelihoods_scramble <- c(likelihoods_scramble, my_lik)
  }
}

unique_scrambles <- unique_scrambles[-1,] # remove original
likelihoods_scramble <- likelihoods_scramble[-1]
print("scrambles before filter")
print(length(likelihoods_scramble))
# remove simulations that failed and have likelihoods worse than broom
unique_scrambles <- unique_scrambles[likelihoods_scramble >= trunc(likelibroom, 3),]
likelihoods_scramble <- likelihoods_scramble[likelihoods_scramble >= trunc(likelibroom, 3)]
print("scrambles after filter")
print(length(likelihoods_scramble))
# take 1000 randomizations max
combined_results$scrambles <- unique_scrambles[1:(min(nrow(unique_scrambles), 1000)), ]
combined_results$lik_scr <- likelihoods_scramble[1:(min(length(likelihoods_scramble), 1000))]

# isomorphic trees for random trees
V <- combined_results$V
labels <- getSigma(V)
likelihoods_rnd <- combined_results$lik_orig
unique_V_random <- array(V, dim = c(dim(V), 1))

for (i in 1:dim(combined_results$random_V)[3]){
  V2 <- combined_results$random_V[,,i]
  Sigma <- getSigma(V2)
  my_lik <- combined_results$lik_rnd[i]
  for (j in 1:length(labels)){
    if (labels[j] == Sigma){
      if (likelihoods_rnd[j] < my_lik) likelihoods_rnd[j] <- my_lik
      break
    }
  }
  if (j == length(labels)){
    # this is a new tree
    labels <- c(labels, Sigma)
    likelihoods_rnd <- c(likelihoods_rnd, my_lik)
    unique_V_random <- abind(unique_V_random, V2)
  }
}

unique_V_random <- unique_V_random[,,-1] # remove original
likelihoods_rnd <- likelihoods_rnd[-1]
print("rnd before filter")
print(length(likelihoods_rnd))

# remove simulations that failed and have likelihoods worse than broom
unique_V_random <- unique_V_random[,,likelihoods_rnd >= trunc(likelibroom, 3)]
likelihoods_rnd <- likelihoods_rnd[likelihoods_rnd >= trunc(likelibroom, 3)]
print("rnd after filter")
print(length(likelihoods_rnd))

# take 1000 randomizations max
combined_results$random_V <- unique_V_random[,,1:min(dim(unique_V_random)[3], 1000)]
combined_results$lik_rnd <- likelihoods_rnd[1:min(1000, length(likelihoods_rnd))]

results <- combined_results
save(results, file = paste0("../organized_results/", results$label, "_", results$modelnum, ".RData"))
}