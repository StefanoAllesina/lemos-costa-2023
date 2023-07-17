library(rstan)

expose_stan_functions("../functions/endpoints_including_broom.stan")
source("maximization.R")
source("get_random_V.R")

run_model <- function(file_obs, modelnum, seed = 1, 
                      num_search = 10, num_optim = 25, num_rand = 100){

  set.seed(seed)
  # for each model, we have a corresponding model with the broom tree
  isbroom <- modelnum > 10
  Obs <- t(as.matrix(read.csv(file_obs)))
  V <- as.matrix(read.csv(paste0(tools::file_path_sans_ext(file_obs), "_V.csv")))
  if (isbroom){
    n <- ncol(V)
    V <- V[c(1, (n):nrow(V)),]
  }
  # create compact version of the data to speed up calculations
  P <- (Obs > 0) * 1
  # communities
  comm <- apply(P, 2, paste, collapse = "")
  replicates <- as.integer(factor(comm, levels = unique(comm)))
  # now take only unique communities
  P <- t(unique(t(P)))
  nsp <- colSums(P)
  set.seed(seed)
  # start parameters from perturbed identity matrix
  n <- nrow(P)
  if (modelnum == 1) pars <- c(0.01 * rnorm(n-1), rep(1, n + 1))
  if (modelnum == 2) pars <- c(0.01 * rnorm(n-1), rep(1, 2 * n + 1))
  if (modelnum == 3) pars <- c(0.01 * rnorm(n-1), rep(1, 2 * n + 1))
  if (modelnum == 4) pars <- c(0.01 * rnorm(n-1), rep(1, 3 * n + 1))
  # models with broom tree
  if (modelnum == 11) {
    pars <- c(0.01 * rnorm(1), rep(1, n + 1))
    modelnum <- 1
  }
  if (modelnum == 12) {
    pars <- c(0.01 * rnorm(1), rep(1, 2 * n + 1))
    modelnum <- 2
  }
  if (modelnum == 13) {
    pars <- c(0.01 * rnorm(1), rep(1, 2 * n + 1))
    modelnum <- 3
  }
  if (modelnum == 14) {
    pars <- c(0.01 * rnorm(1), rep(1, 3 * n + 1))
    modelnum <- 4
  }
  # run original V
  pars_orig <- multi_search(pars, V, P, Obs, replicates, nsp, modelnum, num_search, num_optim)
  lik_orig <- likely_gamma(pars_orig, V, P, Obs, replicates, nsp, modelnum)
  X_orig <- get_prediction(pars_orig, V, P, Obs, replicates, nsp, modelnum)
  plot_results(pars_orig, V, P, Obs, replicates, nsp, modelnum)
  likelihoods_scramble <- numeric(0)
  likelihoods_random <- numeric(0)
  scrambles <- matrix(0, 0, ncol(V))
  random_V <- array(0, c(dim(V), num_rand))
  if (!isbroom){
    for (i in 1:num_rand){
      print(paste("rnd scramble", i))
      scramble <- sample(1:ncol(V))
      scrambles <- rbind(scrambles, scramble)
      V2 <- V[,scramble]
      colnames(V2) <- colnames(V)
      pars_rnd <- multi_search(pars, V2, P, Obs, replicates, nsp, modelnum, num_search, num_optim)
      lik_rnd <- likely_gamma(pars_rnd, V2, P, Obs, replicates, nsp, modelnum)
      likelihoods_scramble <- c(likelihoods_scramble, lik_rnd)
      print(sum(likelihoods_scramble < lik_orig) / length(likelihoods_scramble))
      print(paste("rnd Yule", i))
      V2 <- get_V_random_tree(n)
      colnames(V2) <- colnames(V)
      random_V[,,i] <- V2
      pars_rnd <- multi_search(pars, V2, P, Obs, replicates, nsp, modelnum, num_search, num_optim)
      lik_rnd <- likely_gamma(pars_rnd, V2, P, Obs, replicates, nsp, modelnum)
      likelihoods_random <- c(likelihoods_random, lik_rnd)
      print(sum(likelihoods_random < lik_orig) / length(likelihoods_random))
    }
  } else {
    modelnum <- modelnum + 10
  }
  # save results
  label <- basename(tools::file_path_sans_ext(file_obs))
  results <- list(
    label = label,
    modelnum = modelnum,
    lik_orig = lik_orig,
    Obs = Obs,
    Pred = X_orig,
    pars_orig = pars_orig,
    V = V,
    P = P,
    replicates = replicates,
    nsp = nsp,
    lik_scr = likelihoods_scramble,
    lik_rnd = likelihoods_random,
    scrambles = scrambles,
    random_V = random_V
  )
  save(results, file = paste0("../results/", label, "_", modelnum, "_", seed, ".RData"))
  return(results)
}
