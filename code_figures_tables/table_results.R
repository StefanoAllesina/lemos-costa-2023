rm(list = ls())

# BUILD TABLE RESULTS
# BUILD TIBBLES FOR PLOTTING

library(tidyverse)

get_result <- function(fn){
  load(fn)
  # Part 1: data summary
  nsp <- nrow(results$Obs)
  # number of unique communities
  ncomm <- ncol(results$P)
  # number of measurements
  nmeas <- sum(results$Obs!=0)
  # number of parameters
  npars <- ifelse(results$modelnum == 1, 2 * nsp, 
                  ifelse(results$modelnum == 4, 4 * nsp - 2, 3 * nsp - 1))
  # Part 2: p-values
  likelihoods_scramble <- c(results$lik_orig, results$lik_scr)
  likelihoods_rnd <- c(results$lik_orig, results$lik_rnd)
  pval_scramble <- mean(likelihoods_scramble[-1] >= likelihoods_scramble[1])
  pval_rnd <- mean(likelihoods_rnd[-1] >= likelihoods_rnd[1])
  # Part 3: correlation
  observed <- as.vector(results$Obs)
  predicted <- as.vector(results$Pred)
  predicted <- predicted[observed > 0]
  observed <- observed[observed > 0]
  correlation <- cor(log(observed), log(predicted))
  
  return(tibble(label = results$label, 
                nsp = nsp,
                ncomm = ncomm,
                nmeas = nmeas,
                model = results$model, 
                likelihood = results$lik_orig,
                npars = npars,
                pval_rnd = pval_rnd, 
                pval_scramble = pval_scramble, 
                corr = correlation, 
                alpha = results$pars_orig[length(results$pars_orig)],
                nrnd = length(results$lik_rnd), nscr = length(results$lik_scr)))
}

tb <- tibble()
for (i in list.files("../organized_results/", pattern = ".*RData")){
  tmp <- get_result(paste0("../organized_results/", i))
  print(tmp)
  tb <- rbind(tb, tmp)
}

# change names
tb <- tb %>% inner_join(read_csv("lookup_names.csv")) %>% select(-label) %>% rename(label = short)

write_csv(tb, file = "table_results.csv")


fn <- list.files("../organized_results/", pattern = ".*RData")
observations <- tibble()
likelihoods <- tibble()

# load and process all models
for (i in 1:length(fn)){
  load(paste0("../organized_results/", fn[i]))
  
  # Part 1: Observed vs Predicted
  Obs <- as_tibble(t(results$Obs))
  # assign community
  Obs <- Obs %>% mutate(comm = apply((as.matrix(Obs) > 0) * 1, 1, paste, collapse = ""),
                        experiment = row_number())
  Obs <- Obs %>% mutate(model = results$modelnum, label = results$label)
  Obs <- Obs %>% pivot_longer(names_to = "species", values_to = "observed", 
                              cols = -c(model, label, comm, experiment)) %>% filter(observed !=0)
  Obs <- Obs %>% group_by(model, label, comm, species) %>% 
    mutate(mean_observed = exp(mean(log(observed)))) %>% ungroup() # geometric mean
  Pred <- t(results$Pred)
  colnames(Pred) <- rownames(results$Obs)
  Pred <- as_tibble(Pred)
  # assign community
  Pred <- Pred %>% mutate(comm = apply((as.matrix(Pred) > 0) * 1, 1, paste, collapse = ""),
                          experiment = row_number())
  Pred <- Pred %>% mutate(model = results$modelnum, label = results$label)
  Pred <- Pred %>% pivot_longer(names_to = "species", values_to = "predicted", 
                                cols = -c(model, label, comm, experiment)) %>% filter(predicted !=0)
  observations <- bind_rows(observations, 
                            inner_join(Obs, Pred,
                                       by = join_by(comm, experiment, model, label, species)))
  
  # Part 2: Likelihoods
  all_lik <- c(results$lik_orig, results$broom, results$lik_rnd, results$lik_scr)
  all_lik <- tibble(logLik = all_lik, type = c("original tree",
                                               "broom tree",
                                               rep("random tree", length(results$lik_rnd)),
                                               rep("shuffled tree", length(results$lik_scr))))
  all_lik <- all_lik %>% mutate(model = results$modelnum, label = results$label)
  likelihoods <- bind_rows(likelihoods, all_lik)
}

# change names
observations <- observations %>% inner_join(read_csv("lookup_names.csv")) %>% select(-label) %>% rename(label = short)
likelihoods <- likelihoods %>% inner_join(read_csv("lookup_names.csv")) %>% select(-label) %>% rename(label = short)


save(observations, likelihoods, file = "forplots.RData")