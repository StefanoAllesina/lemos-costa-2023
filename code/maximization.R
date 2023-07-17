tominimize <- function(pars, V, P, Obs, replicates, nsp, modelnum){
  return(-likely_gamma(abs(pars), V, P, Obs, replicates, nsp, modelnum))
}

tominimize_broom <- function(pars, V, P, Obs, replicates, nsp, modelnum){
  # zero the pars for internal branches
  n <- ncol(V)
  pars[2:(n-1)] <- 0
  return(-likely_gamma(abs(pars), V, P, Obs, replicates, nsp, modelnum))
}

hc <- function(pars, V, P, Obs, replicates, nsp, modelnum,
               startpert = 0.1, nn1 = 25, nn2 = 25, myfunc){
  cur_res <- myfunc(pars, V, P, Obs, replicates, nsp, modelnum)
  dev <- startpert
  for (j in 1:nn1){
    tmp_res <- cur_res
    tmp_par <- pars
    for (i in 1:nn2){
      pert_par <- pars * (1 + rnorm(length(pars)) * dev)
      pert_res <- myfunc(pert_par, V, P, Obs, replicates, nsp, modelnum)
      if (pert_res < tmp_res){
        tmp_res <- pert_res
        tmp_par <- pert_par
      }
    }
    if (tmp_res < cur_res){
      pars <- tmp_par
      cur_res <- tmp_res
    }
    dev <- 0.8 * dev
  }
  return(pars)
}

call_optim <- function(pars, V, P, Obs, replicates, nsp, modelnum,
                       maxit = 10000, myfunc){
  tmp <- optim(par = pars, fn = myfunc, V = V, P = P, Obs = Obs, replicates = replicates, nsp = nsp, modelnum = modelnum,
               method = "Nelder-Mead", control = list(trace = FALSE, maxit = maxit))
  tmp <- optim(par = tmp$par, fn = myfunc, V = V, P = P, Obs = Obs, replicates = replicates, nsp = nsp,modelnum = modelnum,
               method = "BFGS", control = list(trace = FALSE, maxit = maxit))
  return(tmp$par)
}

sa <- function(cur_pars, V, P, Obs, replicates, nsp, modelnum,
               b0 = 2, nsteps = 25000, pert = 0.1, myfunc){
  cur_fit <- myfunc(cur_pars, V, P, Obs, replicates, nsp, modelnum)   
  best_pars <- cur_pars
  best_fit <- cur_fit
  changed <- FALSE
  for (i in 1:nsteps){
    accept <- FALSE
    # modify the parameters according to pert
    new_pars <- cur_pars * (1 + rnorm(length(cur_pars), mean = 0, sd = pert ))
    new_fit <- myfunc(new_pars, V, P, Obs, replicates, nsp, modelnum)   
    if (new_fit < cur_fit){
      accept <- TRUE
    } else{
      paccept <- exp(b0 * (cur_fit - new_fit) * log(i + 3.718282))
      if (runif(1) < paccept) {
        accept <- TRUE
      }
    }
    if (accept){
      cur_fit <- new_fit
      cur_pars <- new_pars
      if (cur_fit < best_fit){
        best_fit <- cur_fit
        best_pars <- cur_pars
        changed <- TRUE
      }
    }
    if (i %% 15000 == 0) {
      # cat(best_fit) 
      # cat(" ")
      if (changed == FALSE) {
        cur_pars <- best_pars
        cur_fit <- best_fit
        pert <- 0.9 * pert
      }
      changed <- FALSE
    }
  }
  return(best_pars)
}

plot_results <- function(pars, V, P, Obs, replicates, nsp, modelnum){
  X <- get_prediction(abs(pars), V, P, Obs, replicates, nsp, modelnum)
  plot(as.vector(X[Obs != 0.0]), as.vector(Obs[Obs != 0 ]), log = "xy")
  print(cor(log(as.vector(X[Obs != 0.0])), log(as.vector(Obs[Obs != 0 ]))))
}

simple_search <- function(pars, V, P, Obs, replicates, nsp, modelnum,
                          num_optim){
  if (nrow(V) > ncol(V) + 1){
    # First part: maximize using only broom/star tree
    print("maximization using star tree (100 rounds)")
    for (i in 1:100){
      pars <- sa(pars, V, P, Obs, replicates, nsp, modelnum, pert = 0.25/i, myfunc = tominimize_broom)
      print(c(i, -tominimize_broom(pars, V, P, Obs, replicates, nsp, modelnum)))
    }
    # now fit also internal branches
    n <- ncol(V)
    pars[2:(n-1)] <- rnorm(length(2:(n-1)))*0.1 * max(abs(pars[n+1:(n-1)]))
  }
  # Second part: include internal branches
  print("include internal branches (100 rounds)")
  mypert <- 0.1
  for (i in 1:100){
    pars <- sa(pars, V, P, Obs, replicates, nsp, modelnum, pert = mypert, myfunc = tominimize)
    if (i %% 5 == 0) {
      mypert <- mypert * 0.75
      print(c(i, -tominimize(pars, V, P, Obs, replicates, nsp, modelnum)))
    }
  }
  # third part: numerical optimization
  print("numerical optimization (number of rounds as specified)")
  for (i in 1:num_optim){
    pars <- abs(call_optim(pars, V, P, Obs, replicates, nsp, modelnum, myfunc = tominimize))
    if (i %% 200 == 0) {
      for (j in 1:5){
       pars <- abs(hc(pars, V, P, Obs, replicates, nsp, modelnum, startpert = 0.001/j, nn1 = 20, nn2 = 20, myfunc = tominimize))
       pars <- sa(pars, V, P, Obs, replicates, nsp, modelnum, pert = 0.001/j, myfunc = tominimize)

      }
      print(c(i, -tominimize(pars, V, P, Obs, replicates, nsp, modelnum)))
    }
  }
  return(pars)
}


multi_search <- function(pars, V, P, Obs, replicates, nsp, modelnum,
                        num_search, num_optim){
  cur_pars <- abs(pars)
  cur_likely <- -tominimize_broom(cur_pars, V, P, Obs, replicates, nsp, modelnum)
  print(cur_likely)
  for (i in 1:num_search){
    new_pars <- abs(simple_search(pars, V, P, Obs, replicates, nsp, modelnum, num_optim))
    new_likely <- -tominimize(new_pars, V, P, Obs, replicates, nsp, modelnum)
    if (new_likely > cur_likely){
      cur_likely <- new_likely
      cur_pars <- new_pars
    }
    print(c(i, new_likely, cur_likely))
  }
  return(cur_pars)
}
