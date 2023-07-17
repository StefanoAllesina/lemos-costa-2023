library(tidyverse)
fn <- "../organized_results/cadotte_2013_1.RData"
#fn <- "../organized_results/vr_2_1.RData"
#fn <- "../organized_results/bioII_2014-8_1.RData"

load(fn)
V_orig <- results$V
n <- ncol(V_orig)
pars <- abs(results$pars_orig[1:(2 * n - 1)])
pars <- pars * rowSums(V_orig)
# set root and tips to zero
pars[1] <- 0
pars[n:(2 * n -1)] <- 0
# compute mean
pars <- pars / max(pars[pars > 0])
mu <- mean(pars[pars > 0])
#mu <- 10^-4
# set a label to all pars > mu
labels <- rep("", length(pars))
nlarge <- sum(pars > mu)
labels[pars > mu] <- letters[1:nlarge]

get_lookup <- function(V){
  # associate a number to each possible config
  coderow <- apply(V, 1, paste, collapse = "")
  return(coderow)
}

highlight_V <- get_lookup(V_orig)
highlight_V <- highlight_V[labels != ""]
highlight_labs <- labels[labels != ""]

forbox <- tibble()
# random trees
for (i in 1:dim(results$random_V)[3]){
  Vi_coderow <- get_lookup(results$random_V[,,i])
  mylab <- ""
  for (j in 1:length(highlight_V)){
    if (highlight_V[j] %in% Vi_coderow) mylab <- paste0(mylab, highlight_labs[j])
  }
  forbox <- rbind(forbox, tibble(logLik = results$lik_rnd[i], lab = mylab))
}
# shuffled trees
for (i in 1:length(results$lik_scr)){
  Vi_coderow <- get_lookup(V_orig[,results$scrambles[i,]])
  mylab <- ""
  for (j in 1:length(highlight_V)){
    if (highlight_V[j] %in% Vi_coderow) mylab <- paste0(mylab, highlight_labs[j])
  }
  forbox <- rbind(forbox, tibble(logLik = results$lik_scr[i], lab = mylab))
}

# reorder the labels
uniquelabs <- unique(forbox$lab)
uniquelabs <- uniquelabs[order(nchar(uniquelabs))]
forbox$lab <- factor(forbox$lab, levels = uniquelabs)


pl <- ggplot(forbox) + aes(x = lab, y = logLik) + geom_boxplot() + geom_hline(yintercept = results$lik_orig, lty = 2)
show(pl)
