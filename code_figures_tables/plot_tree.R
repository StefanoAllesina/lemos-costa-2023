#rm(list = ls())
#devtools::install_github("YuLab-SMU/ggtree")
library(tidyverse)
library(ape)
library(ggtree)
library(phytools)
library(stringr)

get_tree <- function(V, pars){
  # Order the columns of V so that they have consecutive ones
  cst <- apply(V, 2, paste, collapse = "")
  V <- V[,order(cst,decreasing = TRUE)]
  A <- t(V[-1,]) %*% diag(pars[-1]) %*% V[-1,]
  tree <- ape::vcv2phylo(A)
  return(reorder(tree))
}

plot_tree <- function(resultsfile, scalefac, label){
load(resultsfile)
V <- as.matrix(results$V)
n <- ncol(V)
prs <- results$pars_orig[1:(2 * n - 1)]
scaleint <- max(t(V) %*% diag(prs) %*% V)
prs <- prs / scaleint
tree <- get_tree(V, prs^(1/scalefac))
# shorten latin names
for (i in 1:length(tree$tip.label)){
  nm <- tree$tip.label[i]
  nm <- str_replace(nm, "_", " ")
  init <- str_extract(nm, "[A-Z]")
  spname <- str_extract_all(nm, "\\s\\w+")
  tree$tip.label[i] <- paste0(init, ".", spname, " ")
}
pl <- ggtree(tree, aes(colour = scalefac * log(branch.length)), 
             size = 1.5, layout = "dendrogram") + 
  geom_tiplab(colour = "black", size = 2.5, geom = "text") + 
  geom_text(aes(x = branch, label=labels), hjust=-.3) +
  scale_colour_gradientn("log(strength)", colours = c("#051923",
    "#003554",  "#006494", "#0582ca", "#00a6fb")) + 
  theme(legend.position = "none") + ggtitle(label)
  return(pl)
}

pl1 <- plot_tree("../organized_results/cadotte_2013_1.RData", 3, "Cadotte (2013)")
pl2 <- plot_tree("../organized_results/bioII_2014-8_1.RData", 3, "Biodiversity II (8/2014)")
pl3 <- plot_tree("../organized_results/vr_2_1.RData", 3, "Wageningen (2001)")
