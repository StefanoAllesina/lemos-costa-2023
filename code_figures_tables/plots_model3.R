# Plots for model 3
rm(list = ls())
library(ggstance)
library(tidyverse)
library(ggtree)
library(ape)

fn <- "../organized_results/bioII_2008-7_3.RData"
fn_name <- "Biodiversity II (2008) model 3"

get_tree <- function(V, pars){
  # Order the columns of V so that they have consecutive ones
  cst <- apply(V, 2, paste, collapse = "")
  V <- V[,order(cst,decreasing = TRUE)]
  A <- t(V[-1,]) %*% diag(pars[-1]) %*% V[-1,]
  tree <- ape::vcv2phylo(A)
  return(reorder(tree))
}

plot_tree_m3 <- function(resultsfile, scalefac, label){
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

  # Change names of columns of V to match tip names of the tree
  sp_names <- paste0(str_extract(str_replace(colnames(V), "_", " "), "[A-Z]"), ".", 
                     str_extract_all(str_replace(colnames(V), "_", " "), "\\s\\w+"), " ")
  
  df_bar_data <- tibble(id = sp_names,
                        bar_value = results$pars_orig[(2*n):(3*n-1)])
  
  
  pl <- ggtree(tree, aes(colour = scalefac * log(branch.length)), 
               size = 1.5) + 
    geom_tiplab(colour = "black", size = 2.5, geom = "text") + 
    scale_colour_gradientn("log(strength)", colours = c("#051923","#003554",  "#006494", "#0582ca", "#00a6fb")) + 
    geom_facet(panel = expression(theta), data = df_bar_data, geom = geom_barh,stat='identity',
               aes(x = bar_value), width = .6, fill="#6C2691", color = NA) +
    theme(legend.position = "none") + ggtitle(label)
  return(pl)
}


pl_tree <- plot_tree_m3(fn, 3, fn_name)
pl_tree


## Creating all the plots
files <- list.files(path = "../organized_results/", pattern = "\\_3.RData$")
fn <- paste0("../organized_results/", files)
fn_name <- tools::file_path_sans_ext(basename(fn))
fn_name <- unlist(strsplit(sub('(^[^_]+_[^_]+)_(.*)$', '\\1', fn_name), ' '))

for(i in 1:length(fn)){
  pl_tree <- plot_tree_m3(fn[i], 3, paste0(fn_name[i], " model 3"))
  ggsave(filename = paste0("figures_models_2-4/", fn_name[i], "_m3.pdf"), plot = pl_tree)
  print(i)
}
