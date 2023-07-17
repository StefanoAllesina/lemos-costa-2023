rm(list=ls())
library(tidyverse)
library(ggtree)
library(ape)

# Some parts of the plot need to be "mannually" fixed -- like the colors and the position of the labels in the phylogeny  
#fn <- "../organized_results/cadotte_2013_1.RData"
#fn <- "../organized_results/vr_9_1.RData"
fn <- "../organized_results/bioII_2008-7_1.RData"

#fn_name <- tools::file_path_sans_ext(basename(fn))
#fn_name <- "Cadotte (2013)"
#fn_name <- "Wageningen (2008)" # vr1 refers to year 2000
fn_name <- "Biodiversity II (2008)"


load(fn)
V_orig <- results$V
n <- ncol(V_orig)
pars <- abs(results$pars_orig[1:(2 * n - 1)])
# set root and tips to zero
pars[1] <- 0
pars[n:(2 * n -1)] <- 0
# compute mean
mu <- mean(pars[pars > 0])
#mu <- 0.001
# set a label to all pars > mu
labels <- rep("", length(pars))
nlarge <- sum(pars > mu)
labels[pars > mu] <- letters[1:nlarge]

get_lookup <- function(V){
  # associate a number to each possible config
  coderow <- apply(V, 1, paste, collapse = "")
  return(coderow)
}

get_tree <- function(V, pars){
  # Order the columns of V so that they have consecutive ones
  cst <- apply(V, 2, paste, collapse = "")
  V <- V[,order(cst,decreasing = TRUE)]
  A <- t(V[-1,]) %*% diag(pars[-1]) %*% V[-1,]
  tree <- ape::vcv2phylo(A)
  return(reorder(tree))
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
#forbox$lab <- factor(forbox$lab, levels = uniquelabs)
forbox$lab <- factor(forbox$lab, levels = c("", "a", "b", "c", "ab", "ac", "bc", "abc"))
#forbox$lab <- factor(forbox$lab, levels = c("", "a", "b", "c", "ab", "ac", "bc", "abc"))
#forbox$lab <- factor(forbox$lab, levels = c("", "b", "c", "d", "bd", "cd")) # bioII 2006
#forbox$lab <- factor(forbox$lab, levels = c("", "c", "d", "cd")) # bioII 2014
#forbox$lab <- factor(forbox$lab, levels = c("", "a", "b", "c", "d", "e", "de", "ae", "be")) # bio2008
#forbox$lab <- factor(forbox$lab, levels = c("", "b", "c", "d", "e", "bc", "bd", "cd", "ce", "de")) # bio2011 greater threshold
## Color gradients
d13 <- colorRampPalette( c( "#FFAB14" , "#622285" ) )( 100 )
d12 <- colorRampPalette( c( "#FFAB14" , "#CC0000" ) )( 100 )
d23 <- colorRampPalette( c( "#622285" , "#CC0000" ) )( 100 )



#my_colors <- c("#F9E0D6", d13[1])
#my_colors <- c("#F9E0D6", d13[1], d13[100], d13[50], d23[50])
my_colors <- c("#F9E0D6", d13[1], d13[100], d12[100], d13[50], d23[50],  d12[50], "#993907")
#my_colors <- c("#F9E0D6", d13[1], d13[100], d12[100], d13[50], d23[50],  d12[50], "#968C20", "#027A96", "#993907")
#names(my_colors) <- levels(forbox$lab)

pl <- ggplot(forbox) + aes(x = lab, y = logLik, fill = lab) + 
  geom_boxplot() +
  #geom_violin() +
  #geom_point(alpha = 0.2) +
  geom_hline(yintercept = results$lik_orig, lty = 2) +
  theme_bw() +
  #scale_fill_linha() +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none")
#show(pl)


## Plot tree
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
  q <- ggtree(tree)
  d <- q$data
  d <- d[!d$isTip,]
  # Change names of columns of V to match tip names of the tree
  sp_names <- paste0(str_extract(str_replace(colnames(V), "_", " "), "[A-Z]"), ".", 
                     str_extract_all(str_replace(colnames(V), "_", " "), "\\s\\w+"), " ")
  # Get the most recent common ancestor of the identified branches
  for (k in 1:length(highlight_labs[highlight_labs %in% uniquelabs])) {
    sp_highlight <- as.logical(as.numeric(strsplit(highlight_V[k], "")[[1]]))
    selected_clade <- sp_names[sp_highlight]
    node_mcra <- getMRCA(tree, selected_clade)
    
    d$label[d$node==node_mcra] <- highlight_labs[highlight_labs %in% uniquelabs][k]
  }
  
  my_new_colors <- rep("#F9E0D6", length(d$label))
  my_new_colors[!is.na(d$label)] <- my_colors[2:6][rank(d$label[!is.na(d$label)])] # we have to fix this depending on # of groups
  pl <- ggtree(tree, aes(colour = scalefac * log(branch.length)), 
               size = 1.5, layout = "dendrogram") + 
    geom_tiplab(colour = "black", size = 2.5, geom = "text") + 
    geom_text(aes(label=label), data = d, hjust=-.5, vjust = -21.5, color = my_new_colors, size = 4)+  # also need to fix position
    scale_colour_gradientn("log(strength)", colours = c("#051923","#003554",  "#006494", "#0582ca", "#00a6fb")) + 
    theme(legend.position = "none") + ggtitle(label)
  return(pl)
}


pl_tree <- plot_tree(fn, 3, fn_name)
#show(pl_tree)

#cadotte_arranged <- ggpubr::ggarrange(pl_tree, pl, ncol = 2)
#ggsave("subtree_plots/fig_subtree_cadotte.pdf", plot = cadotte_arranged, width = 5, height = 4)

#vr_arranged <- ggpubr::ggarrange(pl_tree, pl, ncol = 2)
#ggsave("subtree_plots/fig_subtree_vr2000.pdf", plot = vr_arranged, width = 5, height = 4)
 
bio2_arranged <- ggpubr::ggarrange(pl_tree, pl, ncol = 2)
ggsave("subtree_plots/fig_subtree_bioII_2008_bp.pdf", plot = bio2_arranged, width = 5, height = 4)

