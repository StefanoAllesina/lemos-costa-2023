# Subtree plots for three selected datasets -- figure main text
rm(list = ls())
library(tidyverse)
library(ggtree)
library(ape)

## Color gradients
d13 <- colorRampPalette( c( "#FFAB14" , "#622285" ) )( 100 )
d12 <- colorRampPalette( c( "#FFAB14" , "#CC0000" ) )( 100 )
d23 <- colorRampPalette( c( "#622285" , "#CC0000" ) )( 100 )

get_tree <- function(V, pars){
  # Order the columns of V so that they have consecutive ones
  cst <- apply(V, 2, paste, collapse = "")
  V <- V[,order(cst,decreasing = TRUE)]
  A <- t(V[-1,]) %*% diag(pars[-1]) %*% V[-1,]
  tree <- ape::vcv2phylo(A)
  return(reorder(tree))
}

get_lookup <- function(V){
  # associate a number to each possible config
  coderow <- apply(V, 1, paste, collapse = "")
  return(coderow)
}

res_boxplot <- function(results_file){
  load(results_file)
  V_orig <- results$V
  n <- ncol(V_orig)
  pars <- abs(results$pars_orig[1:(2 * n - 1)])
  # set root and tips to zero
  pars[1] <- 0
  pars[n:(2 * n -1)] <- 0
  # compute mean
  mu <- mean(pars[pars > 0])
  # set a label to all pars > mu
  labels <- rep("", length(pars))
  nlarge <- sum(pars > mu)
  labels[pars > mu] <- letters[1:nlarge]
  
  highlight_V <- get_lookup(V_orig)
  highlight_V <- highlight_V[labels != ""]
  highlight_labs <- labels[labels != ""]
  
  forbox <- tibble(logLik = results$lik_orig, lab = "original_lik")
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
  return(list(forbox = forbox, highlights = list(V = highlight_V, labs = highlight_labs)))
}

plot_tree <- function(resultsfile, scalefac, label, phylo_groups, vjust_label, highlights){
  load(resultsfile)
  V <- as.matrix(results$V)
  n <- ncol(V)
  prs <- results$pars_orig[1:(2 * n - 1)]
  scaleint <- max(t(V) %*% diag(prs) %*% V)
  prs <- prs / scaleint
  tree <- get_tree(V, prs^(1/scalefac))
  highlight_V <-  highlights$V
  highlight_labs <-  highlights$labs
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
  my_new_colors[!is.na(d$label)] <- my_colors[2:(phylo_groups+1)][order(d$label[!is.na(d$label)])] # we have to fix this depending on # of groups
  pl <- ggtree(tree, aes(colour = scalefac * log(branch.length)), 
               size = 1.5, layout = "dendrogram") + 
    geom_tiplab(colour = "black", size = 2.7, geom = "text") + 
    geom_text(aes(label=label), data = d, hjust=-.5, vjust = vjust_label, color = my_new_colors, size = 5)+  # also need to fix position
    scale_colour_gradientn("log(strength)", colours = c("#051923","#003554",  "#006494", "#0582ca", "#00a6fb")) + 
    theme(legend.position = "none") + ggtitle(label)
  return(pl)
}


#__________________________________________________________________________#
#               PLOTS
#__________________________________________________________________________#
## Cadotte
cadotte <- "../organized_results/cadotte_2013_1.RData"
fn_name <- "Cadotte (2013)"

res_fn <- res_boxplot(cadotte)
log_Lik <- res_fn$forbox$logLik[res_fn$forbox$lab=="original_lik"]
forbox <- res_fn$forbox %>% filter(lab != "original_lik")
# reorder the labels
uniquelabs <- unique(forbox$lab)
uniquelabs <- uniquelabs[order(nchar(uniquelabs))]
forbox$lab <- factor(forbox$lab, levels = c("", "a", "b", "ab"))

my_colors <- c("#F9E0D6", d13[1], d13[100], d13[50], d23[50])

bp_cadotte <- ggplot(forbox) + aes(x = lab, y = logLik, fill = lab) + 
  geom_boxplot() +
  #geom_point(alpha = 0.2) +
  geom_hline(yintercept = log_Lik, lty = 2) +
  theme_bw() +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none")


pl_tree_cadotte <- plot_tree(cadotte, 3, fn_name, phylo_groups = 2, vjust_label =-14, highlights = res_fn$highlights)

pl_cadotte <- ggpubr::ggarrange(pl_tree_cadotte, bp_cadotte, ncol = 2)



##########################################################
## Van Ruijven
vr2001 <- "../organized_results/vr_4_1.RData"
fn_name <- "Wageningen (2003)" # vr1 refers to year 2000

res_fn <- res_boxplot(vr2001)
log_Lik <- res_fn$forbox$logLik[res_fn$forbox$lab=="original_lik"]
forbox <- res_fn$forbox %>% filter(lab != "original_lik")
# reorder the labels
uniquelabs <- unique(forbox$lab)
uniquelabs <- uniquelabs[order(nchar(uniquelabs))]
forbox$lab <- factor(forbox$lab, levels = c("", "a", "b", "ab"))

bp_vr2001 <- ggplot(forbox) + aes(x = lab, y = logLik, fill = lab) + 
  geom_boxplot() +
  #geom_point(alpha = 0.2) +
  geom_hline(yintercept = log_Lik, lty = 2) +
  theme_bw() +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none")


pl_tree_vr2001 <- plot_tree(vr2001, 3, fn_name, phylo_groups = 2, vjust_label =-14, highlights = res_fn$highlights)

pl_vr <- ggpubr::ggarrange(pl_tree_vr2001, bp_vr2001, ncol = 2)



##########################################################
## Biodiversity II
bioII2001 <- "../organized_results/bioII_2012-7_1.RData"
fn_name <- "Biodiversity II (2012)"

res_fn <- res_boxplot(bioII2001)
log_Lik <- res_fn$forbox$logLik[res_fn$forbox$lab=="original_lik"]
forbox <- res_fn$forbox %>% filter(lab != "original_lik")
# reorder the labels
uniquelabs <- unique(forbox$lab)
forbox$lab <- factor(forbox$lab, levels = c("", "a", "b", "c", "ac", "bc"))

my_colors <- c("#F9E0D6", d13[1], d13[100], d12[100], d12[50], d23[50])

bp_bioII2001 <- ggplot(forbox) + aes(x = lab, y = logLik, fill = lab) + 
  geom_boxplot() +
  #geom_point(alpha = 0.2) +
  geom_hline(yintercept = log_Lik, lty = 2) +
  theme_bw() +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none")


pl_tree_bioII2001 <- plot_tree(bioII2001, 3, fn_name, phylo_groups = 4, vjust_label =-15.5, highlights = res_fn$highlights)

pl_bioII <- ggpubr::ggarrange(pl_tree_bioII2001, bp_bioII2001, ncol = 2)


#ggpubr::ggarrange(pl_cadotte, pl_vr, pl_bioII, nrow = 3)
