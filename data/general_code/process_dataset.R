library(tidyverse)
# devtools::install_github("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)

get_data <- function(dt, spinfo, min_comm = 1, biom_thresh = 0, min_plots = 0, absent_thresh = 0.01){
  print("organizing data")
  # add info communities
  dt <- dt %>% mutate(spcode = as.numeric(as.factor(species)))
  dt <- dt %>% group_by(plot) %>% 
    mutate(comm = paste(sort(spcode), collapse = "-")) %>% 
    ungroup() %>% select(-spcode)
  
  # set biomass to zero for species falling below absent_thresh 
  dt <- dt %>% group_by(species) %>% mutate(maxbio = max(biomass)) %>% 
    ungroup() %>% mutate(biomass = ifelse(biomass < absent_thresh * maxbio, 0.0, biomass)) %>% 
    filter(biomass > 0) %>% select(-maxbio)
  
  # remove plots with too little biomass
  torem <- dt %>% filter(biomass < biom_thresh) %>% select(plot)
  dt <- anti_join(dt, torem, by = "plot")
  
  # remove plots with rare species
  torem <- tibble(aa = 1)
  while(nrow(torem) > 0){
    torem <- dt %>% 
      select(species, comm) %>% distinct() %>% 
      group_by(species) %>% mutate(nc = n()) %>% 
      ungroup() %>% 
      filter(nc < min_comm)
    torem <- torem %>% select(comm)
    dt <- anti_join(dt, torem, by = "comm")
  }
  
  if (nrow(dt) == 0) return(NULL)
  print(dt)
  # pl1 <- ggplot(dt) + aes(x = as.factor(plot), y = biomass, fill = species) + geom_col(position = "dodge") + facet_wrap(~comm, scales = "free") + scale_y_sqrt()
  # show(pl1)
  # pl2 <- ggplot(dt %>% group_by(species, comm) %>% summarise(mm = mean(biomass), vv = var(biomass), .groups = "drop")) + 
  #   aes(x = mm, y = vv, colour = species) + geom_point() + facet_wrap(~species, scales = "free")
  # show(pl2)
  
  # build the tree
  print("building the tree (this can take a minute)")
  # load the info on spp, genus, family (note: take only the species in dt)
  spinfo <- dt %>% select(species) %>% distinct() %>% 
    inner_join(spinfo, by = "species")
  n <- nrow(spinfo)
  tmp <- phylo.maker(as.data.frame(spinfo))
  tree <- tmp$scenario.3
  G <- igraph::graph_from_edgelist(tree$edge)
  Dist <- as.matrix(igraph::distances(G, mode="in"))
  V <- t((!is.infinite(Dist[1:n, 1:(2 * n - 1)])) * 1)
  V <- V[order(rowSums(V), decreasing = TRUE),]
  colnames(V) <- tree$tip.label
  
  # now make the names the same
  print("organizing output")
  dt <- dt %>% mutate(species = str_replace(species, " ", "_"))
  X <- dt %>% pivot_wider(names_from = species, values_from = biomass, values_fill = 0)
  X <- X %>% arrange(comm) %>% select(-plot, -comm) %>% as.matrix()
  # order the cols of the matrix as V
  X <- X[,tree$tip.label]
  # order the rows of the matrix according to the community
  commX <- apply((X > 0) * 1, 1, function(x) paste(x, collapse = ""))
  X <- X[order(commX), ]
  if (nrow(X) < min_plots) return(NULL)
  return(list(X = X, V = V, tree = tree))
}