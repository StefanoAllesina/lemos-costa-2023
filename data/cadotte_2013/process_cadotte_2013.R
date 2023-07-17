source("../general_code/process_dataset.R")

# read the data
dt <- read.csv("rawdata/dryad_data/BIOMASS.csv", stringsAsFactors = FALSE) %>% 
  as_tibble() %>% 
  select(-TREAT, -RICH, -PHY, -LITTER, -ABOVE_BIO) %>% 
  gather(Species, biomass, -PLOT) %>% 
  filter(!is.na(biomass)) %>% filter(biomass > 0)

# change species names to normalize them
dt <- dt %>% 
  rename(code = Species) %>% 
  inner_join(read_csv("Spp_codes.csv") %>% 
               select(code, species), by = "code") %>% 
  rename(plot = PLOT) %>% 
  select(-code) 

spinfo <- read_csv("Spp_codes.csv") %>% select(-code)

out <- get_data(dt, spinfo, min_comm = 5, biom_thresh = 0.0, absent_thresh = 0.025)

write_csv(as_tibble(out$X), file = "organized_data/cadotte_2013.csv")
write_csv(as_tibble(out$V), file = "organized_data/cadotte_2013_V.csv")
# save the tree
tree <- out$tree
save(tree, file =  "organized_data/cadotte_2013_tree.RData")
pdf("organized_data/cadotte_2013_tree.pdf")
plot(out$tree)
dev.off()
