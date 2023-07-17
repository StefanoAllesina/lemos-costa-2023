source("../general_code/process_dataset.R")
library(readxl)

dt <- as.data.frame(readxl::read_xlsx("vanruijven-specbiomass-2000-2008.xlsx"))
dt[dt == "x"] <- 0
dt <- dt %>% select(-Richness, -total) %>%
  select(-Plotnr, -Block) %>% 
  mutate(plot = row_number()) %>% mutate_all(.funs = as.numeric) %>% 
  pivot_longer(names_to = "species", values_to = "biomass", cols = -c("plot", "year")) %>% 
  filter(biomass > 0)

# change species names to normalize them
dt <- dt %>% 
  rename(code = species) %>% 
  inner_join(read_csv("species_codes.csv") %>% 
               select(code, species), by = "code") %>% 
  select(-code) 


year_s <- dt %>% select(year) %>% distinct()
spinfo <- read_csv("species_codes.csv") %>% select(-code)


for (myyear in 1:nrow(year_s)){
  dt2 <- dt %>% filter(year == year_s$year[myyear]) %>% 
    select(plot, species, biomass) 
  print(year_s$year[myyear])
  out <- get_data(dt2, spinfo, min_comm = 6, biom_thresh = 1.0, absent_thresh = 0.025)
  write_csv(as_tibble(out$X), file = paste0("organized_data/vr_", myyear, ".csv"))
  write_csv(as_tibble(out$V), file = paste0("organized_data/vr_", myyear, "_V.csv"))
  tree <- out$tree
  save(tree, file =  paste0("organized_data/vr_", myyear, "_tree.RData"))
  pdf(paste0("organized_data/vr_", myyear, "_tree.pdf"))
  plot(out$tree)
  dev.off()
}
