library(tidyverse)
source("../general_code/process_dataset.R")

#system("rm selected_data/bioII_2006-7.csv") # causes a parsing error
#system("rm selected_data/bioII_2007-8.csv")
system("rm organized_data/*")

# read the data
for (ffn in list.files("selected_data/")){
ffn <- paste0("selected_data/", ffn)
dt <- read_csv(ffn)

dt <- dt %>% mutate(measure = row_number()) %>% 
  pivot_longer(names_to = "code", values_to = "biomass", cols = -measure) %>% filter(biomass > 0)

# change species names to normalize them
dt <- dt %>% 
  inner_join(read_csv("Spp_codes.csv") %>% 
               select(code, species), by = "code") %>% 
  rename(plot = measure) %>% 
  select(-code) 

spinfo <- read_csv("Spp_codes.csv") %>% select(-code)

out <- get_data(dt, spinfo, min_comm = 6, biom_thresh = 0.25, min_plots = 10, absent_thresh = 0.025)


  if (!is.null(out)){
    labelfile <- tools::file_path_sans_ext(basename(ffn))
    write_csv(as_tibble(out$X), file = paste0("organized_data/",labelfile, ".csv"))
    write_csv(as_tibble(out$V), file = paste0("organized_data/",labelfile, "_V.csv"))
    tree <- out$tree
    save(tree, file =  paste0("organized_data/",labelfile, "_tree.RData"))
    pdf(paste0("organized_data/",labelfile, "_tree.pdf"))
    plot(out$tree)
    dev.off()
  }
}