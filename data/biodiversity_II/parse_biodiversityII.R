library(tidyverse)
dt <- read_tsv("raw_data/e120_Plant aboveground biomass data.txt")
dt <- dt %>% select(Plot, Strip, Substrip, Year, Month, Species, `Biomass (g/m2)`) %>% filter(`Biomass (g/m2)` > 0, `Biomass (g/m2)` < 9999) %>% rowwise() %>% mutate(Date = paste0(Year, "-", Month)) %>% select(-Year, -Month)
# take average biomass per species by plot (i.e., average across Strips)
#dt_count_strips <- dt %>% select(Plot, Date, Strip) %>% distinct() %>% group_by(Plot, Date) %>% tally()
#dt <- dt %>% group_by(Plot, Date, Species) %>% summarise(Biomass = sum(`Biomass (g/m2)`), .groups = "drop") 
#dt <- dt %>% inner_join(dt_count_strips) %>% mutate(Biomass = Biomass / n) %>% select(-n)
dt <- dt %>% rowwise() %>% mutate(Plot = paste0(Plot, "_", Strip)) %>% select(-Strip, -Substrip) %>% 
  rename(Biomass = `Biomass (g/m2)`)
# save the species names for normalization
# dt %>% group_by(Species) %>% tally() %>% write_csv(file = "species_list.csv")
# normalize species names: all misc goes into "non_target"
dt <- dt %>% inner_join(read_csv("species_list.csv")) %>% select(-Species, -n, -species) %>% rename(species = target)
# merge indistinguishable 
dt <- dt %>% mutate(species = ifelse(species == "Solidago rigida", "Solidago nemoralis", species))
dt <- dt %>% mutate(species = ifelse(species == "Petalostemum candidum", "Petalostemum villosum", species))
dt <- dt %>% mutate(species = ifelse(species == "Petalostemum purpureum", "Petalostemum villosum", species))
dt <- dt %>% group_by(Date, Plot, species) %>% summarise(Biomass = sum(Biomass), .groups = "drop")
# only take plots/rings in which target biomass exceeds a threshold
dt <- dt %>% group_by(Date, Plot) %>% 
  mutate(total_biom = sum(Biomass)) %>% ungroup()

# spread species
dt <- dt %>% pivot_wider(names_from = species, values_from = Biomass, values_fill = 0) %>% mutate(prop_nontarget = non_target / total_biom)
THRESH_NONTARGET <- 0.1
# The distribution of proportion of biomass being off-target is peaked at 1;
dt <- dt %>% filter(prop_nontarget < THRESH_NONTARGET)
# select samplingdate and treatments such that at least 90 plots are present
to_select <- dt %>% group_by(Date) %>% tally() %>% arrange(desc(n)) %>% filter(n > 90)

system("rm selected_data/*")

for (i in 1:nrow(to_select)){
  mydate <- to_select$Date[i]
  label_file <- paste0("selected_data/bioII_", mydate, ".csv")
  dt2 <- dt %>% filter(Date == mydate) %>% 
    select(-Date, -Plot, -total_biom, -prop_nontarget, -non_target)
  write_csv(dt2, file = label_file)
}
