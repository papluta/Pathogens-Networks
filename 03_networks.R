source('02_land_extrapolation.R')

library(bipartite)
library(purrr)
library(spaa)
library(ggpubr)

trans2021.2 <- trans2021 %>% filter(Bee_group == 'Honeybee' | Bee_group == 'Bumblebee' | Bee_group == 'Wildbee') %>% 
  filter(Bee_species != 'no_identificiation')
trans2022.2 <- trans2022 %>% filter(Bee_group == 'Honeybee' | Bee_group == 'Bombus' | Bee_group == 'Wildbee') %>% 
  filter(Bee_species != 'no_identificiation')

levels(as.factor(trans2021.2$Bee_group))
levels(as.factor(trans2021.2$Bee_species))
levels(as.factor(trans2022.2$Bee_group))
levels(as.factor(trans2022.2$Bee_species))


adj.matrix.2021 <- frame2webs(trans2021.2, varnames = c("visited_flower_species", "Bee_species",
                                                        "Site"), emptylist = TRUE)
adj.matrix.2022 <- frame2webs(trans2022.2, varnames = c("visited_flower_species", "Bee_species",
                                            "Site"), emptylist = TRUE)


net.met2021 <- adj.matrix.2021 %>% lapply(networklevel, index=c("connectance", "NODF"), level = 'higher', 
                                            weighted = TRUE)
net.met2022 <- adj.matrix.2022 %>% lapply(networklevel, index=c("connectance", "NODF"), level = 'higher', 
                                          weighted = TRUE)

net.met.long2021 <- net.met2021 %>% bind_rows(.id = 'Site') %>% mutate(Year = 2021)
net.met.long2022 <- net.met2022 %>% bind_rows(.id = 'Site') %>% mutate(Year = 2022)
net.met.long.both <- rbind(net.met.long2021, net.met.long2022)

net.size2021 <- map(adj.matrix.2021, ~ data.frame(bee_rich = nrow(.x), plant_rich = ncol(.x))) %>% 
  bind_rows(.id = "Site") %>% mutate(Year = 2021)
net.size2022 <- map(adj.matrix.2022, ~ data.frame(bee_rich = nrow(.x), plant_rich = ncol(.x)))%>% 
  bind_rows(.id = "Site") %>% mutate(Year = 2022)

net.size <- rbind(net.size2021, net.size2022) %>% 
  rowwise() %>%
  mutate(sum_nodes = sum(bee_rich, plant_rich))

network.metrics.both <- data.frame(Site = net.met.long.both$Site, 
                                   Connectance = net.met.long.both$connectance,
                                   NODF = net.met.long.both$NODF,
                                   Year = net.met.long.both$Year) %>% 
  left_join(net.size, by = c("Site", "Year"))

#### CALCULATING INDIVIDUAL NICHE OVERLAP

morisita.2021 <- map(adj.matrix.2021, function(x) as.matrix(niche.overlap(x, method = 'morisita')))
morisita.2022 <- map(adj.matrix.2022, function(x) as.matrix(niche.overlap(x, method = 'morisita')))

# raw morisita

morisita2021.hbbl <- map(morisita.2021, ~ as.data.frame(.x) %>% 
                           select(any_of(c("Apis mellifera", "Bombus lapidarius", "Lasioglossum pauxillum", "Andrena minutula"))) %>%
                                    rownames_to_column(var = "Species"))
morisita2022.hbbl <- map(morisita.2022, ~ as.data.frame(.x) %>% 
                           select(any_of(c("Apis mellifera", "Bombus lapidarius", "Lasioglossum pauxillum", "Andrena minutula"))) %>%
                                           rownames_to_column(var = "Species"))

morisita2021.hbbl.2 <- morisita2021.hbbl %>% bind_rows(.id = "Site") %>% mutate(`Andrena minutula` = NA, Year = 2021)
morisita2022.hbbl.2 <- morisita2022.hbbl %>% bind_rows(.id = "Site") %>% mutate(Year = 2022)

morisita.raw <- rbind(morisita2021.hbbl.2, morisita2022.hbbl.2)



