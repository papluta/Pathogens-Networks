library(readr)
library(tidyverse)
source('01_custom_functions.R')

dens <- read.csv('Data/density.csv')

trans2021 <- read_csv('Data/Transects2021.csv')
trans2022 <- read_csv('Data/Transects2022.csv') 

### LANDUSE

land500.2021 <- read.csv("Data/Landuse2021_500m.csv")

land500.2022 <- read.csv('Data/Landuse2022_500m.csv')


### FLOWER COVER
hm <- trans2021 %>% distinct(Site, Transect_ID, Run, Transect_type)
fl.cv.2021 <- read_csv('Data/Flowercover2021.csv')
fl.cv.2022 <- read_csv('Data/Flowercover2022.csv') 


distinct(fl.cv.2021, Transect_type) # the different habitat types
distinct(fl.cv.2022, Transect_type) # the different habitat types

## Calculating flower cover

fl.cv2.2021 <- fl.cv.2021 %>% group_by(Run, Site, Transect_type) %>% summarise(flcv.m = mean(cover_perc)) %>% 
  mutate(Transect_type2 = case_when(
    Transect_type == "grassy_stripe" ~ "Grassy_str",
    Transect_type == "fallow" ~ "Fallow",
    Transect_type == "flower_fieBS12" ~ "Flower_fieBS12",
    Transect_type == "flower_fieBS11" ~ "Flower_fieBS11",
    Transect_type == "flower_fieBS2" ~ "Flower_fieBS2",
    Transect_type == "organic" ~ "CropBV1",
    Transect_type == "grassland" ~ "semi_natur",
    Transect_type == "aum" ~ "Other_AUM",
    TRUE ~ NA)) %>% 
  group_by(Site, Transect_type2) %>% 
  summarise(flower_mean = mean(flcv.m)) %>% 
  pivot_wider(names_from = Transect_type2, values_from = flower_mean) 

fl.cv2.2022 <- fl.cv.2022 %>% 
  group_by(Run, Site, Transect_type) %>% 
  summarise(flcv.m = mean(Total_flower_cover_percentage)) %>% 
  mutate(Transect_type2 = case_when(
    Transect_type == "gs" ~ "Grassy_str",
    Transect_type == "fallow" ~ "Fallow",
    Transect_type == "bs12" ~ "Flower_fieBS12",
    Transect_type == "bs11" ~ "Flower_fieBS11",
    Transect_type == "bs2" ~ "Flower_fieBS2",
    Transect_type == "organic" ~ "CropBV1",
    Transect_type == "grassland/gs" ~ "semi_natur",
    Transect_type == "grassland" ~ "semi_natur",
    Transect_type == "aum" ~ "Other_AUM",
    TRUE ~ NA)) %>%   
  group_by(Site, Transect_type2) %>% 
  summarise(flower_mean = mean(flcv.m))%>% 
  pivot_wider(names_from = Transect_type2, values_from = flower_mean) 

fl.cv2.both <- rbind(fl.cv2.2021, fl.cv2.2022)



## EXTRAPOLATION

# normalisation 100 beacause flower cover is a proportion
flower_cover2021 <- extrapolation_function(fl.cv2.2021, land500.2021, normalization = 100)
flower_cover2021$Year = 2021
flower_cover2022 <- extrapolation_function(fl.cv2.2022, land500.2022, normalization = 100)
flower_cover2022$Year = 2022


flower.both <- rbind(flower_cover2021, flower_cover2022) %>% rename(flower_dens = sum)

flower.both %>%
  ggplot(aes(Site, flower_dens)) +
  geom_bar(stat = 'identity', col = 'black', fill = 'white')+
  facet_wrap(~Year, nrow = 2, ncol = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.98))

### Bee abundance

ab.2021 <- trans2021 %>% 
  # excluding non-bees
  filter(Bee_group == 'Honeybee' | Bee_group == 'Bumblebee' | Bee_group == 'Wildbee') %>% 
  filter(Bee_species != 'no_identificiation') %>% 
  group_by(Run, Site, Bee_species, Transect_ID) %>% 
  summarise(n = n()) 

ab.2021.list <- split(ab.2021, ab.2021$Bee_species)

ab.2021.list2 <- list()
abundance.2021.list <- list()
# calculate bee density per bee species, considering 0 when a species is missing from a site and then extrapolate
for (i in 1:length(ab.2021.list)) {
  ab.2021.list2[[i]] <- trans2021 %>% distinct(Site, Transect_ID, Run, Transect_type) %>% 
    left_join(ab.2021.list[[i]], by = join_by('Site', 'Transect_ID', 'Run')) %>% 
    mutate(n = ifelse(is.na(n),0,n)) %>%
    group_by(Site, Transect_type) %>% 
    summarise(Ab.m = mean(n)) %>%
    mutate(Transect_type2 = case_when(
      # same as flower cover
      Transect_type == "grassy_stripe" ~ "Grassy_str",
      Transect_type == "fallow" ~ "Fallow",
      Transect_type == "flower_fieBS12" ~ "Flower_fieBS12",
      Transect_type == "flower_fieBS11" ~ "Flower_fieBS11",
      Transect_type == "flower_fieBS2" ~ "Flower_fieBS2",
      Transect_type == "organic" ~ "CropBV1",
      Transect_type == "grassland" ~ "semi_natur",
      Transect_type == "aum" ~ "Other_AUM",
      TRUE ~ NA)) %>% 
    select(-Transect_type) %>% 
    group_by(Site, Transect_type2) %>% 
    summarise(ab_mean = mean(Ab.m)) %>%
    pivot_wider(names_from = Transect_type2, values_from = ab_mean)
  # normalisation 0.02 because transect walks were 200m2
  abundance.2021.list[[i]] <- suppressWarnings(extrapolation_function(ab.2021.list2[[i]], land500.2021, normalization = 0.02))
}

names(abundance.2021.list) <- names(ab.2021.list)
abundance.2021 <- abundance.2021.list %>% 
  bind_rows(.id = 'Species') %>% 
  mutate(Species = sub(' agg.', '', Species)) %>% 
  mutate(Year = 2021) %>%
  rename(bee_abundance = sum)

honeybees.2021 <- abundance.2021 %>% filter(Species == 'Apis mellifera')
lapidarius.2021 <- abundance.2021 %>% filter(Species == 'Bombus lapidarius') 


ab.2022 <- trans2022 %>% 
  filter(Bee_group == 'Honeybee' | Bee_group == 'Bombus' | Bee_group == 'Wildbee') %>% 
  group_by(Run, Site, Bee_species, Transect_ID) %>% 
  summarise(n = n()) %>% 
  filter(Bee_species != 0 & !is.na(Bee_species))
ab.2022.list <- split(ab.2022, ab.2022$Bee_species)

ab.2022.list2 <- list()
abundance.2022.list <- list()
for (i in 1:length(ab.2022.list)) {
  ab.2022.list2[[i]] <- trans2022 %>% 
    distinct(Site, Transect_ID, Run, Transect_type) %>% 
    left_join(ab.2022.list[[i]], by = join_by('Site', 'Transect_ID', 'Run')) %>% 
    mutate(n = ifelse(is.na(n),0,n)) %>%
    group_by(Site, Transect_type) %>% 
    summarise(Ab.m = mean(n)) %>%
    mutate(Transect_type2 = case_when(
      Transect_type == "gs" ~ "Grassy_str",
      Transect_type == "fallow" ~ "Fallow",
      Transect_type == "bs12" ~ "Flower_fieBS12",
      Transect_type == "bs11" ~ "Flower_fieBS11",
      Transect_type == "bs2" ~ "Flower_fieBS2",
      Transect_type == "organic" ~ "CropBV1",
      Transect_type == "grassland/gs" ~ "semi_natur",
      Transect_type == "grassland" ~ "semi_natur",
      Transect_type == "aum" ~ "Other_AUM",
      TRUE ~ NA)) %>%     
    select(-Transect_type) %>% 
    group_by(Site, Transect_type2) %>% 
    summarise(ab_mean = mean(Ab.m)) %>%
    pivot_wider(names_from = Transect_type2, values_from = ab_mean)
    abundance.2022.list[[i]] <- suppressWarnings(extrapolation_function(ab.2022.list2[[i]], land500.2022, normalization = 0.02))
}
names(abundance.2022.list) <- names(ab.2022.list)
abundance.2022 <- abundance.2022.list %>% 
  bind_rows(.id = 'Species') %>% 
  mutate(Species = sub(' agg.', '', Species)) %>% 
  mutate(Year = 2022) %>%
  rename(bee_abundance = sum)

honeybees.2022 <- abundance.2022 %>% filter(Species == 'Apis mellifera') 
lapidarius.2022 <- abundance.2022 %>% filter(Species == 'Bombus lapidarius') 

honeybees.both <- rbind(honeybees.2021, honeybees.2022)
lapidarius.both <- rbind(lapidarius.2021, lapidarius.2022)
abundance.both <- rbind(abundance.2021, abundance.2022)

