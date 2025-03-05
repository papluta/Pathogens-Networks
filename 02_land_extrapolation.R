library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)
dens <- read.csv('Data/density.csv')

trans2021 <- read_csv('Data/Transects2021_20230713.csv') %>% filter(run == 2) %>% 
  rename(Site = Landscape_ID, Run = run, Bee_species = bee_species, Bee_group = bee_group) %>%
  mutate(Site = sub('Wm', 'WM', Site)) 
trans2022 <- read_csv('Data/Transects2022_20231110.csv') %>% filter(run == 2) %>% 
  rename(Site = Landscape_ID, Transect_type = Transect_tpe, Run = run) %>%
  mutate(Site = sub('Wm', 'WM', Site)) 
### LANDUSE

landuse2021 <- read.csv("Data/Landuse2021_500m_Combee_fieldedge_240930.csv")

land500.2021 <- landuse2021 %>%
  mutate(area_ha = area/10000) %>% 
  mutate(Landscape3 = sub('Wm', 'WM', Landscape3)) %>%
  group_by(Landscape3,Combi6) %>%
  summarise(sum=sum(area_ha)) %>%
  drop_na(Landscape3) %>%
  filter(Combi6 != '') %>%
  pivot_wider(names_from = Combi6, values_from = sum) %>%
  replace(is.na(.), 0) %>%
  rename(Site = Landscape3) %>%
  mutate(semi_natur = semi_natur + Grassy_str) %>%
  select(-Grassy_str) %>%
  rename(Grassy_str = Fieldedge) # changing name Fieldedge to grassy strip to match the transect file



source('01_landscape_processing.R')

land500.2022 <- radii[['500m']] %>% select(-`NA`)


### FLOWER COVER
hm <- trans2021 %>% distinct(Site, Transect_ID, Run, Transect_type)
fl.cv.2021 <- read_csv('Data/Flowercover2021_20230713.csv') %>% rename(Run = run, Site = Landscape_ID) %>%  
  mutate(Site = sub('Wm', 'WM', Site)) %>% 
  full_join(hm, by = c('Site', 'Run', 'Transect_ID')) %>% 
  filter(Run == 2) 
fl.cv.2022 <- read_csv('Data/Flowercover2022.csv') %>% filter(Run == 2)

distinct(fl.cv.2021, Transect_type) # the different habitat types
distinct(fl.cv.2022, Transect_type) # the different habitat types

# assigning the habitat types and AES labels
AES.2021 <- data.frame(Transect_type = distinct(fl.cv.2021, Transect_type),
                       Transect_type2 = c('Grassy_str', 'Fallow', 'Flower_fieBS12', 'Flower_fieBS11', 'CropBV1', 'Flower_fieBS2','semi_natur', 'Other_AUM'))

AES.2022 <- data.frame(Transect_type = distinct(fl.cv.2022, Transect_type),
                 Transect_type2 = c('Other_AUM', 'Flower_fieBS11', 'Flower_fieBS12', 'Flower_fieBS2', 'Grassy_str', 'CropBV1', 'Fallow', 'semi_natur'))


fl.cv2.2021 <- fl.cv.2021 %>% group_by(Run, Site, Transect_type) %>% summarise(flcv.m = mean(cover_perc)) %>% 
  left_join(AES.2021, by = 'Transect_type') %>% group_by(Site, Transect_type2) %>% summarise(flower_mean = mean(flcv.m))%>% 
  pivot_wider(names_from = Transect_type2, values_from = flower_mean) 

fl.cv2.2022 <- fl.cv.2022 %>% group_by(Run, Site, Transect_type) %>% summarise(flcv.m = mean(Total_flower_cover_percentage)) %>% 
  left_join(AES.2022, by = 'Transect_type') %>% group_by(Site, Transect_type2) %>% summarise(flower_mean = mean(flcv.m))%>% 
  pivot_wider(names_from = Transect_type2, values_from = flower_mean) 

fl.cv2.both <- rbind(fl.cv2.2021, fl.cv2.2022)

ggplot(fl.cv2.both, aes(Transect_type2, flcv.m, fill = as.factor(Transect_type2)))+
  geom_boxplot()+
  facet_wrap(~Year)+
  theme(legend.position = 'none')


## EXTRAPOLATION


extrapolation_function <- function(relative.data, landuse, normalization) {
  extrapolated <- relative.data %>% ungroup() %>%
    #taking mean across all sites in case of NAs (here even if the habitat was absent from a landscape, the mean will be calculated, but then it will be multiplied by zero anyway)
    mutate(Flower_fieBS2 = ifelse(is.na(Flower_fieBS2), mean(Flower_fieBS2, na.rm = T), Flower_fieBS2),
           Flower_fieBS12 = ifelse(is.na(Flower_fieBS12), mean(Flower_fieBS12, na.rm = T), Flower_fieBS12),
           Flower_fieBS11 = ifelse(is.na(Flower_fieBS11), mean(Flower_fieBS11, na.rm = T), Flower_fieBS11),
           Other_AUM = ifelse(is.na(Other_AUM), mean(Other_AUM, na.rm = T), Other_AUM),
           Fallow = ifelse(is.na(Fallow), mean(Fallow, na.rm = T), Fallow),
           CropBV1 = ifelse(is.na(CropBV1), mean(CropBV1, na.rm = T), CropBV1),
           semi_natur = ifelse(is.na(semi_natur), mean(semi_natur, na.rm = T), semi_natur), 
           Grassy_str = ifelse(is.na(Grassy_str), mean(Grassy_str, na.rm = T), Grassy_str)) %>%
    left_join(landuse, by = 'Site') %>% 
    #extrapolating (.x are the flower estimates, .y is the area in hectares)
    mutate(Org.ex = CropBV1.x * CropBV1.y/normalization,
           Fl_BS11_ex = Flower_fieBS11.x * Flower_fieBS11.y/normalization,
           Fl_BS12_ex = Flower_fieBS12.x * Flower_fieBS12.y/normalization,
           Fl_BS2_ex = Flower_fieBS2.x * Flower_fieBS2.y/normalization,
           Fallow_ex = Fallow.x * Fallow.y/normalization,
           Other_AUM_ex = Other_AUM.x * Other_AUM.y/normalization,
           Grassy_str_ex = Grassy_str.x * Grassy_str.y/normalization,
           semi_natur_ex = semi_natur.x * semi_natur.y/normalization) %>%
    mutate(sum = Org.ex + Fl_BS11_ex + Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex + Grassy_str_ex,
           SNH.ex = 0.5*Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex + Grassy_str_ex,
           Ann.fl.ex = Fl_BS11_ex + 0.5*Fl_BS12_ex ) %>%
    mutate(EX_per_agr = sum/(Crop + CropBV1.y + Flower_fieBS2.y + Flower_fieBS12.y + Flower_fieBS11.y + Fallow.y +
                                    Other_AUM.y + semi_natur.y + Grassy_str.y)) %>%
    select(Site, ends_with('.ex'), sum, EX_per_agr)
  return(extrapolated)
}

flower_cover2021 <- extrapolation_function(fl.cv2.2021, land500.2021, normalization = 100)
flower_cover2021$Year = 2021
flower_cover2022 <- extrapolation_function(fl.cv2.2022, land500.2022, normalization = 100)
flower_cover2022$Year = 2022


flower.both <- rbind(flower_cover2021, flower_cover2022) %>% rename(FL_per_agr = EX_per_agr)

flower.both %>%
  ggplot(aes(Site, FL_per_agr)) +
  geom_bar(stat = 'identity', col = 'black', fill = 'white')+
  facet_wrap(~Year, nrow = 2, ncol = 1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.98))

### hb abundance
ab.2021 <- trans2021 %>% filter(Bee_group == 'Honeybee' | Bee_group == 'Bumblebee' | Bee_group == 'Wildbee') %>% 
  filter(Bee_species != 'no_identificiation') %>% group_by(Run, Site, Bee_species, Transect_ID) %>% summarise(n = n()) 
ab.2021.list <- split(ab.2021, ab.2021$Bee_species)

ab.2021.list2 <- list()
abundance.2021.list <- list()
for (i in 1:length(ab.2021.list)) {
  ab.2021.list2[[i]] <- trans2021 %>% distinct(Site, Transect_ID, Run, Transect_type) %>% left_join(ab.2021.list[[i]], by = join_by('Site', 'Transect_ID', 'Run')) %>% 
    mutate(n = ifelse(is.na(n),0,n)) %>%
    group_by(Site, Transect_type) %>% summarise(Ab.m = mean(n)) %>%
    left_join(AES.2021, by = 'Transect_type') %>% mutate(Transect_type2 = ifelse(Transect_type == 'grassland', 'semi_natur', Transect_type2)) %>%
    select(-Transect_type) %>% group_by(Site, Transect_type2) %>% summarise(ab_mean = mean(Ab.m)) %>%
    pivot_wider(names_from = Transect_type2, values_from = ab_mean)
  abundance.2021.list[[i]] <- extrapolation_function(ab.2021.list2[[i]], land500.2021, normalization = 0.02)
}
names(abundance.2021.list) <- names(ab.2021.list)
abundance.2021 <- abundance.2021.list %>% bind_rows(.id = 'Species') %>% mutate(Species = sub(' agg.', '', Species)) %>% mutate(Year = 2021)

honeybees.2021 <- abundance.2021 %>% filter(Species == 'Apis mellifera')
lapidarius.2021 <- abundance.2021 %>% filter(Species == 'Bombus lapidarius') 


ab.2022 <- trans2022 %>% filter(Bee_group == 'Honeybee' | Bee_group == 'Bombus' | Bee_group == 'Wildbee') %>% 
  group_by(Run, Site, Bee_species, Transect_ID) %>% summarise(n = n()) %>% filter(Bee_species != 0 & !is.na(Bee_species))
ab.2022.list <- split(ab.2022, ab.2022$Bee_species)

ab.2022.list2 <- list()
abundance.2022.list <- list()
for (i in 1:length(ab.2022.list)) {
  ab.2022.list2[[i]] <- trans2022 %>% distinct(Site, Transect_ID, Run, Transect_type) %>% left_join(ab.2022.list[[i]], by = join_by('Site', 'Transect_ID', 'Run')) %>% 
  mutate(n = ifelse(is.na(n),0,n)) %>%
  group_by(Site, Transect_type) %>% summarise(Ab.m = mean(n)) %>%
  left_join(AES.2022, by = 'Transect_type') %>% mutate(Transect_type2 = ifelse(Transect_type == 'grassland', 'semi_natur', Transect_type2)) %>%
  select(-Transect_type) %>% group_by(Site, Transect_type2) %>% summarise(ab_mean = mean(Ab.m)) %>%
  pivot_wider(names_from = Transect_type2, values_from = ab_mean)
  abundance.2022.list[[i]] <- extrapolation_function(ab.2022.list2[[i]], land500.2022, normalization = 0.02)
}
names(abundance.2022.list) <- names(ab.2022.list)
abundance.2022 <- abundance.2022.list %>% bind_rows(.id = 'Species') %>% mutate(Species = sub(' agg.', '', Species)) %>% mutate(Year = 2022)

honeybees.2022 <- abundance.2022 %>% filter(Species == 'Apis mellifera') 
lapidarius.2022 <- abundance.2022 %>% filter(Species == 'Bombus lapidarius') 

honeybees.both <- rbind(honeybees.2021, honeybees.2022)
lapidarius.both <- rbind(lapidarius.2021, lapidarius.2022)
abundance.both <- rbind(abundance.2021, abundance.2022)

save(honeybees.both, lapidarius.both, abundance.both, flower.both, file = 'Data/land_extrapolated.RData')
