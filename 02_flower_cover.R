library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)

fl.cv <- read_csv('Data/Flowercover2022.csv')


distinct(fl.cv, Transect_type) # the different habitat types

# assigning the habitat types AES labels
hm <- data.frame(Transect_type = distinct(fl.cv, Transect_type),
                 Transect_type2 = c('Other_AUM', 'Flower_fieBS11', 'Flower_fieBS12', 'Flower_fieBS2', 'Fallow', 'Grassland','Crop', 'CropBV1', 'GrasslandBV1'),
                 Edge_only = c(0, 0, 0, 0, 0, 0, 1, 0, 0)) # when calculating area treat grassy strips as an edge
hm2 <- data.frame(Transect_type = distinct(fl.cv, Transect_type),
                  AES = c('SNH', 'Flower', 'Flower', 'Flower','SNH','SNH','SNH','Crop','Org.farm'))
fl.cv2 <- fl.cv %>% mutate(Date = as.Date(Date, '%m/%d/%Y')) %>% group_by(Run, Site, Transect_type) %>% summarise(flcv.m = mean(Total_flower_cover_percentage), Date = max(Date)) %>% 
  left_join(hm, by = 'Transect_type') %>% left_join(hm2, by = 'Transect_type') %>% 
  mutate(Month = ifelse(Run ==1, 'June', ifelse(Run==2, 'July', 'August')))


ggplot(fl.cv2, aes(Transect_type2, flcv.m, fill = as.factor(Run)))+
  geom_boxplot()


load('Data/radii.RData')
load('Data/landscape_metrics2022.RData')

## CALCULATING TOTAL FLOWER COVER AREA (EXTRAPOLATION)

# 500 m
land500 <- radii[['500m']]
fl.area500 <- land500 %>% pivot_longer(cols = Crop:Plantation,names_to = 'Transect_type2', values_to = 'area') %>% filter(Transect_type2 %in% unique(fl.cv2$Transect_type2)) %>%
  right_join(fl.cv2, by = c('Site', 'Transect_type2')) %>% 
  # extrapolating the flower cover to the habitat proportion of each habitat type except for grassy strips in 2000 m radius
  mutate(flcv.area = ifelse(Edge_only == 0, ((flcv.m/100) * area)/0.2, NA)) %>%
  left_join(Het.lm500 %>% dplyr::select(Site, starts_with('Edge')) %>% dplyr::select(!c(Edge.dens)) %>%
              pivot_longer(cols = starts_with('Edge'), names_to = 'Transect_type2', values_to = 'Edge') %>%
              mutate(Transect_type2 = gsub('Edge.', '', Transect_type2)), by = c('Site', 'Transect_type2')) %>%
  # extrapolating the flower cover to length of the crop field edge (grassy strips) in 2000 m radius
  mutate(flcv.edge = ifelse(Edge_only == 1, ((flcv.m/100) * (Edge*0.0001))/0.2, NA)) %>% 
  mutate(flcv.edge = ifelse(is.na(Edge), NA, flcv.edge)) %>% mutate(flcv.area = ifelse(is.na(area), NA, flcv.area)) %>%
  # summing up the flower cover of the area and the edge
  mutate(flower.cover = ifelse(is.na(flcv.area), flcv.edge, flcv.area)) %>% drop_na(flower.cover)

land1000 <- radii[['1000m']] ## in 1000 m radius there are additional habitats: Allotments and BW1 so the pivoting changes a bit
fl.area1000 <- land1000 %>% pivot_longer(cols = Allotments:Plantation,names_to = 'Transect_type2', values_to = 'area') %>% filter(Transect_type2 %in% unique(fl.cv2$Transect_type2)) %>%
  right_join(fl.cv2, by = c('Site', 'Transect_type2')) %>% mutate(flcv.area = ifelse(Edge_only == 0, ((flcv.m/100) * area)/0.2, NA)) %>%
  left_join(Het.lm1000 %>% dplyr::select(Site, starts_with('Edge')) %>%
              pivot_longer(cols = starts_with('Edge'), names_to = 'Transect_type2', values_to = 'Edge') %>%
              mutate(Transect_type2 = gsub('Edge.', '', Transect_type2)), by = c('Site', 'Transect_type2')) %>%
  mutate(flcv.edge = ifelse(Edge_only == 1, ((flcv.m/100) * (Edge*0.0001))/0.2, NA)) %>% 
  mutate(flcv.edge = ifelse(is.na(Edge), NA, flcv.edge)) %>% mutate(flcv.area = ifelse(is.na(area), NA, flcv.area)) %>%
  mutate(flower.cover = ifelse(is.na(flcv.area), flcv.edge, flcv.area)) %>% drop_na(flower.cover)

land1500 <- radii[['1500m']] ## before  pivot_longer(cols = Crop:Plantation,
fl.area1500 <- land1500 %>% pivot_longer(cols = Allotments:Plantation,names_to = 'Transect_type2', values_to = 'area') %>% filter(Transect_type2 %in% unique(fl.cv2$Transect_type2)) %>%
  right_join(fl.cv2, by = c('Site', 'Transect_type2')) %>% mutate(flcv.area = ifelse(Edge_only == 0, ((flcv.m/100) * area)/0.2, NA)) %>%
  left_join(Het.lm1500 %>% dplyr::select(Site, starts_with('Edge')) %>%
              pivot_longer(cols = starts_with('Edge'), names_to = 'Transect_type2', values_to = 'Edge') %>%
              mutate(Transect_type2 = gsub('Edge.', '', Transect_type2)), by = c('Site', 'Transect_type2')) %>%
  mutate(flcv.edge = ifelse(Edge_only == 1, ((flcv.m/100) * (Edge*0.0001))/0.2, NA)) %>% 
  mutate(flcv.edge = ifelse(is.na(Edge), NA, flcv.edge)) %>% mutate(flcv.area = ifelse(is.na(area), NA, flcv.area)) %>%
  mutate(flower.cover = ifelse(is.na(flcv.area), flcv.edge, flcv.area)) %>% drop_na(flower.cover)

land2000 <- radii[['2000m']]
fl.area2000 <- land2000 %>% pivot_longer(cols = Allotments:Plantation,names_to = 'Transect_type2', values_to = 'area') %>% filter(Transect_type2 %in% unique(fl.cv2$Transect_type2)) %>%
  right_join(fl.cv2, by = c('Site', 'Transect_type2')) %>% mutate(flcv.area = ifelse(Edge_only == 0, ((flcv.m/100) * area)/0.2, NA)) %>%
  left_join(Het.lm2000 %>% dplyr::select(Site, starts_with('Edge')) %>% dplyr::select(!c(Edge.dens)) %>%
              pivot_longer(cols = starts_with('Edge'), names_to = 'Transect_type2', values_to = 'Edge') %>%
              mutate(Transect_type2 = gsub('Edge.', '', Transect_type2)), by = c('Site', 'Transect_type2')) %>%
  mutate(flcv.edge = ifelse(Edge_only == 1, ((flcv.m/100) * (Edge*0.0001))/0.2, NA)) %>% 
  mutate(flcv.edge = ifelse(is.na(Edge), NA, flcv.edge)) %>% mutate(flcv.area = ifelse(is.na(area), NA, flcv.area))%>%
  mutate(flower.cover = ifelse(is.na(flcv.area), flcv.edge, flcv.area)) %>% drop_na(flower.cover)

save(fl.area500, fl.area1000, fl.area1500, fl.area2000, file = 'Data/fl_cover_extr.RData')

dens <- read_csv('Data/density.csv')

hb.ab <- read_csv('Data/hb_abundance2022.csv') %>% rename(Date = date, Site = Landscape_ID, Transect_type = Transect_tpe) %>% mutate(Date = as.Date(Date, '%m/%d/%Y')) %>% group_by(Run, Site, Transect_type) %>% summarise(hb.m = mean(Honeybee), Date = max(Date)) %>% 
  left_join(hm, by = 'Transect_type') %>% left_join(hm2, by = 'Transect_type') %>% mutate(Month = ifelse(Run ==1, 'June', ifelse(Run==2, 'July', 'August'))) %>%
  mutate(hb.m = round(hb.m, digits = 2)) %>% left_join(dens, by = 'Site')

ggplot(hb.ab, aes(as.factor(Run), log10(hb.m+1), fill = as.factor(Density)))+
  geom_boxplot()

ggplot(hb.ab %>% filter(Run != 1), aes(as.factor(Density), log10(hb.m+1), fill = as.factor(Density)))+
  geom_boxplot()

ggplot(hb.ab %>% filter(Run != 1), aes(Transect_type, log10(hb.m+1), fill = as.factor(Density)))+
  geom_boxplot()
