library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)

fl.cv <- read_csv('Data/Flowercover2022.csv')


distinct(fl.cv, Transect_type) # the different habitat types

# assigning the habitat types AES labels
hm <- data.frame(Transect_type = distinct(fl.cv, Transect_type),
                 Transect_type2 = c('Other_AUM', 'Flower_fieBS11', 'Flower_fieBS12', 'Flower_fieBS2', 'Fallow', 'Grassland','Grassy_str', 'CropBV1', 'GrasslandBV1'),
                 Edge_only = c(0, 0, 0, 0, 0, 0, 1, 0, 0)) # when calculating area treat grassy strips as an edge
hm2 <- data.frame(Transect_type = distinct(fl.cv, Transect_type),
                  AES = c('SNH', 'Flower', 'Flower', 'Flower','SNH','SNH','SNH','Org.farm','SNH'))
fl.cv2 <- fl.cv %>% mutate(Date = as.Date(Date, '%m/%d/%Y')) %>% group_by(Run, Site, Transect_type) %>% summarise(flcv.m = mean(Total_flower_cover_percentage), Date = max(Date)) %>% 
  left_join(hm, by = 'Transect_type') %>% left_join(hm2, by = 'Transect_type') %>% 
  mutate(Month = ifelse(Run ==1, 'June', ifelse(Run==2, 'July', 'August')))


ggplot(fl.cv2, aes(Transect_type2, flcv.m, fill = as.factor(Run)))+
  geom_boxplot()

ggplot(fl.cv2 %>% filter(Run == 2), aes(Transect_type2, flcv.m))+
  geom_boxplot()


load('Data/240617_landuse2022.RData')
load('Data/240617landscape_metrics2022.RData')

## CALCULATING TOTAL FLOWER COVER AREA (EXTRAPOLATION)
#flower.AES <- fl.cv2 %>% filter(Run == 2) %>% group_by(AES, Transect_type2) %>% summarise(mean.flower = mean(flcv.m)) %>% group_by(AES) %>% summarise(mean.flower = mean(mean.flower)) %>% 
#  pivot_wider(names_from = AES, values_from = mean.flower)
flower_cover <- fl.cv2 %>% filter(Run == 2) %>% group_by(Site, Transect_type2, AES) %>% summarise(flower_mean = mean(flcv.m)) %>% 
  group_by(Site, AES) %>% summarise(flower_mean = mean(flower_mean)) %>% pivot_wider(names_from = AES, values_from = flower_mean) %>% ungroup() %>%
  mutate(Flower = ifelse(is.na(Flower), mean(Flower, na.rm = T), Flower),
         Org.farm = ifelse(is.na(Org.farm), mean(Org.farm, na.rm = T), Org.farm),
         SNH = ifelse(is.na(SNH), mean(SNH, na.rm = T), SNH))

mean(flower_cover$Flower, na.rm = T)


## CANNOT CALCULATE IT AT THE TRANSECT TYPE LEVEL BECAUSE GRASSY STRIPS ARE NOT MAPPED OUT COMPLETELY ##

#flower_cover <- fl.cv2 %>% filter(Run == 2) %>% group_by(Site, Transect_type2) %>% summarise(flower_mean = mean(flcv.m))%>% 
#  pivot_wider(names_from = Transect_type2, values_from = flower_mean) %>% 
#  left_join(radii[['1000m']], by = 'Site') %>%
#  mutate(Org_ex = CropBV1.x * CropBV1.y,
#         Fl_BS11_ex = Flower_fieBS11.x * Flower_fieBS11.y,
#         Fl_BS12_ex = Flower_fieBS12.x * Flower_fieBS12.y,
#         Fl_BS2_ex = Flower_fieBS2.x * Flower_fieBS2.y,
#         Fallow_ex = Fallow.x * Fallow.y,
#         Other_AUM_ex = Other_AUM.x * Other_AUM.y,
#         Grassy
#  )


# 500 m
land500 <- subset(land.all, radius == '500m')

flower500 <- land500 %>% select(Site, Org.farm, Total.fl, SNH.nf) %>% 
  rename(Org.farm.h = Org.farm) %>%
  cbind(flower.AES) %>%
  mutate(Org.fl = Org.farm.h * Org.farm,
         Flower.fl = Total.fl * Flower,
         SNH.fl = SNH.nf * SNH)%>%
  select(Site, Org.fl, Flower.fl, SNH.fl)

flower500$sum.fl <- rowSums(flower500[2:4], na.rm = T)


land1000 <- subset(land.all, radius == '1000m')

flower1000 <- land1000 %>% select(Site, Total.fl, Org.farm, SNH.nf, Crop) %>% 
  rename(Org.farm.h = Org.farm, Total.fl.h = Total.fl, SNH.nf.h = SNH.nf) %>% left_join(flower_cover, by = 'Site') %>%
  mutate(Org.fl = Org.farm.h * Org.farm/100,
         Flower.fl = Total.fl.h * Flower/100,
         SNH.fl = SNH.nf.h * SNH/100)
flower1000$sum.fl <- rowSums(flower100_AES[9:11], na.rm = T)
flower1000$FL_per_agr <- with(flower100_AES, sum.fl/(Org.farm.h + SNH.nf.h + Total.fl.h + Crop))


### hb abundance
hb.ab <- read_csv('Data/hb_abundance2022.csv') %>% rename(Date = date, Site = Landscape_ID, Transect_type = Transect_tpe) %>% mutate(Date = as.Date(Date, '%m/%d/%Y')) %>% group_by(Run, Site, Transect_type) %>% summarise(hb.m = mean(Honeybee), Date = max(Date)) %>% 
  left_join(hm, by = 'Transect_type') %>% left_join(hm2, by = 'Transect_type') %>% mutate(Month = ifelse(Run ==1, 'June', ifelse(Run==2, 'July', 'August'))) %>% mutate(AES = ifelse(Transect_type == 'grassland', 'SNH', AES), Transect_type2 = ifelse(Transect_type == 'grassland', 'grassland', Transect_type2))
hb.ab2 <- hb.ab %>% filter(Run == 2) %>% group_by(Transect_type2, Site, AES) %>% summarise(hb.m = mean(hb.m)) %>% group_by(Site, AES) %>%
  summarise(hb.m = mean(hb.m)) %>% pivot_wider(names_from = AES, values_from = hb.m) %>% ungroup() %>%
  mutate(Flower = ifelse(is.na(Flower), mean(Flower, na.rm = T), Flower),
         Org.farm = ifelse(is.na(Org.farm), mean(Org.farm, na.rm = T), Org.farm),
         SNH = ifelse(is.na(SNH), mean(SNH, na.rm = T), SNH))


honeybees500 <- land500 %>% select(Site, Org.farm, Total.fl, SNH.nf) %>% 
  rename(Org.farm.h = Org.farm) %>%
  cbind(hb.ab2) %>%
  mutate(Org.hb = Org.farm.h * Org.farm / 0.02,
         Flower.hb = Total.fl * Flower / 0.02,
         SNH.hb = SNH.nf * SNH / 0.02) %>%
  select(Site, Org.hb, Flower.hb, SNH.hb)

honeybees500$sum.hb <- rowSums(honeybees500[2:4], na.rm = T)

honeybees1000 <- land1000 %>% select(Site, Org.farm, Total.fl, SNH.nf, Crop) %>% 
  rename(Org.farm.h = Org.farm) %>%
  left_join(hb.ab2, by = 'Site') %>%
  mutate(Org.hb = Org.farm.h * Org.farm / 0.02,
         Flower.hb = Total.fl * Flower / 0.02,
         SNH.hb = SNH.nf * SNH / 0.02) 
honeybees1000$sum.hb <- rowSums(honeybees1000[9:11], na.rm = T)
honeybees1000$HB_per_agr <- with(honeybees1000, sum.hb/(Org.farm.h + SNH.nf + Total.fl + Crop))
