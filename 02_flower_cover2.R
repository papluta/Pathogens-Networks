library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)

fl.cv <- read_csv('Data/Flowercover2022.csv')


distinct(fl.cv, Transect_type) # the different habitat types

# assigning the habitat types and AES labels
hm <- data.frame(Transect_type = distinct(fl.cv, Transect_type),
                 Transect_type2 = c('Other_AUM', 'Flower_fieBS11', 'Flower_fieBS12', 'Flower_fieBS2', 'Fallow', 'semi_natur','semi_natur', 'CropBV1', 'semi_natur'))

## not using this categorization in the end:
hm2 <- data.frame(Transect_type = distinct(fl.cv, Transect_type),
                  AES = c('SNH', 'Flower', 'Flower', 'Flower','SNH','SNH','SNH','Org.farm','SNH'))

fl.cv2 <- fl.cv %>% mutate(Date = as.Date(Date, '%m/%d/%Y')) %>% group_by(Run, Site, Transect_type) %>% summarise(flcv.m = mean(Total_flower_cover_percentage), Date = max(Date)) %>% 
  left_join(hm, by = 'Transect_type') %>% left_join(hm2, by = 'Transect_type') %>% 
  mutate(Month = ifelse(Run ==1, 'June', ifelse(Run==2, 'July', 'August')))


ggplot(fl.cv2, aes(Transect_type2, flcv.m, fill = as.factor(Run)))+
  geom_boxplot()

ggplot(fl.cv2 %>% filter(Run == 2), aes(Transect_type2, flcv.m))+
  geom_boxplot()


land1000 <- read.csv('Data/Landuse_Combee2022_1km.csv')
land1000_2 <- land1000 %>%
  mutate(area_ha = area1/10000) %>% 
  group_by(Landscape3,Combi6) %>%
  summarise(sum=sum(area_ha)) %>%
  drop_na(Landscape3) %>%
  pivot_wider(names_from = Combi6, values_from = sum) %>%
  replace(is.na(.), 0)%>%
  rename(Site = Landscape3) %>%
  mutate(semi_natur = semi_natur + Grassy_str)

## CALCULATING TOTAL FLOWER COVER AREA (EXTRAPOLATION)

flower_cover1000 <- fl.cv2 %>% filter(Run == 2) %>% group_by(Site, Transect_type2) %>% summarise(flower_mean = mean(flcv.m))%>% 
  pivot_wider(names_from = Transect_type2, values_from = flower_mean) %>% ungroup() %>%
  #taking mean across all sites in case of NAs
  mutate(Flower_fieBS2 = ifelse(is.na(Flower_fieBS2), mean(Flower_fieBS2, na.rm = T), Flower_fieBS2),
         Flower_fieBS12 = ifelse(is.na(Flower_fieBS12), mean(Flower_fieBS12, na.rm = T), Flower_fieBS12),
         Flower_fieBS11 = ifelse(is.na(Flower_fieBS11), mean(Flower_fieBS11, na.rm = T), Flower_fieBS11),
         Other_AUM = ifelse(is.na(Other_AUM), mean(Other_AUM, na.rm = T), Other_AUM),
         Fallow = ifelse(is.na(Fallow), mean(Fallow, na.rm = T), Fallow),
         CropBV1 = ifelse(is.na(CropBV1), mean(CropBV1, na.rm = T), CropBV1),
         semi_natur = ifelse(is.na(semi_natur), mean(semi_natur, na.rm = T), semi_natur)) %>%
  left_join(land1000_2, by = 'Site') %>% 
  #extrapolating (.x are the flower estimates, .y is the area in hectares)
  mutate(Org_ex = CropBV1.x * CropBV1.y/100,
         Fl_BS11_ex = Flower_fieBS11.x * Flower_fieBS11.y/100,
         Fl_BS12_ex = Flower_fieBS12.x * Flower_fieBS12.y/100,
         Fl_BS2_ex = Flower_fieBS2.x * Flower_fieBS2.y/100,
         Fallow_ex = Fallow.x * Fallow.y/100,
         Other_AUM_ex = Other_AUM.x * Other_AUM.y/100,
         semi_natur_ex = semi_natur.x * semi_natur.y/100) %>%
  mutate(sum.fl = Org_ex + Fl_BS11_ex + Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex,
         SNH.ex = 0.5*Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex,
         Ann.fl.ex = Fl_BS11_ex + 0.5*Fl_BS12_ex ) %>%
  mutate(FL_per_agr = sum.fl / (Crop + CropBV1.y + Flower_fieBS2.y + Flower_fieBS12.y + Flower_fieBS11.y + Fallow.y +
                                  Other_AUM.y + semi_natur.y)) %>%
  select(Site, Org_ex, ends_with('.ex'), sum.fl, FL_per_agr)

plot(density(flower_cover1000$FL_per_agr))
dotchart(flower_cover1000$FL_per_agr, labels = flower_cover1000$Site)


land500 <- read.csv('Data/Landuse_Combee2022_500m.csv')
land500_2 <- land500 %>%
  mutate(area_ha = area1/10000) %>% 
  group_by(Landscape3,Combi6) %>%
  summarise(sum=sum(area_ha)) %>%
  drop_na(Landscape3) %>%
  pivot_wider(names_from = Combi6, values_from = sum) %>%
  replace(is.na(.), 0)%>%
  rename(Site = Landscape3) %>%
  mutate(semi_natur = semi_natur + Grassy_str)




## CALCULATING TOTAL FLOWER COVER AREA (EXTRAPOLATION)

flower_cover500 <- fl.cv2 %>% filter(Run == 2) %>% group_by(Site, Transect_type2) %>% summarise(flower_mean = mean(flcv.m))%>% 
  pivot_wider(names_from = Transect_type2, values_from = flower_mean) %>% ungroup() %>%
  #taking mean across all sites in case of NAs
  mutate(Flower_fieBS2 = ifelse(is.na(Flower_fieBS2), mean(Flower_fieBS2, na.rm = T), Flower_fieBS2),
         Flower_fieBS12 = ifelse(is.na(Flower_fieBS12), mean(Flower_fieBS12, na.rm = T), Flower_fieBS12),
         Flower_fieBS11 = ifelse(is.na(Flower_fieBS11), mean(Flower_fieBS11, na.rm = T), Flower_fieBS11),
         Other_AUM = ifelse(is.na(Other_AUM), mean(Other_AUM, na.rm = T), Other_AUM),
         Fallow = ifelse(is.na(Fallow), mean(Fallow, na.rm = T), Fallow),
         CropBV1 = ifelse(is.na(CropBV1), mean(CropBV1, na.rm = T), CropBV1),
         semi_natur = ifelse(is.na(semi_natur), mean(semi_natur, na.rm = T), semi_natur)) %>%
  left_join(land500_2, by = 'Site') %>% 
  #extrapolating (.x are the flower estimates, .y is the area in hectares)
  mutate(Org_ex = CropBV1.x * CropBV1.y/100,
         Fl_BS11_ex = Flower_fieBS11.x * Flower_fieBS11.y/100,
         Fl_BS12_ex = Flower_fieBS12.x * Flower_fieBS12.y/100,
         Fl_BS2_ex = Flower_fieBS2.x * Flower_fieBS2.y/100,
         Fallow_ex = Fallow.x * Fallow.y/100,
         Other_AUM_ex = Other_AUM.x * Other_AUM.y/100,
         semi_natur_ex = semi_natur.x * semi_natur.y/100) %>%
  mutate(sum.fl = Org_ex + Fl_BS11_ex + Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex,
         SNH.ex = 0.5*Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex,
         Ann.fl.ex = Fl_BS11_ex + 0.5*Fl_BS12_ex ) %>%
  mutate(FL_per_agr = sum.fl / (Crop + CropBV1.y + Flower_fieBS2.y + Flower_fieBS12.y + Flower_fieBS11.y + Fallow.y +
                                  Other_AUM.y + semi_natur.y)) %>%
  select(Site, Org_ex, ends_with('.ex'), sum.fl, FL_per_agr)


### hb abundance
hb.ab <- read_csv('Data/hb_abundance2022.csv') %>% rename(Date = date, Site = Landscape_ID, Transect_type = Transect_tpe) %>% mutate(Date = as.Date(Date, '%m/%d/%Y')) %>% group_by(Run, Site, Transect_type) %>% summarise(hb.m = mean(Honeybee), Date = max(Date)) %>% 
  left_join(hm, by = 'Transect_type') %>% left_join(hm2, by = 'Transect_type') %>% mutate(AES = ifelse(Transect_type == 'grassland', 'semi_natur', AES), Transect_type2 = ifelse(Transect_type == 'grassland', 'semi_natur', Transect_type2))
hb.ab2 <- hb.ab %>% filter(Run == 2) %>% group_by(Transect_type2, Site) %>% summarise(hb.m = mean(hb.m)) %>% 
  pivot_wider(names_from = Transect_type2, values_from = hb.m) 

honeybees1000 <- hb.ab2 %>% ungroup() %>%
  mutate(Flower_fieBS2 = ifelse(is.na(Flower_fieBS2), mean(Flower_fieBS2, na.rm = T), Flower_fieBS2),
         Flower_fieBS12 = ifelse(is.na(Flower_fieBS12), mean(Flower_fieBS12, na.rm = T), Flower_fieBS12),
         Flower_fieBS11 = ifelse(is.na(Flower_fieBS11), mean(Flower_fieBS11, na.rm = T), Flower_fieBS11),
         Other_AUM = ifelse(is.na(Other_AUM), mean(Other_AUM, na.rm = T), Other_AUM),
         Fallow = ifelse(is.na(Fallow), mean(Fallow, na.rm = T), Fallow),
         CropBV1 = ifelse(is.na(CropBV1), mean(CropBV1, na.rm = T), CropBV1),
         semi_natur = ifelse(is.na(semi_natur), mean(semi_natur, na.rm = T), semi_natur)) %>%
  left_join(land1000_2, by = 'Site') %>% 
  mutate(Org_ex = CropBV1.x * CropBV1.y/0.02,
         Fl_BS11_ex = Flower_fieBS11.x * Flower_fieBS11.y/0.02,
         Fl_BS12_ex = Flower_fieBS12.x * Flower_fieBS12.y/0.02,
         Fl_BS2_ex = Flower_fieBS2.x * Flower_fieBS2.y/0.02,
         Fallow_ex = Fallow.x * Fallow.y/0.02,
         Other_AUM_ex = Other_AUM.x * Other_AUM.y/0.02,
         semi_natur_ex = semi_natur.x * semi_natur.y/0.02
  ) %>%
  mutate(sum.hb = Org_ex + Fl_BS11_ex + Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex,
         SNH.ex = 0.5*Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex,
         Ann.fl.ex = Fl_BS11_ex + 0.5*Fl_BS12_ex ) %>%
  mutate(HB_per_agr = sum.hb / (Crop + CropBV1.y + Flower_fieBS2.y + Flower_fieBS12.y + Flower_fieBS11.y + Fallow.y +
                                  Other_AUM.y + semi_natur.y)) %>%
  select(Site, Org_ex, ends_with('.ex'), sum.hb, HB_per_agr)

plot(density(honeybees1000$HB_per_agr))
dotchart(honeybees1000$HB_per_agr, labels = honeybees1000$Site)
dotchart(honeybees1000$sum.hb, labels = honeybees1000$Site)

honeybees500 <- hb.ab2 %>% ungroup() %>%
  mutate(Flower_fieBS2 = ifelse(is.na(Flower_fieBS2), mean(Flower_fieBS2, na.rm = T), Flower_fieBS2),
         Flower_fieBS12 = ifelse(is.na(Flower_fieBS12), mean(Flower_fieBS12, na.rm = T), Flower_fieBS12),
         Flower_fieBS11 = ifelse(is.na(Flower_fieBS11), mean(Flower_fieBS11, na.rm = T), Flower_fieBS11),
         Other_AUM = ifelse(is.na(Other_AUM), mean(Other_AUM, na.rm = T), Other_AUM),
         Fallow = ifelse(is.na(Fallow), mean(Fallow, na.rm = T), Fallow),
         CropBV1 = ifelse(is.na(CropBV1), mean(CropBV1, na.rm = T), CropBV1),
         semi_natur = ifelse(is.na(semi_natur), mean(semi_natur, na.rm = T), semi_natur)) %>%
  left_join(land500_2, by = 'Site') %>% 
  mutate(Org_ex = CropBV1.x * CropBV1.y/0.02,
         Fl_BS11_ex = Flower_fieBS11.x * Flower_fieBS11.y/0.02,
         Fl_BS12_ex = Flower_fieBS12.x * Flower_fieBS12.y/0.02,
         Fl_BS2_ex = Flower_fieBS2.x * Flower_fieBS2.y/0.02,
         Fallow_ex = Fallow.x * Fallow.y/0.02,
         Other_AUM_ex = Other_AUM.x * Other_AUM.y/0.02,
         semi_natur_ex = semi_natur.x * semi_natur.y/0.02
  ) %>%
  mutate(sum.hb = Org_ex + Fl_BS11_ex + Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex,
         SNH.ex = 0.5*Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex,
         Ann.fl.ex = Fl_BS11_ex + 0.5*Fl_BS12_ex ) %>%
  mutate(HB_per_agr = sum.hb / (Crop + CropBV1.y + Flower_fieBS2.y + Flower_fieBS12.y + Flower_fieBS11.y + Fallow.y +
                                  Other_AUM.y + semi_natur.y)) %>%
  select(Site, Org_ex, ends_with('.ex'), sum.hb, HB_per_agr)

### different approach with pooling at the AES level

#flower_cover2 <- fl.cv2 %>% filter(Run == 2) %>% group_by(Site, Transect_type2, AES) %>% summarise(flower_mean = mean(flcv.m)) %>% 
#  group_by(Site, AES) %>% summarise(flower_mean = mean(flower_mean)) %>% pivot_wider(names_from = AES, values_from = flower_mean) %>% ungroup() %>%
#  mutate(Flower = ifelse(is.na(Flower), mean(Flower, na.rm = T), Flower),
#         Org.farm = ifelse(is.na(Org.farm), mean(Org.farm, na.rm = T), Org.farm),
#         SNH = ifelse(is.na(SNH), mean(SNH, na.rm = T), SNH))

#mean(flower_cover$Flower, na.rm = T)

#land1000 <- subset(land.all, radius == '1000m')

#flower1000 <- land1000 %>% select(Site, Total.fl, Org.farm, SNH.nf, Crop) %>% 
#  rename(Org.farm.h = Org.farm, Total.fl.h = Total.fl, SNH.nf.h = SNH.nf) %>% left_join(flower_cover2, by = 'Site') %>%
#  mutate(Org.fl = Org.farm.h * Org.farm/100,
#         Flower.fl = Total.fl.h * Flower/100,
#         SNH.fl = SNH.nf.h * SNH/100)
#flower1000$sum.fl <- rowSums(flower100_AES[9:11], na.rm = T)
#flower1000$FL_per_agr <- with(flower100_AES, sum.fl/(Org.farm.h + SNH.nf.h + Total.fl.h + Crop))


#honeybees1000 <- land1000 %>% select(Site, Org.farm, Total.fl, SNH.nf, Crop) %>% 
#  rename(Org.farm.h = Org.farm) %>%
#  left_join(hb.ab2, by = 'Site') %>%
#  mutate(Org.hb = Org.farm.h * Org.farm / 0.02,
#         Flower.hb = Total.fl * Flower / 0.02,
#         SNH.hb = SNH.nf * SNH / 0.02) 
#honeybees1000$sum.hb <- rowSums(honeybees1000[9:11], na.rm = T)
#honeybees1000$HB_per_agr <- with(honeybees1000, sum.hb/(Org.farm.h + SNH.nf + Total.fl + Crop))
