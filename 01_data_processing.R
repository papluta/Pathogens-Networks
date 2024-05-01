library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)

### DATA ###
hind <- read.csv('Data/hind_2022.csv') # honey bee pathogens
bind <- read.csv('Data/bind_2022.csv') # B. lapidarius bee pathogens
#pind <- read.csv('Data/pind_2022.csv') # B. pascuorum pathogens - not finished!
wind <- read.csv('Data/wind_2022.csv') # other wild bee pathogens
dens <- read.csv('Data/density.csv') # which site is density high
hb.ab <- read.csv('Data/hb_abundance2022.csv') # Kathrin's honey bee abundance 
load('Data/landuse2022.RData') # landuse composition at scales from 100 to 2000 m
load('Data/landscape_metrics2022.RData')
wind.sp <- read.csv('Data/WIND_species.csv') %>% select(Sample_ID, Species) # updated barcodes for wild bees
load('Data/abs_pathogen240212.RData') # normalized absolute quantification 

wind <- wind %>% select(-Species) %>% left_join(wind.sp, by = join_by('Sample.ID' == 'Sample_ID'))
q.norm[is.na(q.norm)] <- 0 # changing undetected Ct to 0

## current issues:
# - check absolute quantification for all samples
# - Shannon from vegan (shp) and Shannon from landscapemetrics (rasters) give different results - recheck the rasters
# - not sure if I calculated correctly the total flower cover (grassy strips!)

## cleaning up the molecular analysis

hind2 <- hind %>% filter(ACTIN < 27) %>% dplyr::select(-ACTIN) %>% ## housekeeping gene below 27
  mutate(dwvb = ifelse(!is.na(DWV.B) & !is.na(DWV.B.SD ), 1, 0), # changing undetected Ct to 0
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0),
         sbv = ifelse(!is.na(SBV) & !is.na(SBV.SD ), 1, 0)) %>% 
  dplyr::select(Site, Sample.ID, dwvb, bqcv, abpv, sbv) %>%
  mutate(Group = 'hb')

bind2 <- bind %>% filter(ACTIN < 27) %>% dplyr::select(-ACTIN) %>% 
  mutate(dwvb = ifelse(!is.na(DWV.B) & !is.na(DWV.B.SD ), 1, 0),
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0),
         sbv = ifelse(!is.na(SBV) & !is.na(SBV.SD ), 1, 0)) %>%
  dplyr::select(Site, Sample.ID, dwvb, bqcv, abpv, sbv) %>%
  mutate(Group = 'bb')

wind2 <- wind %>% filter(X28S < 27) %>% dplyr::select(-X28S) %>% 
  mutate(dwvb = ifelse(!is.na(DWV.B) & !is.na(DWV.B.SD ), 1, 0),
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0),
         sbv = ifelse(!is.na(SBV) & !is.na(SBV.SD ), 1, 0)) %>%
  dplyr::select(Site, Sample.ID, dwvb, bqcv, abpv, sbv) %>%
  mutate(Group = 'wb')

all2 <- rbind(hind2, bind2, wind2)

data <- all2 %>% rename(Sample = Sample.ID) %>% left_join(q.norm, by = 'Sample') %>% 
  left_join(land.all %>% filter(radius == '2000m'), by = 'Site') %>%
  left_join(Het.lm2000 %>% select(Iji, Ennd.Total_flowerstrip, Edge.dens), by = 'Site') %>%
  

                           