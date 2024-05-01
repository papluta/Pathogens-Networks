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
wind.sp <- read.csv('Data/WIND_species.csv') %>% select(Sample_ID, Species) # updated barcodes for wild bees
load('Data/abs_pathogen240212.RData') # normalized absolute quantification
dens <- read.csv('Data/density.csv') # which site is density high
hb.ab <- read.csv('Data/hb_abundance2022.csv') # Kathrin's honey bee abundance 
 
load('Data/landuse2022.RData') # landuse composition at scales from 200 to 2000 m
load('Data/landscape_metrics2022.RData') # Shannon index, Edge density, patch edge length, patch ENND, IJI, flower strips ENND at 500, 1000, 1500, and 2000 m radii
load('Data/fl_cover_extr.RData') # extrapolated flower cover

wind <- wind %>% select(-Species) %>% left_join(wind.sp, by = join_by('Sample.ID' == 'Sample_ID')) # updating the barcodes
q.norm[is.na(q.norm)] <- 0 # changing undetected Ct to 0

## current issues:
# - check absolute quantification for all samples
# - Shannon from vegan (shp) and Shannon from landscapemetrics (rasters) give different results - recheck the rasters
# - not sure if I calculated correctly the total flower cover (grassy strips!)
# - which flower cover estimation run use for pathogen data? currently using an average of 1st and 2nd
# - sites missing total flowerstrip Ennd - add 2000 m distance?

## cleaning up the molecular analysis

hind2 <- hind %>% filter(ACTIN < 27) %>% dplyr::select(-ACTIN) %>% ## housekeeping gene below 27
  mutate(dwvb = ifelse(!is.na(DWV.B) & !is.na(DWV.B.SD ), 1, 0), # changing undetected Ct to 0
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0),
         sbv = ifelse(!is.na(SBV) & !is.na(SBV.SD ), 1, 0)) %>% 
  dplyr::select(Site, Sample.ID, dwvb, bqcv, abpv, sbv) %>%
  mutate(Species = 'Apis mellifera', Group = 'hb')

bind2 <- bind %>% filter(ACTIN < 27) %>% dplyr::select(-ACTIN) %>% 
  mutate(dwvb = ifelse(!is.na(DWV.B) & !is.na(DWV.B.SD ), 1, 0),
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0),
         sbv = ifelse(!is.na(SBV) & !is.na(SBV.SD ), 1, 0)) %>%
  dplyr::select(Site, Sample.ID, dwvb, bqcv, abpv, sbv) %>%
  mutate(Species = 'Bombus lapidarius', Group = 'bb')

wind2 <- wind %>% filter(X28S < 27) %>% dplyr::select(-X28S) %>% 
  mutate(dwvb = ifelse(!is.na(DWV.B) & !is.na(DWV.B.SD ), 1, 0),
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0),
         sbv = ifelse(!is.na(SBV) & !is.na(SBV.SD ), 1, 0)) %>%
  dplyr::select(Site, Sample.ID, dwvb, bqcv, abpv, sbv, Species) %>%
  mutate(Group = 'wb')

all2 <- rbind(hind2, bind2, wind2)

data <- all2 %>% rename(Sample = Sample.ID) %>% left_join(q.norm, by = 'Sample') %>% 
  left_join(dens, by = 'Site') %>%
  left_join(land.all %>% filter(radius == '2000m'), by = 'Site') %>%
  left_join(Het[['2000m']] %>% select(Site, SH, Iji, Ennd.Total_flowerstrip, Edge.dens), by = 'Site') %>%
  left_join(flower.extr[['2000m']] %>% filter(Run !=3) %>% group_by(Site, Run) %>% summarise(fl.cv = sum(flower.cover)) %>%
              group_by(Site) %>% summarise(fl.cv = mean(fl.cv)), by = 'Site') %>%
  mutate(Density = as.factor(Density))
                           