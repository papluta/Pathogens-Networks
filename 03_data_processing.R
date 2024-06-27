library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
source('02_flower_cover.R')

### DATA ###
hind <- read.csv('Data/hind_2022.csv') # honey bee pathogens
bind <- read.csv('Data/bind_2022.csv') # B. lapidarius bee pathogens
pind <- read.csv('Data/pind_2022.csv') # B. pascuorum pathogens - not finished!
wind <- read.csv('Data/wind_2022.csv') # other wild bee pathogens
wind.sp <- read.csv('Data/WIND_species.csv') %>% dplyr::select(Sample_ID, Species) # updated barcodes for wild bees
load('Data/abs_pathogen240212.RData') # normalized absolute quantification
dens <- read.csv('Data/density.csv') # which site is density high
hb.ab <- read.csv('Data/hb_abundance2022.csv') # Kathrin's honey bee abundance
load('Data/240621network_parameters.RData')


load('Data/240617_landuse2022.RData') # landuse composition at scales from 200 to 2000 m
load('Data/240617landscape_metrics2022.RData') # Shannon index, Edge density, patch edge length, patch ENND, IJI, flower strips ENND at 500, 1000, 1500, and 2000 m radii
#load('Data/fl_cover_extr.RData') # extrapolated flower cover


wind <- wind %>% dplyr::select(-Species) %>% left_join(wind.sp, by = join_by('Sample.ID' == 'Sample_ID')) # updating the barcodes
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

pind2 <- pind %>% filter(B_ACTIN < 22) %>% dplyr::select(-B_ACTIN, -B_ACTIN.SD) %>% 
  mutate(dwvb = ifelse(!is.na(DWV.B) & !is.na(DWV.B.SD ), 1, 0),
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0),
         sbv = ifelse(!is.na(SBV) & !is.na(SBV.SD ), 1, 0)) %>%
  dplyr::select(Site, Sample.ID, dwvb, bqcv, abpv, sbv) %>%
  mutate(Species = 'Bombus pascuorum', Group = 'bp')

wind2 <- wind %>% filter(X28S < 27) %>% dplyr::select(-X28S) %>% 
  mutate(dwvb = ifelse(!is.na(DWV.B) & !is.na(DWV.B.SD ), 1, 0),
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0),
         sbv = ifelse(!is.na(SBV) & !is.na(SBV.SD ), 1, 0)) %>%
  dplyr::select(Site, Sample.ID, dwvb, bqcv, abpv, sbv, Species) %>%
  mutate(Group = 'wb')

all2 <- rbind(hind2, bind2, wind2, pind2)

hb.prev <- hind2 %>% group_by(Site) %>% summarise(dwvb.hb = mean(dwvb),
                                                  bqcv.hb = mean(bqcv),
                                                  abpv.hb = mean(abpv),
                                                  sbv.hb = mean(sbv))

##### KATHRINS DATA######

#flower_Kat <- read.csv('Data/fc_extrapol_2022_20240625.csv')
#flower_Kat2 <- flower_Kat %>% filter(Run == 2)
#flower_Kat2$totalFC_h <- flower_Kat2$totalFC / 10000

#flower_Kat2
#ann.fl <- data.frame(flower_Kat2$Landscape_ID, flower_Kat2$fcflowerf_stand/flower_Kat2$Annual_flower/100)
#org <- data.frame(flower_Kat2$Landscape_ID, flower_Kat2$fcorganic_stand/land_Kat$CropBV1/100)

#flower_cover <- fl.cv2 %>% filter(Run == 2) %>% group_by(Site, AES) %>% summarise(flower_mean = mean(flcv.m)) %>% pivot_wider(names_from = AES, values_from = flower_mean)
#print(flower_cover, n = Inf)

#flower_Pat <- flower1000 %>% mutate(totalFC_h_P = round(sum.fl / ((3.14 * 1000^2)/10000), digits = 3))
#hb_Kat <- read.csv('Data/hb_run2_extrapolated_2022_20240625.csv')
#land_Kat <- read.csv('Data/Landscape_variables_ComBee_2022_20230203.csv')

#comp <- land_Kat %>% select(Land_ID, CropBV1, Annual_flower, semi) %>% mutate(Land_ID = sub('Wm', 'WM', Land_ID)) %>% left_join(land1000 %>% select(Site, Org.farm, Ann.fl, SNH), by = join_by('Land_ID' == 'Site'))

#data_Kat <- flower_Kat2 %>% select(Landscape_ID, totalFC, totalFC_h, totalagrfcm2) %>% rename(Site = Landscape_ID, FL_per_agr = totalagrfcm2) %>%
#  left_join(hb_Kat %>% mutate(HB_per_agr = Total_honeybee_abu / agr_total_m2, totalHB_h = Total_honeybee_abu/10000) %>% select(Landscape_ID, Total_honeybee_abu, totalHB_h, HB_per_agr) %>% rename(Site = Landscape_ID, totalHB = Total_honeybee_abu), by = 'Site')


data <- all2 %>% rename(Sample = Sample.ID) %>% left_join(q.norm, by = 'Sample') %>% 
  left_join(dens, by = 'Site') %>%
  left_join(land.all %>% filter(radius == '1000m'), by = 'Site') %>%
  mutate(across(Ann.fl:AES, function(x) round(x / ((3.14 * 1000^2)/10000), digits = 3), .names = '{col}.p')) %>%
  left_join(land_metrics1000, by = 'Site') %>%
  mutate(Density = as.factor(Density)) %>%
  filter(Species != "" & Species != 'NA' & Species != 'Sipha flava' & Species != 'Bombus sylvarum' & Species != 'Oedogonium sp. BN3'
         & Species != 'Megalocoleus molliculus' & Species != 'Orasema occidentalis' & Species != 'Orisarma intermedium' & Species != 'Lindenius albilabris') %>%
  left_join(network_parameters2, by = 'Site') %>% rename(weighted.nested = `weighted NODF`) %>%
  #left_join(fl.cv2, by = 'Site') %>%
  left_join(network_parameters_species3, by = join_by('Site', 'Species')) %>%
  left_join(flower1000 %>% select(Site, sum.fl, FL_per_agr), by = 'Site') %>% left_join(honeybees1000%>% select(Site, sum.hb, HB_per_agr), by = 'Site') %>%
  left_join(hb.prev, by = 'Site') %>%
  mutate(across(dwvb.hb:sbv.hb, function(x) x * sum.hb, .names = '{col}.f')) %>%
  mutate(across(dwvb.hb:sbv.hb, function(x) x * HB_per_agr, .names = '{col}.agr'))
  
  

dif <- data %>% select(Site, Species, dwvb) %>% full_join(network_parameters_species3, by = join_by('Site', 'Species'))

no.match.Pat <- subset(dif, is.na(closeness))
no.match.Kat <- subset(dif, is.na(dwvb))

  
library(cowplot)

# prevalence

h <- hind2 %>% left_join(dens, by = 'Site') %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Virus, Density) %>% summarise(Prev = mean(Presence)) %>% mutate(Virus = factor(Virus, levels = c('dwvb', 'bqcv', 'abpv', 'sbv'))) %>%
  ggplot(aes(Virus, Prev, fill = as.factor(Density))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  labs(title = 'Honey bee', x = NULL, y = NULL)+
  theme_minimal()+
  theme(legend.position = 'none')+
  ylim(0,1)
b <- bind2 %>% left_join(dens, by = 'Site') %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Virus, Density) %>% summarise(Prev = mean(Presence)) %>% mutate(Virus = factor(Virus, levels = c('dwvb', 'bqcv', 'abpv', 'sbv'))) %>%
  ggplot(aes(Virus, Prev, fill = as.factor(Density))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  labs(title = 'B. lapidarius', x = NULL, y = NULL)+
  theme_minimal()+
  theme(legend.position = 'none')+
  ylim(0,1)
p <- pind2 %>% left_join(dens, by = 'Site') %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Virus, Density) %>% summarise(Prev = mean(Presence)) %>% mutate(Virus = factor(Virus, levels = c('dwvb', 'bqcv', 'abpv', 'sbv'))) %>%
  ggplot(aes(Virus, Prev, fill = as.factor(Density))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  labs(title = 'B. pascuorum', x = NULL, y = NULL)+
  theme_minimal()+
  theme(legend.position = 'none')+
  ylim(0,1)
w <- wind2 %>% left_join(dens, by = 'Site') %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Virus, Density) %>% summarise(Prev = mean(Presence)) %>% mutate(Virus = factor(Virus, levels = c('dwvb', 'bqcv', 'abpv', 'sbv'))) %>%
  ggplot(aes(Virus, Prev, fill = as.factor(Density))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  theme_minimal()+
  theme(legend.position = 'none')+
  labs(title = 'Wild bee', x = NULL, y = NULL)+
  ylim(0,1)

plot_grid(h,b,p,w)

# viral load
data[1:8] %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Virus, Group, Site) %>% summarise(Prev = mean(Presence)) %>% 
  mutate(Virus = factor(Virus, levels = c('dwvb', 'bqcv', 'abpv', 'sbv'))) %>%
  filter(Virus == 'abpv') %>%
  ggplot(aes(Site, Prev, fill = as.factor(Group)))+
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  #scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  labs(title = NULL, x = NULL, y = NULL)+
  theme_minimal()+
  theme(legend.position = 'none')+
  ylim(0,1)

data %>% filter(Group == 'hb') %>% filter(Site == 'Gos2') %>% dplyr::select(Sample, ABPV.abs)
data %>% filter(Group == 'bb') %>% filter(Site == 'Gos2') %>% dplyr::select(Sample, ABPV.abs)
data %>% filter(Group == 'bp') %>% filter(Site == 'Gos2') %>% dplyr::select(Sample, ABPV.abs)
data %>% filter(Group == 'wb') %>% filter(Site == 'Gos2') %>% dplyr::select(Sample, ABPV.abs)

h <- hind2 %>% left_join(dens, by = 'Site') %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Virus, Density) %>% summarise(Prev = mean(Presence)) %>% mutate(Virus = factor(Virus, levels = c('dwvb', 'bqcv', 'abpv', 'sbv'))) %>%
  ggplot(aes(Virus, Prev, fill = as.factor(Density))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  labs(title = 'Honey bee', x = NULL, y = NULL)+
  theme_minimal()+
  theme(legend.position = 'none')+
  ylim(0,1)
b <- bind2 %>% left_join(dens, by = 'Site') %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Virus, Density) %>% summarise(Prev = mean(Presence)) %>% mutate(Virus = factor(Virus, levels = c('dwvb', 'bqcv', 'abpv', 'sbv'))) %>%
  ggplot(aes(Virus, Prev, fill = as.factor(Density))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  labs(title = 'B. lapidarius', x = NULL, y = NULL)+
  theme_minimal()+
  theme(legend.position = 'none')+
  ylim(0,1)
p <- pind2 %>% left_join(dens, by = 'Site') %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Virus, Density) %>% summarise(Prev = mean(Presence)) %>% mutate(Virus = factor(Virus, levels = c('dwvb', 'bqcv', 'abpv', 'sbv'))) %>%
  ggplot(aes(Virus, Prev, fill = as.factor(Density))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  labs(title = 'B. pascuorum', x = NULL, y = NULL)+
  theme_minimal()+
  theme(legend.position = 'none')+
  ylim(0,1)
w <- wind2 %>% left_join(dens, by = 'Site') %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% 
  group_by(Virus, Density) %>% summarise(Prev = mean(Presence)) %>% mutate(Virus = factor(Virus, levels = c('dwvb', 'bqcv', 'abpv', 'sbv'))) %>%
  ggplot(aes(Virus, Prev, fill = as.factor(Density))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  theme_minimal()+
  theme(legend.position = 'none')+
  labs(title = 'Wild bee', x = NULL, y = NULL)+
  ylim(0,1)

plot_grid(h,b,p,w)

data[1:13] %>% pivot_longer(cols = DWVB.abs:SBV.abs, names_to = 'Virus', values_to = 'Load') %>% mutate(Virus = sub('.abs', '', Virus)) %>%
  mutate(Group = factor(Group, levels = c('hb', 'bb', 'bp', 'wb'))) %>% mutate(Group = recode_factor(Group, 'hb' = 'Apis mellifera', 
                                                                                                     'bb' = 'B. lapidarius', 
                                                                                                     'bp' = 'B. pascuorum', 
                                                                                                     'wb' = 'Wild bee'),
                                                                               Virus = factor(Virus, levels = c('DWVB', 'BQCV', 'ABPV', 'SBV'))) %>%
  filter(Load != 0) %>%
  ggplot(aes(Virus, log10(Load+1), fill = as.factor(Density))) +
  geom_boxplot() +
  facet_wrap(~Group)+
  theme_minimal(base_size = 16)+
  theme(legend.position = 'none')+
  scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  labs(x = NULL, y = 'Log-viral load')

