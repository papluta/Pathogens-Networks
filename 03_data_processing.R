library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)

### DATA ###
hind <- read.csv('Data/hind_2022.csv') # honey bee pathogens
bind <- read.csv('Data/bind_2022.csv') # B. lapidarius bee pathogens
pind <- read.csv('Data/pind_2022.csv') # B. pascuorum pathogens - not finished!
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

data <- all2 %>% rename(Sample = Sample.ID) %>% left_join(q.norm, by = 'Sample') %>% 
  left_join(dens, by = 'Site') %>%
  left_join(land.all %>% filter(radius == '2000m'), by = 'Site') %>%
  left_join(Het[['2000m']] %>% select(Site, SH, Iji, Ennd.Total_flowerstrip, Edge.dens), by = 'Site') %>%
  left_join(flower.extr[['2000m']] %>% filter(Run !=3) %>% group_by(Site, Run) %>% summarise(fl.cv = sum(flower.cover)) %>%
              group_by(Site) %>% summarise(fl.cv = mean(fl.cv)), by = 'Site') %>%
  mutate(Density = as.factor(Density)) %>%
  filter(Species != "" & Species != 'NA' & Species != 'Sipha flava' & Species != 'Bombus sylvarum' & Species != 'Oedogonium sp. BN3'
         & Species != 'Megalocoleus molliculus' & Species != 'Orasema occidentalis' & Species != 'Orisarma intermedium' & Species != 'Lindenius albilabris')

hm <- pind2 %>% left_join(dens, by = 'Site') %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% group_by(Virus, Density) %>% summarise(Prev = mean(Presence))

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

data %>% filter(Group == 'hb') %>% filter(Site == 'Gos2') %>% select(Sample, ABPV.abs)
data %>% filter(Group == 'bb') %>% filter(Site == 'Gos2') %>% select(Sample, ABPV.abs)
data %>% filter(Group == 'bp') %>% filter(Site == 'Gos2') %>% select(Sample, ABPV.abs)
data %>% filter(Group == 'wb') %>% filter(Site == 'Gos2') %>% select(Sample, ABPV.abs)

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
