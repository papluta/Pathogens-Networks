library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)
source('02_flower_honeybee_extrapolation.R')

### DATA ###
hind <- read.csv('Data/hind_2022.csv') # honey bee pathogens
bind <- read.csv('Data/bind_2022.csv') # B. lapidarius bee pathogens
pind <- read.csv('Data/pind_2022.csv') # B. pascuorum pathogens - not finished!
wind <- read.csv('Data/wind_2022.csv') # other wild bee pathogens
q.norm <- read.csv('Data/abs_quant_parasites_2022.csv') # normalized absolute quantification
dens <- read.csv('Data/density.csv') # which site is density high
load('Data/networks_z.RData')
#load('Data/240617_landuse2022.RData') 
#load('Data/240617landscape_metrics2022.RData')


bb.load <- q.norm %>% filter(!grepl('H', Sample)) %>% filter(!grepl('S', Sample)) %>% mutate(Site = ifelse(grepl('B', Sample), sub('(B.*)','', Sample), sub('(P.*)','', Sample))) %>% group_by(Site) %>% summarise(
  dwvb.bb.q = mean(DWVB.abs, na.rm = T),
  bqcv.bb.q = mean(BQCV.abs, na.rm = T),
  abpv.bb.q = mean(ABPV.abs, na.rm = T),
  sbv.bb.q = mean(SBV.abs, na.rm = T),
)
bb.load[is.na(bb.load)] <- 0

hb.load <- q.norm %>% filter(grepl('H', Sample)) %>% mutate(Site = sub('(H.*)','', Sample)) %>% group_by(Site) %>% summarise(
  dwvb.hb.q = mean(DWVB.abs, na.rm = T),
  bqcv.hb.q = mean(BQCV.abs, na.rm = T),
  abpv.hb.q = mean(ABPV.abs, na.rm = T),
  sbv.hb.q = mean(SBV.abs, na.rm = T),
)
hb.load[is.na(hb.load)] <- 0

q.norm[is.na(q.norm)] <- 0 # changing undetected Ct to 0

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

bb.prev <- bind2 %>% group_by(Site) %>% summarise(dwvb.bb = mean(dwvb),
                                                  bqcv.bb = mean(bqcv),
                                                  abpv.bb = mean(abpv),
                                                  sbv.bb = mean(sbv))

data_500 <- all2 %>% rename(Sample = Sample.ID) %>% left_join(q.norm, by = 'Sample') %>% 
  left_join(dens, by = 'Site') %>%
  #left_join(land.all %>% filter(radius == '1000m'), by = 'Site') %>%
  #mutate(across(Ann.fl:AES, function(x) round(x / ((3.14 * 1000^2)/10000), digits = 3), .names = '{col}.p')) %>%
  #left_join(land_metrics1000, by = 'Site') %>%
  mutate(Density = as.factor(Density)) %>%
  filter(Species != "" & Species != 'NA' & Species != 'Sipha flava' & Species != 'Bombus sylvarum' & Species != 'Oedogonium sp. BN3'
         & Species != 'Megalocoleus molliculus' & Species != 'Orasema occidentalis' & Species != 'Orisarma intermedium' & Species != 'Lindenius albilabris') %>%
  left_join(networks.z3, by = 'Site') %>% 
  #left_join(fl.cv2, by = 'Site') %>%
  #left_join(network_parameters_species3, by = join_by('Site', 'Species')) %>%
  left_join(flower_cover500 %>% select(Site, sum.fl, FL_per_agr), by = 'Site') %>% 
  left_join(honeybees500 %>% select(Site, sum.hb, HB_per_agr), by = 'Site') %>%
  left_join(hb.prev, by = 'Site') %>% left_join(hb.load, by = 'Site') %>% 
  # calculating infection exposure 1st way: prevalence * abundance
  mutate(across(dwvb.hb:sbv.hb, function(x) log(x * HB_per_agr + 1), .names = '{col}.agr')) %>%
  # calculating infection exposure 2st way: prevalence * abundance * viral load
  mutate(dwvb.f = log(dwvb.hb * HB_per_agr * dwvb.hb.q + 1),
         bqcv.f = log(bqcv.hb * HB_per_agr * bqcv.hb.q+1),
         abpv.f = log(abpv.hb * HB_per_agr * abpv.hb.q+1),
         sbv.f = log(sbv.hb * HB_per_agr * sbv.hb.q+1)) %>% 
  #filter(Site != 'WM630') %>%
  mutate(HB_per_agr = log(HB_per_agr), FL_per_agr = FL_per_agr * 100) %>%
  left_join(bumblebees500 %>% select(Site, sum.bb, BB_per_agr), by = 'Site') %>% left_join(bb.load, by = 'Site') %>% left_join(bb.prev, by = 'Site') %>%
  mutate(dwvb.bb.f = log(dwvb.bb * BB_per_agr * dwvb.bb.q + 1),
         bqcv.bb.f = log(bqcv.bb * BB_per_agr * bqcv.bb.q+1),
         abpv.bb.f = log(abpv.bb * BB_per_agr * abpv.bb.q+1),
         sbv.bb.f = log(sbv.bb * BB_per_agr * sbv.bb.q+1)) %>%
  mutate(Closeness.z = ifelse(is.na(Closeness.z), -3, Closeness.z))
  
data_1000 <- all2 %>% rename(Sample = Sample.ID) %>% left_join(q.norm, by = 'Sample') %>% 
  left_join(dens, by = 'Site') %>%
  #left_join(land.all %>% filter(radius == '1000m'), by = 'Site') %>%
  #mutate(across(Ann.fl:AES, function(x) round(x / ((3.14 * 1000^2)/10000), digits = 3), .names = '{col}.p')) %>%
  #left_join(land_metrics1000, by = 'Site') %>%
  mutate(Density = as.factor(Density)) %>%
  filter(Species != "" & Species != 'NA' & Species != 'Sipha flava' & Species != 'Bombus sylvarum' & Species != 'Oedogonium sp. BN3'
         & Species != 'Megalocoleus molliculus' & Species != 'Orasema occidentalis' & Species != 'Orisarma intermedium' & Species != 'Lindenius albilabris') %>%
  left_join(networks.z3, by = 'Site') %>% 
  #left_join(fl.cv2, by = 'Site') %>%
  #left_join(network_parameters_species3, by = join_by('Site', 'Species')) %>%
  left_join(flower_cover1000 %>% select(Site, sum.fl, FL_per_agr), by = 'Site') %>% 
  left_join(honeybees1000 %>% select(Site, sum.hb, HB_per_agr), by = 'Site') %>%
  left_join(hb.prev, by = 'Site') %>% left_join(hb.load, by = 'Site') %>% 
  # calculating infection exposure 1st way: prevalence * abundance
  mutate(across(dwvb.hb:sbv.hb, function(x) log(x * HB_per_agr + 1), .names = '{col}.agr')) %>%
  # calculating infection exposure 2st way: prevalence * abundance * viral load
  mutate(dwvb.f = log(dwvb.hb * HB_per_agr * dwvb.hb.q + 1),
         bqcv.f = log(bqcv.hb * HB_per_agr * bqcv.hb.q+1),
         abpv.f = log(abpv.hb * HB_per_agr * abpv.hb.q+1),
         sbv.f = log(sbv.hb * HB_per_agr * sbv.hb.q+1)) %>% 
  #filter(Site != 'WM630') %>%
  mutate(HB_per_agr = log(HB_per_agr), FL_per_agr = FL_per_agr * 100) %>%
  left_join(bumblebees500 %>% select(Site, sum.bb, BB_per_agr), by = 'Site') %>% left_join(bb.load, by = 'Site') %>% left_join(bb.prev, by = 'Site') %>%
  mutate(dwvb.bb.f = log(dwvb.bb * BB_per_agr * dwvb.bb.q + 1),
         bqcv.bb.f = log(bqcv.bb * BB_per_agr * bqcv.bb.q+1),
         abpv.bb.f = log(abpv.bb * BB_per_agr * abpv.bb.q+1),
         sbv.bb.f = log(sbv.bb * BB_per_agr * sbv.bb.q+1)) %>%
  mutate(Closeness.z = ifelse(is.na(Closeness.z), -3, Closeness.z))

#dif <- data %>% select(Site, Species, dwvb) %>% full_join(network_parameters_species3, by = join_by('Site', 'Species'))

#no.match.Pat <- subset(dif, is.na(closeness))
#no.match.Kat <- subset(dif, is.na(dwvb))

  
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


data %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% mutate(Group = sub('bp','bb',Group)) %>% 
  group_by(Virus, Density, Group) %>% summarise(Prev = mean(Presence)) %>% mutate(Group = recode_factor(Group, 'hb' = 'Honey bee', 
                                                                                                        'bb' = 'Bumble bee', 
                                                                                                        'wb' = 'Other wild bee')) %>%
  mutate(Virus = as.factor(Virus)) %>% 
  ggplot(aes(Virus, Prev, fill = as.factor(Density))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  scale_fill_manual(values = c('#ebb93a', '#c9183a'))+
  labs(title = NULL, x = NULL, y = NULL)+
  theme_bw(base_size = 18)+
  facet_wrap(~Group, nrow = 1, ncol = 3)+
  theme(legend.position = 'none')+
  ylim(0,1)

data %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% mutate(Group = sub('bp','bb',Group)) %>% 
  group_by(Virus, Group) %>% summarise(Prev = mean(Presence)) %>% mutate(Group = recode_factor(Group, 'hb' = 'Honey bee', 
                                                                                                        'bb' = 'Bumble bee', 
                                                                                                        'wb' = 'Other wild bee')) %>%
  mutate(Virus = as.factor(Virus)) %>% 
  ggplot(aes(Virus, Prev, fill = as.factor(Group))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  scale_fill_manual(values = c('#d55e00', '#009e73', '#f0e442'))+
  labs(title = NULL, x = NULL, y = NULL, fill = 'Bee group')+
  theme_bw(base_size = 18)+
  scale_x_discrete(labels = c('abpv' = 'ABPV', 'bqcv' = 'BQCV','dwvb'= 'DWV-B', 'sbv' = 'SBV'))+
  #facet_wrap(~Group, nrow = 1, ncol = 3)+
  theme(axis.text = element_text(color = 'black'))+
  ylim(0,1)

ggsave(filename = 'Data/Models/prev2.png', height = 3, width = 8)

pprev <- data %>% pivot_longer(cols = dwvb:sbv, names_to = 'Virus', values_to = 'Presence') %>% mutate(Group = sub('bp','bb',Group)) %>% 
  group_by(Virus, Group, Site) %>% summarise(Prev = mean(Presence)) %>% mutate(Group = recode_factor(Group, 'hb' = 'Honey bee', 
                                                                                               'bb' = 'Bumble bee', 
                                                                                               'wb' = 'Other wild bee')) %>%
  mutate(Group = factor(Group, levels = c('Honey bee', 'Bumble bee', 'Other wild bee'))) %>%
  mutate(Virus = as.factor(Virus))
  
pprev %>% group_by(Virus, Group) %>% summarise(Prev2 = mean(Prev), sd = sd(Prev)) %>%
  ggplot(aes(Virus, Prev2, fill = as.factor(Group))) +
  geom_bar(stat = 'identity', position = position_dodge2(width = 0.3), col = 'black')+
  #geom_point(data = pprev, aes(Virus, Prev, group = Group), col = 'black', fill = 'white', shape = 21, position = position_dodge2(width = 0.9))+
  geom_errorbar(aes(ymin = Prev2 - sd, ymax = Prev2 + sd), position = position_dodge2(width = 0.3))+
  scale_fill_manual(values = c('#d55e00', '#009e73', '#f0e442'))+
  labs(title = NULL, x = NULL, y = NULL, fill = 'Bee group')+
  theme_bw(base_size = 18)+
  scale_x_discrete(labels = c('abpv' = 'ABPV', 'bqcv' = 'BQCV','dwvb'= 'DWV-B', 'sbv' = 'SBV'))+
  #facet_wrap(~Group, nrow = 1, ncol = 3)+
  theme(axis.text = element_text(color = 'black'))#+
  ylim(0,1)





data[1:13] %>% pivot_longer(cols = DWVB.abs:SBV.abs, names_to = 'Virus', values_to = 'Load') %>% mutate(Virus = sub('.abs', '', Virus)) %>% mutate(Group = sub('bp','bb',Group)) %>% 
  mutate(Group = factor(Group, levels = c('hb', 'bb', 'wb'))) %>% mutate(Group = recode_factor(Group, 'hb' = 'Honey bee', 
                                                                                                     'bb' = 'Bumble bee', 
                                                                                                     'wb' = 'Other wild bee')) %>%
  filter(Load != 0) %>%
  ggplot(aes(Virus, log10(Load+1), fill = as.factor(Group))) +
  geom_boxplot() +
  theme_bw(base_size = 16)+
  theme()+
  scale_fill_manual(values = c('#d55e00', '#009e73', '#f0e442'))+
  labs(x = NULL, y = 'Log-viral load', fill = 'Bee group')

ggsave(filename = 'Data/Models/load2.png', height = 3, width = 8)
