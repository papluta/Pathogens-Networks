library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)

# Load these files or re-run previous scripts (01 - 04)
# load('Data/land_extrapolated.RData')
# load('Data/morisita_zscore_both.RData') # bl means B. lapidarius
# load('Data/network_metrics_both_years.RData')
# data.pathogen.both = read.csv('Data/data_pathogen_both.csv')

dens <- read.csv('Data/density.csv')

## CALCULATING THE MEAN VALUE FOR NETWORK METRICS TO SUBSTITUTE IN CASE OF NAS

morisita.mean <- morisita.zscore.both %>% group_by(Species, Year) %>% summarise(Morisita.mean = mean(Morisita.z, na.rm = T))
morisita.mean2 <- morisita.zscore.both %>% group_by(Species) %>% summarise(Morisita.mean2 = mean(Morisita.z, na.rm = T))
morisita.bl.mean <- morisita.zscore.bl.both %>% group_by(Species, Year) %>% summarise(Morisita.bl.mean = mean(Morisita.z, na.rm = T))
morisita.bl.mean2 <- morisita.zscore.bl.both %>% group_by(Species) %>% summarise(Morisita.bl.mean2 = mean(Morisita.z, na.rm = T))

## CALCULATING PATHOGEN EXPOSURE TO THE MAIN HOST (A. MELLIFERA AND B. LAPIDARIUS)
hb.pat <- data.pathogen.both %>% filter(Species == 'Apis mellifera') %>% mutate(across(DWVB.abs:ABPV.abs, function(x) ifelse(x == 0, NA, x))) %>% 
  group_by(Site, Year) %>% 
  summarise(dwvb.p = mean(dwvb), dwvb.q.p = mean(DWVB.abs, na.rm = T),
            bqcv.p = mean(bqcv), bqcv.q.p = mean(BQCV.abs, na.rm = T),
            abpv.p = mean(abpv), abpv.q.p = mean(ABPV.abs, na.rm = T)) %>%
  mutate(across(everything(), function(x) ifelse(is.na(x), 0, x))) %>% 
  left_join(abundance.both %>% filter(Species == 'Apis mellifera'), by = c('Year','Site')) %>% rowwise() %>%
  left_join(flower.both %>% select(Site, Year, sum) %>% rename(FL.sum = sum), by = c('Site', 'Year')) %>% 
  mutate(dwvb.f = log(dwvb.p * dwvb.q.p * sum/FL.sum/10000 + 1),
         bqcv.f = log(bqcv.p * bqcv.q.p * sum/FL.sum/10000 + 1),
         abpv.f = log(abpv.p * abpv.q.p * sum/FL.sum/10000 + 1),
         hb.dens = sum/FL.sum/10000) %>% 
  select(Site, Year, dwvb.f, bqcv.f, hb.dens)


bl.pat <- data.pathogen.both %>% filter(Species == 'Bombus lapidarius') %>% mutate(across(DWVB.abs:ABPV.abs, function(x) ifelse(x == 0, NA, x))) %>% 
  group_by(Site, Year) %>% 
  summarise(dwvb.p = mean(dwvb), dwvb.q.p = mean(DWVB.abs, na.rm = T),
            bqcv.p = mean(bqcv), bqcv.q.p = mean(BQCV.abs, na.rm = T),
            abpv.p = mean(abpv), abpv.q.p = mean(ABPV.abs, na.rm = T)) %>%
  mutate(across(everything(), function(x) ifelse(is.na(x), 0, x))) %>% 
  left_join(flower.both %>% select(Site, Year, sum) %>% rename(FL.sum = sum), by = c('Site', 'Year')) %>% 
  left_join(abundance.both %>% filter(Species == 'Bombus lapidarius'), by = c('Year','Site')) %>% rowwise() %>%
  mutate(dwvb.bl.f = log(dwvb.p * sum/FL.sum/10000 * dwvb.q.p + 1),
         bqcv.bl.f = log(bqcv.p * sum/FL.sum/10000 * bqcv.q.p + 1),
         abpv.bl.f = log(abpv.p * sum/FL.sum/10000 * abpv.q.p + 1)) %>% 
  select(Site, Year, abpv.bl.f)


### COMBINING THE DATA
ab.site = abundance.both %>% group_by(Site, Year) %>% summarise(sum.ab = sum(sum))

data.both <- data.pathogen.both %>% 
  left_join(dens, by = 'Site') %>% mutate(Density = ifelse(Year == 2021, 0, Density)) %>%
  left_join(flower.both %>% select(Site, Year, FL.sum = sum), by = c('Year', 'Site')) %>%
  left_join(hb.pat, by = c('Site', 'Year')) %>%
  left_join(bl.pat, by = c('Site', 'Year')) %>%
  mutate(abpv.bl.f = ifelse(is.na(abpv.bl.f), mean(abpv.bl.f, na.rm = T), abpv.bl.f)) %>%
  left_join(ab.site, by = c('Site', 'Year')) %>%
  mutate(across(DWVB.abs:ABPV.abs, function(x) x/BUFFER)) %>%
  mutate(sum.ab.fl = log(sum.ab/FL.sum/10000)) %>%
  left_join(network.metrics.both %>% select(Connectance, Site, Year), by = c('Site', 'Year')) %>%
  left_join(morisita.zscore.both, by = c('Species', 'Site', 'Year')) %>%
  left_join(morisita.mean, by = c('Species', 'Year')) %>%
  left_join(morisita.mean2, by = 'Species') %>%
  mutate(Morisita.z = ifelse(is.na(Morisita.z), Morisita.mean, Morisita.z)) %>%
  mutate(Morisita.z = ifelse(is.na(Morisita.z), Morisita.mean2, Morisita.z)) %>%
  left_join(morisita.zscore.bl.both %>% rename(Morisita.bl.z = Morisita.z), by = c('Species', 'Site', 'Year')) %>%
  left_join(morisita.bl.mean, by = c('Species', 'Year')) %>%
  left_join(morisita.bl.mean2, by = 'Species') %>% 
  mutate(Morisita.bl.z = ifelse(is.na(Morisita.bl.z), Morisita.bl.mean, Morisita.bl.z)) %>%
  mutate(Morisita.bl.z = ifelse(is.na(Morisita.bl.z), Morisita.bl.mean2, Morisita.bl.z)) %>%
  select(-Morisita.mean, -Morisita.mean2, -Morisita.bl.mean, -Morisita.bl.mean2) %>%
  mutate(Year = as.factor(Year), Density = as.factor(Density)) %>%
  left_join(bee_rich %>% mutate(Year = as.factor(Year)), by = c('Site', 'Year')) %>%
  left_join(Flower_rich %>% mutate(Year = as.factor(Year)), by = c('Site', 'Year'))

#save(data.both, file = 'Data/final_dataframe.RData')


####  SPLITTING DATA INTO BEE GROUPS, SCALE AND FILTER THE OUTLIER
data.hb <- data.both %>% filter(Species == 'Apis mellifera') %>% ungroup() %>%
  filter(Morisita.bl.z > -20) %>% ## remove outlier site Nor1145
  mutate(across(FL.sum:Morisita.bl.z, ~scale(.)[,1], .names = "{.col}.s"))

data.bb.forhb <- data.both %>% filter(Group == 'bb') %>% ungroup() %>%
  filter(Morisita.z > -10) %>%  ## remove outlier site Nor1145 (B. lapidarius)
  mutate(across(FL.sum:Morisita.bl.z, ~scale(.)[,1], .names = "{.col}.s"))

data.bb.nolp <- data.both %>% filter(Group == 'bb') %>% filter(Species != 'Bombus lapidarius') %>% ungroup() %>%
  filter(Morisita.z > -10) %>%  ## remove outlier site Nor1145 (B. lapidarius)
  mutate(across(FL.sum:Morisita.bl.z, ~scale(.)[,1], .names = "{.col}.s"))

data.wb <- data.both %>% filter(Group == 'wb') %>% ungroup() %>%
  mutate(across(FL.sum:Morisita.bl.z, ~scale(.)[,1], .names = "{.col}.s"))  %>% drop_na(Morisita.bl.z) ## some


data.site <- data.both %>% ungroup() %>% 
  distinct(Site, Year, Density, Connectance, FL.sum, FL_rich, Bee_rich, sum.ab, sum.ab.fl, hb.dens) %>%
  mutate(FL.sum.s = scale(FL.sum)[,1], sum.ab.fl.s = scale(sum.ab.fl)[,1]) 

data.morisita.bl <- data.both %>% ungroup() %>% distinct(Site, Year, Density, Species, Morisita.bl.z, FL.sum, sum.ab, sum.ab.fl, FL_rich, Bee_rich) %>% filter(Species != 'Bombus lapidarius') %>%
  filter(Morisita.bl.z > -12) %>% mutate(across(Morisita.bl.z:sum.ab.fl, ~scale(.)[,1], .names = "{.col}.s"))
data.morisita.hb <- data.both %>% ungroup() %>% distinct(Site, Year, Density, Species, Morisita.z, FL.sum, sum.ab, sum.ab.fl, FL_rich, Bee_rich) %>% filter(Species != 'Apis mellifera') %>%
  filter(Morisita.z > -8) %>% mutate(across(Morisita.z:sum.ab.fl, ~scale(.)[,1], .names = "{.col}.s"))

### SOME SIMPLE RAW DATA VISUALISATION

# dotchart(data.morisita$Morisita.z)
# dotchart(data.site$Connectance)
# dotchart(data.site$abpv.bl.f)
# dotchart(unique(data.bb$Morisita.z), xlab = 'Niche overlap between\nA. mellifera and B. lapidarius', ylab = 'Site')
