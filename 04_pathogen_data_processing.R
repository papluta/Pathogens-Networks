library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)

### DATA ###
hind.2022 <- read.csv('Data/hb_virus_2022.csv') # honey bee pathogens
lind.2022 <- read.csv('Data/bl_virus_2022.csv') # B. lapidarius bee pathogens
pind.2022 <- read.csv('Data/bp_virus_2022.csv') # B. pascuorum pathogens
wind.2022 <- read.csv('Data/wb_virus_2022.csv') # other wild bee pathogens

hind.2022.2 = hind.2022 %>% filter(!is.na(ACTIN) & !is.na(Species)) %>% mutate(Year = 2022)
lind.2022.2 = lind.2022 %>% filter(!is.na(ACTIN) & !is.na(Species)) %>% mutate(Year = 2022)
pind.2022.2 = pind.2022 %>% filter(!is.na(ACTIN) & !is.na(Species)) %>% mutate(Year = 2022)  
wind.2022.2 = wind.2022 %>% filter(!is.na(X28S) & !is.na(Species) & Analysis_done == 'Y') %>% mutate(Year = 2022) %>% select(!Analysis_done)

hind.2021 <- read.csv('Data/hb_virus_2021.csv') # honey bee pathogens
bind.2021 <- read.csv('Data/bb_virus_2021.csv') # B. lapidarius bee pathogens
wind.2021 <- read.csv('Data/wb_virus_2021.csv') # other wild bee pathogens

hind.2021.2 = hind.2021 %>% filter(!is.na(ACTIN) & !is.na(Species) & Analysis_done == 'Y') %>% mutate(Year = 2021) %>% select(!Analysis_done)
bind.2021.2 = bind.2021 %>% filter(!is.na(ACTIN) & !is.na(Species) & Analysis_done == 'Y') %>% mutate(Year = 2021) %>% select(!Analysis_done)
wind.2021.2 = wind.2021 %>% filter(!is.na(X28S) & !is.na(Species) & Analysis_done == 'Y') %>% mutate(Year = 2021) %>% select(!Analysis_done)


## assessing the cut-out value for Actin for honeybees:
hind.both = rbind(hind.2021.2, hind.2022.2)

#plot(density(hind.both$ACTIN))
h.m = mean(hind.both$ACTIN, na.rm = T)
h.sd = sd(hind.both$ACTIN, na.rm = T)
cutout.hb = h.m + 2*h.sd

## assessing the cut-out value for Actin for bumblebees 2021:
bind.both = bind.2021.2[bind.2021.2$Species != 'Bombus pascuorum',]
#plot(density(bind.both$ACTIN, na.rm = T))
b.m = mean(bind.both$ACTIN, na.rm = T)
b.sd = sd(bind.both$ACTIN, na.rm = T)
cutout.bb = b.m + 2*b.sd

## assessing the cut-out value for Actin for bumblebees 2022:
lind.both = lind.2022.2
#plot(density(lind.both$ACTIN, na.rm = T))
l.m = mean(lind.both$ACTIN, na.rm = T)
l.sd = sd(lind.both$ACTIN, na.rm = T)
cutout.bl = l.m + 2*l.sd

## assessing the cut-out value for Actin for bumblebees:
pind.both = rbind(bind.2021.2[bind.2021.2$Species == 'Bombus pascuorum',], pind.2022.2)
#plot(density(pind.both$ACTIN, na.rm = T))
p.m = mean(pind.both$ACTIN, na.rm = T)
p.sd = sd(pind.both$ACTIN, na.rm = T)
cutout.bp = p.m + 2*p.sd

## assessing the cut-out value for 28S for wild bees:
wind.both = rbind(wind.2021.2, wind.2022.2)

plot(density(wind.both$X28S, na.rm = T))
w.m = mean(wind.both$X28S, na.rm = T)
w.sd = sd(wind.both$X28S, na.rm = T)
cutout.wb = w.m + 2*w.sd

## cleaning up the molecular analysis

hind.both2 <- hind.both %>% filter(ACTIN < cutout.hb) %>% ## housekeeping gene below the cutout value
  mutate(dwvb = ifelse(!is.na(DWVB) & !is.na(DWVB.SD ), 1, 0), # changing undetected Ct to 0
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0)) %>% 
  mutate(DWVB.abs = ifelse(dwvb == 0 , 0, DWVB.abs), # removing load of samples that did not meet the postive criteria
         BQCV.abs = ifelse(bqcv == 0, 0, BQCV.abs),
         ABPV.abs = ifelse(abpv == 0, 0, ABPV.abs)) %>% 
  dplyr::select(-c(ACTIN:ABPV.SD)) %>%
  mutate(Group = 'hb')

bind.both2 <- bind.both %>% filter(ACTIN < cutout.bb) %>% ## housekeeping gene below the cutout value
  mutate(dwvb = ifelse(!is.na(DWVB) & !is.na(DWVB.SD ), 1, 0), # changing undetected Ct to 0
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0)) %>% 
  mutate(DWVB.abs = ifelse(dwvb == 0 , 0, DWVB.abs), # removing load of samples that did not meet the postive criteria
         BQCV.abs = ifelse(bqcv == 0, 0, BQCV.abs),
         ABPV.abs = ifelse(abpv == 0, 0, ABPV.abs)) %>% 
  dplyr::select(-c(ACTIN:ABPV.SD)) %>%
  mutate(Group = 'bb')

lind.both2 <- lind.both %>% filter(ACTIN < cutout.bl) %>% ## housekeeping gene below the cutout value
  mutate(dwvb = ifelse(!is.na(DWVB) & !is.na(DWVB.SD ), 1, 0), # changing undetected Ct to 0
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0)) %>% 
  mutate(DWVB.abs = ifelse(dwvb == 0 , 0, DWVB.abs), # removing load of samples that did not meet the postive criteria
         BQCV.abs = ifelse(bqcv == 0, 0, BQCV.abs),
         ABPV.abs = ifelse(abpv == 0, 0, ABPV.abs)) %>% 
  dplyr::select(-c(ACTIN:ABPV.SD)) %>%
  mutate(Group = 'bb')

pind.both2 <- pind.both %>% filter(ACTIN < cutout.bp) %>% ## housekeeping gene below the cutout value
  mutate(dwvb = ifelse(!is.na(DWVB) & !is.na(DWVB.SD ), 1, 0), # changing undetected Ct to 0
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0)) %>% 
  mutate(DWVB.abs = ifelse(dwvb == 0 , 0, DWVB.abs), # removing load of samples that did not meet the postive criteria
         BQCV.abs = ifelse(bqcv == 0, 0, BQCV.abs),
         ABPV.abs = ifelse(abpv == 0, 0, ABPV.abs)) %>% 
  dplyr::select(-c(ACTIN:ABPV.SD)) %>%
  mutate(Group = 'bb')


wind.both2 <- wind.both %>% filter(X28S < cutout.wb) %>% ## housekeeping gene below the cutout value
  mutate(dwvb = ifelse(!is.na(DWVB) & !is.na(DWVB.SD ), 1, 0), # changing undetected Ct to 0
         bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
         abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0)) %>% 
  mutate(DWVB.abs = ifelse(dwvb == 0 , 0, DWVB.abs), # removing load of samples that did not meet the postive criteria
         BQCV.abs = ifelse(bqcv == 0, 0, BQCV.abs),
         ABPV.abs = ifelse(abpv == 0, 0, ABPV.abs)) %>% 
  dplyr::select(-c(X28S:ABPV.SD)) %>%
  mutate(Group = 'wb') %>% mutate(Species = sub('Seladonia tumulorum', 'Halictus tumulorum', Species)) %>% group_by(Species) %>%
  mutate(n_Species = n()) %>% ungroup() %>% filter(n_Species > 2) %>% select(-n_Species)

data.pathogen.both <- rbind(hind.both2, bind.both2, lind.both2, pind.both2, wind.both2) #%>% mutate(across(DWVB.abs:SBV.abs, function(x) x/BUFFER))
levels(as.factor(data.pathogen.both$Species))

#write.csv(data.pathogen.both, file = 'Data/data_pathogen_both.csv', row.names = F)
