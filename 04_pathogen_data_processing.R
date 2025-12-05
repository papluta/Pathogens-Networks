library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(cowplot)

### DATA ###
hind.2021 <- read.csv('Data/hb_virus_2021.csv') # honey bee pathogens
bind.2021 <- read.csv('Data/bb_virus_2021.csv') # B. lapidarius, B. pascuorum & B. terrestris pathogens
wind.2021 <- read.csv('Data/wb_virus_2021.csv') # other wild bee pathogens

hind.2021.2 = hind.2021 %>% filter(!is.na(ACTIN) & !is.na(Species) & Analysis_done == 'Y') %>% mutate(Year = 2021) %>% select(!Analysis_done)
bind.2021.2 = bind.2021 %>% filter(!is.na(ACTIN) & !is.na(Species) & Analysis_done == 'Y') %>% mutate(Year = 2021) %>% select(!Analysis_done)
wind.2021.2 = wind.2021 %>% filter(!is.na(X28S) & !is.na(Species) & Analysis_done == 'Y') %>% mutate(Year = 2021) %>% select(!Analysis_done)

lind.2021.2 <- bind.2021.2 %>% filter(Species == "Bombus lapidarius")
pind.2021.2 <- bind.2021.2 %>% filter(Species == "Bombus pascuorum")
tind.2021.2 <- bind.2021.2 %>% filter(Species == "Bombus terrestris")

hind.2022 <- read.csv('Data/hb_virus_2022.csv') # honey bee pathogens
lind.2022 <- read.csv('Data/bl_virus_2022.csv') # B. lapidarius  pathogens
pind.2022 <- read.csv('Data/bp_virus_2022.csv') # B. pascuorum pathogens
wind.2022 <- read.csv('Data/wb_virus_2022.csv') # other wild bee pathogens

hind.2022.2 = hind.2022 %>% filter(!is.na(ACTIN) & !is.na(Species)) %>% mutate(Year = 2022)
lind.2022.2 = lind.2022 %>% filter(!is.na(ACTIN) & !is.na(Species)) %>% mutate(Year = 2022)
pind.2022.2 = pind.2022 %>% filter(!is.na(ACTIN) & !is.na(Species)) %>% mutate(Year = 2022)  
wind.2022.2 = wind.2022 %>% filter(!is.na(X28S) & !is.na(Species) & Analysis_done == 'Y') %>% 
  mutate(Year = 2022) %>% select(!Analysis_done)

## cleaning up the molecular analysis

remove_low_quality <- function(data, bee_group, gene) {
  data %>% 
    rename(qual = gene) %>%
    filter(qual < mean(qual, na.rm = T) + 2*sd(qual, na.rm = T)) %>% #removing samples that are above 2x SD Ct
    mutate(dwvb = ifelse(!is.na(DWVB) & !is.na(DWVB.SD ), 1, 0), # changing undetected Ct to 0
           bqcv = ifelse(!is.na(BQCV) & !is.na(BQCV.SD ), 1, 0),
           abpv = ifelse(!is.na(ABPV) & !is.na(ABPV.SD ), 1, 0)) %>% 
    mutate(DWVB.abs = ifelse(dwvb == 0 , 0, DWVB.abs), # removing load of samples that did not meet the postive criteria
           BQCV.abs = ifelse(bqcv == 0, 0, BQCV.abs),
           ABPV.abs = ifelse(abpv == 0, 0, ABPV.abs)) %>% 
    dplyr::select(-c(qual:ABPV.SD)) %>%
    mutate(Group = bee_group)
}

hind.2021.3 <- remove_low_quality(hind.2021.2, "hb", "ACTIN")
hind.2022.3 <- remove_low_quality(hind.2022.2, "hb", "ACTIN")

pind.2021.3 <- remove_low_quality(pind.2021.2, "bb", "ACTIN")
pind.2022.3 <- remove_low_quality(pind.2022.2, "bb", "ACTIN")

tind.2021.3 <- remove_low_quality(tind.2021.2, "bb", "ACTIN")


lind.2021.3 <- remove_low_quality(lind.2021.2, "bb", "ACTIN")
lind.2022.3 <- remove_low_quality(lind.2022.2, "bb", "ACTIN")


wind.2021.3 <- remove_low_quality(wind.2021.2, "wb", "X28S") 
wind.2022.3 <- remove_low_quality(wind.2022.2, "wb", "X28S")


data.pathogen.both <- rbind(hind.2021.3, hind.2022.3, 
                            pind.2021.3, pind.2022.3, 
                            tind.2021.3, 
                            lind.2021.3, lind.2022.3,
                            wind.2021.3, wind.2022.3) %>%
  mutate(Species = sub('Seladonia tumulorum', 'Halictus tumulorum', Species)) %>% #harmonizing species names
  group_by(Species) %>%
  mutate(n_Species = n()) %>% 
  ungroup() %>% 
  filter(n_Species > 2) %>% #removing species with less than 3 individuals in the dataset
  select(-n_Species)

#write.csv(data.pathogen.both, file = 'Data/data_pathogen_both.csv', row.names = F)

## calculating number of low quality samples
total <- list(hind.2021.2, hind.2022.2, 
                            pind.2021.2, pind.2022.2, 
                            tind.2021.2, 
                            lind.2021.2, lind.2022.2,
                            wind.2021.2, wind.2022.2) %>% map(., nrow) %>%
  unlist()

filtered <- list(hind.2021.3, hind.2022.3, 
                  pind.2021.3, pind.2022.3, 
                  tind.2021.3, 
                  lind.2021.3, lind.2022.3,
                  wind.2021.3, wind.2022.3) %>% map(., nrow) %>%
  unlist()

sum(total) - sum(filtered)
