library(sf)

load('Data/traits.RData')

source("03_networks.R") # reruns scripts 01-03
source("04_pathogen_data_processing.R")

dens <- read.csv('Data/density.csv')
coord <- read.csv('Data/coordinates.csv')


## transforming degree coordinates to km

coord_sf <- st_as_sf(coord, coords = c("Long", "Lat"), crs = 4326)
coord_sf_utm <- st_transform(coord_sf, crs = 32632)
dist_km <- st_coordinates(coord_sf_utm)  / 1000
coord2 <- data.frame(Site = coord_sf_utm$Site, x_km = dist_km[,1], y_km = dist_km[,2])

morisita.raw.clean <- morisita.raw %>% mutate(Species = sub(" agg.", "", Species)) %>% #harmonizing names
  select(-c("Lasioglossum pauxillum", "Andrena minutula")) # not needed for the stat models

## calculating species mean for resource overlap
morisita.raw.mean <- morisita.raw.clean %>% group_by(Species, Year) %>% summarise(Morisita.mean.hb = mean(`Apis mellifera`, na.rm = T),
                                                                            Morisita.mean.bl = mean(`Bombus lapidarius`, na.rm = T))
morisita.raw.mean2 <- morisita.raw.clean %>% group_by(Species) %>% summarise(Morisita.mean.hb2 = mean(`Apis mellifera`, na.rm = T),
                                                                       Morisita.mean.bl2 = mean(`Bombus lapidarius`, na.rm = T))

## CALCULATING PATHOGEN EXPOSURE TO THE KEY HOSTS (A. MELLIFERA AND B. LAPIDARIUS)
hb.exposure <- data.pathogen.both %>% filter(Species == 'Apis mellifera') %>% 
  mutate(across(DWVB.abs:ABPV.abs, function(x) ifelse(x == 0, NA, x))) %>% 
  group_by(Site, Year) %>% 
  summarise(dwvb.p = mean(dwvb), dwvb.q.p = mean(DWVB.abs, na.rm = T),
            bqcv.p = mean(bqcv), bqcv.q.p = mean(BQCV.abs, na.rm = T),
            abpv.p = mean(abpv), abpv.q.p = mean(ABPV.abs, na.rm = T)) %>%
  mutate(across(everything(), function(x) ifelse(is.na(x), 0, x))) %>% 
  left_join(abundance.both %>% filter(Species == 'Apis mellifera'), by = c('Year','Site')) %>% 
  rowwise() %>%
  left_join(flower.both %>% select(Site, Year, flower_dens), by = c('Site', 'Year')) %>% 
  mutate(dwvb.f = log(dwvb.p * dwvb.q.p * (bee_abundance/flower_dens/10000 +1) + 1), # from m2 to km2
         bqcv.f = log(bqcv.p * bqcv.q.p *  (bee_abundance/flower_dens/10000 +1) + 1),
         abpv.f = log(abpv.p * abpv.q.p *  (bee_abundance/flower_dens/10000 +1) + 1),
         hb_dens = bee_abundance/flower_dens/10000 +1) %>% 
  select(Site, Year, dwvb.f, bqcv.f, hb_dens)


bl.exposure <- data.pathogen.both %>% filter(Species == 'Bombus lapidarius') %>% 
  mutate(across(DWVB.abs:ABPV.abs, function(x) ifelse(x == 0, NA, x))) %>%
  # adding the two missing sites filled with NAs
  right_join(hb.exposure %>% select(Site, Year), by = c("Site", "Year")) %>%
  group_by(Site, Year) %>% 
  summarise(dwvb.p = mean(dwvb), dwvb.q.p = mean(DWVB.abs, na.rm = T),
            bqcv.p = mean(bqcv), bqcv.q.p = mean(BQCV.abs, na.rm = T),
            abpv.p = mean(abpv), abpv.q.p = mean(ABPV.abs, na.rm = T), .groups = "drop") %>%
  mutate(dwvb.q.p = ifelse(dwvb.p == 0, 0, dwvb.q.p),
         bqcv.q.p = ifelse(bqcv.p == 0, 0, bqcv.q.p),
         abpv.q.p = ifelse(abpv.p == 0, 0, abpv.q.p)) %>% 
  # calculating mean prevalence and load for the two sites where B. lapidarius was not screened
  mutate(across(ends_with(".p"), ~ifelse(is.na(.x), mean(.x, na.rm = T), .x))) %>%
  left_join(flower.both %>% select(Site, Year, flower_dens), by = c('Site', 'Year')) %>% 
  left_join(abundance.both %>% filter(Species == 'Bombus lapidarius'), by = c('Year','Site')) %>% 
  rowwise() %>%
  mutate(dwvb.bl.f = log(dwvb.p *  (bee_abundance/flower_dens/10000 +1) * dwvb.q.p + 1),
         bqcv.bl.f = log(bqcv.p *  (bee_abundance/flower_dens/10000 +1) * bqcv.q.p + 1),
         abpv.bl.f = log(abpv.p *  (bee_abundance/flower_dens/10000 +1) * abpv.q.p + 1)) %>% 
  select(Site, Year, abpv.bl.f)


### COMBINING THE DATA
abundance.site = abundance.both %>% 
  group_by(Site, Year) %>% 
  summarise(total_bee_abundance = sum(bee_abundance))

data.both <- data.pathogen.both %>% 
  left_join(dens, by = 'Site') %>% 
  mutate(Density = ifelse(Year == 2021, 0, Density),
         Year = as.numeric(Year)) %>%
  left_join(traits, by = "Species") %>%
  left_join(flower.both %>% select(Site, Year, flower_dens), by = c('Year', 'Site')) %>%
  left_join(hb.exposure, by = c('Site', 'Year')) %>%
  left_join(bl.exposure, by = c('Site', 'Year')) %>%
  left_join(abundance.site, by = c('Site', 'Year')) %>%
  mutate(across(DWVB.abs:ABPV.abs, function(x) x/BUFFER)) %>%
  mutate(total_bee_dens = log(total_bee_abundance/flower_dens/10000 + 1)) %>%
  left_join(network.metrics.both, by = c('Site', 'Year')) %>%
  left_join(morisita.raw %>% 
              rename(Morisita.hb = `Apis mellifera`, 
                     Morisita.bl = `Bombus lapidarius`) %>%
              mutate(Species = sub(" agg.", "", Species)), by = join_by("Site", "Year", "Species")) %>%
  left_join(morisita.raw.mean, by = join_by( "Year", "Species")) %>%
  left_join(morisita.raw.mean2, by = join_by("Species")) %>%
  mutate(Morisita.hb = coalesce(Morisita.hb, Morisita.mean.hb, Morisita.mean.hb2),
         Morisita.bl = coalesce(Morisita.bl, Morisita.mean.bl, Morisita.mean.bl2)) %>%
  select(-c(Morisita.mean.bl2, Morisita.mean.bl, Morisita.mean.hb, Morisita.mean.hb2)) %>%
  left_join(coord2, by = 'Site') %>%
  mutate(Year = as.factor(Year), Density = as.factor(Density)) %>%
  #removing bees with no interaction recorded
  filter(!(is.na(Morisita.hb) & is.na(Morisita.bl)))
  

#save(data.both, file = 'Data/final_dataframe.RData')


####  SPLITTING DATA INTO BEE GROUPS, SCALE AND FILTER THE OUTLIER
data.hb <- data.both %>% filter(Species == 'Apis mellifera') %>% 
  ungroup() %>%
  mutate(across(flower_dens:Morisita.bl, ~scale(.)[,1], .names = "{.col}.s"))

data.bb <- data.both %>% filter(Group == 'bb') %>% 
  ungroup() %>%
  mutate(across(flower_dens:Morisita.bl, ~scale(.)[,1], .names = "{.col}.s"))

data.bb.nolp <- data.both %>% filter(Group == 'bb') %>% filter(Species != 'Bombus lapidarius') %>% 
  ungroup() %>%
  mutate(across(flower_dens:Morisita.bl, ~scale(.)[,1], .names = "{.col}.s"))

data.wb <- data.both %>% filter(Group == 'wb') %>% ungroup() %>%
  mutate(across(flower_dens:Morisita.bl, ~scale(.)[,1], .names = "{.col}.s"))


data.site <- data.both %>% ungroup() %>% 
  distinct(Site, Year, Density, Connectance, flower_dens, hb_dens, bee_rich, plant_rich, total_bee_abundance, total_bee_dens, sum_nodes, x_km, y_km) %>%
  mutate(flower_dens.s = scale(flower_dens)[,1], total_bee_dens.s = scale(total_bee_dens)[,1])

data.morisita.bl <- data.both %>% ungroup() %>% 
  distinct(Site, Year, Density, Species, Morisita.bl, flower_dens, total_bee_abundance, total_bee_dens, plant_rich, bee_rich, sum_nodes, x_km, y_km) %>% 
  filter(Species != 'Bombus lapidarius') %>%
  mutate(across(Morisita.bl:total_bee_dens, ~scale(.)[,1], .names = "{.col}.s")) %>%
  # skewing the 0s and 1s by a marginal number to fit the Beta distribution in brms
  mutate(Morisita.bl.ns = case_when(
    Morisita.bl == 0 ~ 0 + 1e-4,
    Morisita.bl == 1 ~ 1 - 1e-4,
    TRUE ~ Morisita.bl),
    Morisita.bl.l.s = scale(log(Morisita.bl+1e-4))[,1])

data.morisita.hb <- data.both %>% ungroup() %>% distinct(Site, Year, Density, Species, Morisita.hb,flower_dens, total_bee_abundance, total_bee_dens, plant_rich, bee_rich, sum_nodes, x_km, y_km) %>% filter(Species != 'Apis mellifera') %>%
  mutate(across(Morisita.hb:total_bee_dens, ~scale(.)[,1], .names = "{.col}.s")) %>%
  # skewing the 0s and 1s by a marginal number to fit the Beta distribution in brms
  mutate(Morisita.hb.ns = case_when(
    Morisita.hb == 0 ~ 0 + 1e-4,
    Morisita.hb == 1 ~ 1 - 1e-4,
    TRUE ~ Morisita.hb),
    Morisita.hb.l.s = scale(log(Morisita.hb+1e-4))[,1])

### SOME SIMPLE RAW DATA VISUALISATION
# dotchart(data.morisita$Morisita.hb)
# dotchart(data.site$Connectance)
# dotchart(data.site$abpv.bl.f)
