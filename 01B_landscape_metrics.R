library(terra)
library(landscapemetrics)
library(dplyr)
library(readr)
library(tidyverse)
library(sf)


rastlist <- list.files(path = "~/rasters shannon/", pattern='rast500m', 
                       all.files=TRUE, full.names=TRUE)
allrasters <- lapply(rastlist, rast)
names <- list.files(path = "~/rasters shannon/", pattern='rast500m', 
                    all.files=TRUE, full.names=F)

names(allrasters) <- sub("rast500m_Landscape3_","", names) 
names(allrasters) <- sub(".tif","", names(allrasters)) 

## metrics
shannon1000 <- lapply(allrasters, lsm_l_shdi) %>% bind_rows(.id = 'Site') %>% dplyr::select(Site, SH = value)
shannon500 <- lapply(allrasters, lsm_l_shdi) %>% bind_rows(.id = 'Site') %>% dplyr::select(Site, SH = value)


iji1000 <- lapply(allrasters, lsm_l_iji) %>% bind_rows(.id = 'Site') %>% dplyr::select(Site, Iji = value) #%>% select(-Site) %>% cbind(Site) 
iji500 <- lapply(allrasters, lsm_l_iji) %>% bind_rows(.id = 'Site') %>% dplyr::select(Site, Iji = value)


#class.ennd <- lapply(allrasters, lsm_c_enn_mn) %>% bind_rows(.id = 'Site') %>% rename(Ennd = value) %>% select(Site, class, Ennd)
#Ennd.flowers1000 <- class.ennd %>% filter(class == 6) %>% rename(Ennd.flowers = Ennd)
#Ennd.org1000 <- class.ennd %>% filter(class == 4) %>% rename(Ennd.organic = Ennd) 

edge.dens1000 <- lapply(allrasters, lsm_l_ed) %>% bind_rows(.id = 'Site') %>% dplyr::select(Site, Edge.dens = value)
edge.dens500 <- lapply(allrasters, lsm_l_ed) %>% bind_rows(.id = 'Site') %>% dplyr::select(Site, Edge.dens = value)

land_metrics1000 <- shannon1000 %>% left_join(iji1000, by = 'Site') %>% left_join(edge.dens1000, by = 'Site')
land_metrics500 <- shannon500 %>% left_join(iji500, by = 'Site') %>% left_join(edge.dens500, by = 'Site')

save(land_metrics500, land_metrics1000, file = 'Data/240617landscape_metrics2022.RData')


## for flower aggregation ##
rastlist <- list.files(path = "~/rasters flowers/", pattern='rast500m', 
                       all.files=TRUE, full.names=TRUE)
allrasters <- lapply(rastlist, rast)
names <- list.files(path = "~/rasters flowers/", pattern='rast500m', 
                    all.files=TRUE, full.names=F)

names(allrasters) <- sub("rast500m_Landscape3_","", names) 
names(allrasters) <- sub(".tif","", names(allrasters)) 

## metrics
class.ennd <- lapply(allrasters, lsm_c_enn_mn) %>% bind_rows(.id = 'Site') %>% rename(Ennd = value) %>% select(Site, class, Ennd)
Ennd.flowers500 <- class.ennd %>% filter(class == 6) %>% rename(Ennd.flowers = Ennd)
Ennd.org500 <- class.ennd %>% filter(class == 4) %>% rename(Ennd.organic = Ennd) 

