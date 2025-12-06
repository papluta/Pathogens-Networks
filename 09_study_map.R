library(sf)
library(osmdata)
library(readr)
library(dplyr)
library(tidyverse)
library(cowplot)

coord = read.csv('Data/coordinates.csv')
dens = read.csv('Data/density.csv')
subset <- c('Goe1392', 'Goe1425', 'Goe235', 'Goe288', 'Goe47', 'Goe595',
            'Gos1', 'Gos2', 'Nor1', 'Nor1070', 'Nor1145', 'Nor264', 'Nor508', 'Nor918', 'WM1249', 'WM630')

coord2 = coord %>% mutate(Subset = ifelse(Site %in% subset, TRUE, FALSE)) %>% 
  left_join(dens, by = 'Site')

# get the box
lat.max = max(coord$Lat) + 0.3
lat.min = min(coord$Lat) - 0.3
lon.max = max(coord$Long) + 0.3
lon.min = min(coord$Long) - 0.3

bb <- matrix(c(lon.min, lat.min, lon.max, lat.max), nrow = 2, ncol = 2, dimnames = list(c("x","y"), c("min","max")))

# takes a long time, might time out
osm.map <- bb %>%
  opq() %>%
  add_osm_feature(key = "boundary", value = c("administrative")) %>%
  osmdata_sf()


map = ggplot() +
  geom_sf(data = osm.map$osm_multipolygons %>% 
            filter(border_type == 'county' & name %in% c('Landkreis Göttingen','Landkreis Goslar',
                                                         'Landkreis Northeim')), col = 'black', fill = 'gray',alpha = 0.4)+
  geom_sf(data = osm.map$osm_multipolygons %>% 
            filter(name %in% c('Werra-Meißner-Kreis')), col = 'black', fill = 'gray', alpha = 0.4)+
  theme_nothing()+
  labs(x = NULL, y = NULL)+
  geom_point(data = coord2, aes(x = Long, y = Lat, fill = as.factor(Density), shape = Subset), size = 4)+
  scale_fill_manual(values = c('#ddcd4b','#4b5bdd'))+
  scale_shape_manual(values = c(21,22))+
  coord_sf(xlim = c(lon.min, lon.max), ylim = c(lat.min, lat.max), expand = FALSE)+
  geom_text(aes(x = 9.9165, y = 51.54131), label = 'Göttingen')

# ggsave(file = 'Data/fig/map_sites.pdf', width = 8, height = 8)
