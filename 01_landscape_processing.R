library(dplyr)
library(readr)
library(tidyverse)

#radii.raw <- list.files(path = 'Data/raw/', pattern = '240116', full.names = T) %>% 
# set_names() %>%
#  map(~read_csv(.))

r500 <- read_csv('Data/Landuse_Combee2022_500m_updated.csv')
r1000 <- read_csv('Data/Landuse_Combee2022_1km_updated.csv')
radii.raw <- list(`500m` = r500, `1000m` = r1000)

qgis_fun <- function(x) { x %>%
    mutate(area_ha = area1/10000) %>% 
    mutate(Landscape3 = sub('Wm', 'WM', Landscape3)) %>%
    group_by(Landscape3,Combi9) %>%
    summarise(sum=sum(area_ha)) %>%
    drop_na(Landscape3) %>%
    pivot_wider(names_from = Combi9, values_from = sum) %>%
    replace(is.na(.), 0) %>%
    rename(Site = Landscape3) %>%
    mutate(semi_natur = semi_natur + Grassy_str) %>%
    select(-Grassy_str) %>%
    rename(Grassy_str = Fieldedge)
}

radii <- lapply(radii.raw, qgis_fun)
#n <- list.files(path = 'landuse_recent/raw/', pattern = '240116')
#nn <- sub(".csv","m", n)
#nnn <- sub("240116_Fareaclipbuffer","", nn)

#names(radii) <- nnn

#save(radii, file = "Data/radii.RData")


land_fun <- function(x){ x %>% 
    mutate(SNH = semi_natur + Fallow + Other_AUM + Flower_fieBS2 + Flower_fieBS12/2 + Grassy_str,
           Ann.fl = Flower_fieBS11 + Flower_fieBS12/2,
           SNH.nf = semi_natur + Fallow + Other_AUM,
           Total.fl = Flower_fieBS2 + Flower_fieBS12 + Flower_fieBS11,
           #Urban.green = Allotments + Park + Grass + Orchard,
           Grassland = Grassland + GrasslandBV1,
           AES = Other_AUM + Flower_fieBS2 + Flower_fieBS12 + Flower_fieBS11 + CropBV1) %>%
    rename(Org.farm = CropBV1) %>% 
    dplyr::select(Site, Ann.fl, SNH, Org.farm, SNH.nf, Total.fl, Grassland, Forest, Urban, Crop, AES)
}

#land.all <- radii %>% lapply(land_fun) %>% 
#  bind_rows(.id = "radius")

#save(land.all, file = 'Data/240617_landuse2022.RData')
