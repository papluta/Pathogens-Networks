library(dplyr)
library(readr)
library(tidyverse)

radii.raw <- list.files(path = 'Data/raw/', pattern = '240116', full.names = T) %>% 
  set_names() %>%
  map(~read_csv(.))


qgis_fun <- function(x) { x %>%
    mutate(area_ha = area_m/10000) %>% 
    group_by(Name,Combee9) %>%
    summarise(sum=sum(area_ha)) %>%
    drop_na(Name) %>%
    pivot_wider(names_from = Combee9, values_from = sum) %>%
    replace(is.na(.), 0)%>%
    rename(Site = Name)
}

radii <- lapply(radii.raw, qgis_fun)
n <- list.files(path = 'landuse_recent/raw/', pattern = '240116')
nn <- sub(".csv","m", n)
nnn <- sub("240116_Fareaclipbuffer","", nn)

names(radii) <- nnn

#save(radii, file = "Data/radii.RData")

radii[['200m']]$Grassy_str <- 0
radii[['200m']]$Allotments <- 0
radii[['300m']]$Allotments <- 0

land_fun <- function(x){ x %>% 
    mutate(SNH = semi_natur + Fallow + Other_AUM + Flower_fieBS2 + Flower_fieBS12/2 + Grassy_str,
           Ann.fl = Flower_fieBS11 + Flower_fieBS12/2,
           SNH.nf = semi_natur + Fallow + Other_AUM,
           Total.fl = Flower_fieBS2 + Flower_fieBS12 + Flower_fieBS11,
           Urban.green = Allotments + Park + Grass + Orchard,
           Grassland = Grassland + GrasslandBV1,
           Other = Other + Plantation) %>%
    rename(Org.farm = CropBV1) %>% 
    dplyr::select(Site, Ann.fl, SNH, Org.farm, SNH.nf, Total.fl, Urban.green, Grassland, Forest, Scrub, Urban, Crop, Other)
}

land.all <- radii %>% lapply(land_fun) %>% 
  bind_rows(.id = "radius")