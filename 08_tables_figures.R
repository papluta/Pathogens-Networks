source('05_main_data_file.R') 

library(openxlsx2)
library(cowplot)
library(bayestestR)
library(tidybayes)
library(viridis)
library(rphylopic)
library(brms)
library(tidybayes)
library(rstan)
library(bayesplot)

### get phylopics of the bee groups

img_AM <- get_phylopic(uuid = get_uuid(name = 'Apis mellifera', n = 1))
img_BT <- get_phylopic(uuid = get_uuid(name = 'Bombus', n = 2)[2])
img_LS <- get_phylopic(uuid = get_uuid(name = 'Lasioglossum sordidum', n = 1))

###--------------------------------------------------------------------------###

#################
###   TABLES  ###
#################

# create each table and then combine into the supplementary_tables.xslx 

# Table S1

t1 <- list(hind.2021.2, hind.2022.2, 
           pind.2021.2, pind.2022.2, 
           tind.2021.2, 
           lind.2021.2, lind.2022.2,
           wind.2021.2, wind.2022.2) %>% 
  map(., ~ .x %>% select(Site, Year, Species)) %>%
  bind_rows(.id = "Group") %>%
  mutate(Group = case_when(Species == "Apis mellifera" ~"hb",
                           grepl("Bombus", Species) ~ "bb",
                           TRUE ~ "wb")) %>%
  count(Group, Site, Year) %>%
  pivot_wider(names_from = c(Year, Group), values_from = n) %>%
  left_join(coord, by = "Site") %>%
  select(Site, Long, Lat, starts_with("2021"), starts_with("2022"))

#write.csv(t1, file = paste0("Data/Results/", date, "_t1.csv"), row.names = F)

# Table S2

t2 <- data.pathogen.both  %>% 
  mutate(across(ends_with('.abs'), ~ ifelse(.x == 0, NA, .x))) %>%
  group_by(Species) %>%
  summarise(n = n(), 
            across(all_of(c('dwvb', 'bqcv', 'abpv')), ~ mean(.x, na.rm = T), .names = "{.col}.p"),
            across(all_of(c('DWVB.abs', 'BQCV.abs', 'ABPV.abs')), ~ mean(.x, na.rm = T), .names = "{.col}.p"),
            across(all_of(c('DWVB.abs', 'BQCV.abs', 'ABPV.abs')), ~ sd(.x, na.rm = T), .names = "{.col}.sd")
  ) %>%
  mutate(across(ends_with('.p'), ~ ifelse(is.na(.x), 0, .x))) %>%
  left_join(data.pathogen.both  %>% distinct(Site, Species) %>% count(Species) %>%
              rename(n_site = n), by = "Species") %>%
  rename_with(., tolower) %>%
  select(species, n, n_site, starts_with("dwvb"), starts_with("bqcv"), starts_with("abpv"))

# write.csv(t2, file = paste0("Data/Results/", date, "_t2.csv"), row.names = F)

# Table S3 - not calculable in R (primers)

# Table S4

t4 <- r0.results.sensitivity %>% bind_rows(.id = "Virus") %>%
  mutate(Sociality = factor(Sociality, levels = c("1", "2", "5", "10"))) %>%
  group_by(Sociality, Virus, Species) %>%
  summarise(r0.m = mean(r0), r0.sd = sd(r0)) %>%
  transmute(Sociality,
            Virus = toupper(Virus),
            Species,
            R0 = paste0(round(r0.m,2), " ± ", round(r0.sd,2))
            ) %>%
  pivot_wider(names_from = c(Sociality, Virus), values_from = R0)

# write.csv(t4, file = paste0("Data/Results/", date, "_t4.csv"), row.names = F)

# Table S5

t5 <- bind_rows(true.prev.2021 %>% mutate(Year = 2021), true.prev.2022 %>% mutate(Year = 2022)) %>% 
  distinct(Site, Species, Year, .keep_all = T) %>%
  group_by(Species) %>%
  summarise(across(ends_with("true"), ~ mean(.x), .names = "{.col}.m"),
            across(ends_with("true"), ~ sd(.x), .names = "{.col}.sd"),
            n = sum(n)) %>%
  left_join(t2 %>% select(Species = species, dwvb.p, bqcv.p, abpv.p), by = "Species") %>%
  select(Species, n, starts_with("dwvb"), starts_with("bqcv"), starts_with("abpv"))

# write.csv(t5, file = paste0("Data/Results/", date, "_t4.csv"), row.names = F)

# Table S6

r0.comb <- left_join(r0.sim.prev.dwvb, r0.sim.prev.bqcv, by = "Species", suffix = c(".dwvb", ".bqcv")) %>%
  left_join(r0.sim.prev.abpv %>% rename_with(~ paste0(.x, ".abpv"), .cols = 2:last_col()), by = "Species")

t6 <- r0.comb %>%
  select(Species, starts_with("n_n"), starts_with("r0")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  # only one number of networks because it's the same across viruses
  select(Species, n_networks = n_networks.dwvb, starts_with("r0"))

# write.csv(t6, file = paste0('Data/Results/', date, '_t6.csv'))

# Table S7

t7 <- r0.comb %>% select(-starts_with("n_n"), -starts_with("r0"))

# write.csv(t7, file = paste0('Data/Results/', date, '_t7.csv'))

# Table S8

t8 <- bf.list %>%
  map(~ .x[["bf"]]) %>% rbind() %>% as.data.frame() %>%
  pivot_longer(everything()) %>%
  mutate(model1 = paste0(sub(".*\\.", "", name), " removal"),
         model2 = paste0(sub(".*\\.", "", name), " intercept-only"),
         virus = sub("...$","",sub("...","",name))) %>%
  select(virus, model1, model2, `BF in support of Model 1` = value)

# write.csv(t8, file = paste0('Data/Results/', date, '_t8.csv'))


# Table S9
t9 <- lapply(WB_models_combined_raw, custom_summary) %>% bind_rows(.id = 'id') 

preview = t9 %>% filter(pd > 0.94) %>% 
  filter(Predictor != 'Intercept' & Predictor != 'hu_Intercept'& Predictor != 'SpeciesBombuspascuorum' & 
           Predictor != 'hu_SpeciesBombuspascuorum' & Predictor != 'SpeciesBombusterrestris' & Predictor != 'hu_SpeciesBombusterrestris') %>% 
  arrange(Response, Predictor) %>% rownames_to_column() %>% select(-rowname) %>% filter(Predictor != 'Year2022' & Predictor != 'hu_Year2022')

t9 <- t9 %>% mutate(term = ifelse(grepl('hu_', Predictor), 'Probability of presence', 'Viral load')) %>% 
  mutate(Predictor = sub('hu_', '', Predictor), Bee_group = ifelse(grepl('.h$', id), 'Honey bee', 'Bumble bee')) %>% 
  mutate(Bee_group = ifelse(grepl('(.w)$', id), 'Other wild bee', Bee_group)) %>% rename(Virus = Response) %>%
  mutate(Predictor = recode_factor(as.factor(Predictor), 'dwvb.f.s' = 'Viral exposure', 'bqcv.f.s' = 'Viral exposure','bqcv.f.s:Morisita.hb.s' = 'Viral exposure : HB resource overlap', 'abpv.bl.f.s' = 'Viral exposure', 
                                   'abpv.bl.f.s:Morisita.bl.s' = 'Viral exposure : BL resource overlap','Bee_rich' = 'Bee richness', 'Year2022' = 'Year 2022', 'Morisita.hb.s' = 'HB resource overlap', 'Morisita.bl.s' = 'BL resource overlap',
                                   'SpeciesBombuspascuorum' = 'B. pascuorum', 'SpeciesBombusterrestris' = 'B. terrestris', "Connectance.s" = "Connectance", "total_bee_dens.s" = "Total bee density",
                                   "sum_nodes" = "Sum of nodes"),
         Virus = recode_factor(as.factor(Virus), 'DWVBabs' = 'DWV-B', 'BQCVabs' = 'BQCV', 'ABPVabs' = 'ABPV', 'SBVabs' = 'SBV')) %>% select(-id) %>% 
  mutate(Virus = factor(Virus, levels = c('DWV-B', 'BQCV', 'ABPV', 'SBV')), Bee_group = factor(Bee_group, levels = c('Honey bee', 'Bumble bee', 'Other wild bee')),
         Predictor = factor(Predictor, levels = c('Intercept', 'Viral exposure : HB resource overlap', 'Viral exposure : BL resource overlap', 'Viral exposure', 'HB resource overlap', 'BL resource overlap', 'Total bee density','Connectance','Bee richness', 
                                                  'Flower cover', "Sum of nodes",'Year 2022', 'B. pascuorum', 'B. terrestris'))) %>%
  relocate(Virus, Bee_group, term, Predictor) %>% arrange(Virus, Bee_group, term, Predictor)


# write.csv(t9, file = paste0('Data/Results/', date, '_t9.csv'), row.names = F)

# Table S10
t10 <- lapply(WB_models_combined_raw_size, custom_summary) %>% bind_rows(.id = 'id') 

t10 <- t10 %>% mutate(term = ifelse(grepl('hu_', Predictor), 'Probability of presence', 'Viral load')) %>% 
  mutate(Predictor = sub('hu_', '', Predictor), Bee_group = ifelse(grepl('.h$', id), 'Honey bee', 'Bumble bee')) %>% 
  mutate(Bee_group = ifelse(grepl('(.w)$', id), 'Other wild bee', Bee_group)) %>% rename(Virus = Response) %>%
  mutate(Predictor = recode_factor(as.factor(Predictor), 'dwvb.f.s' = 'Viral exposure', 'bqcv.f.s' = 'Viral exposure','bqcv.f.s:Morisita.hb.s' = 'Viral exposure : HB resource overlap', 'abpv.bl.f.s' = 'Viral exposure', 
                                   'abpv.bl.f.s:Morisita.bl.s' = 'Viral exposure : BL resource overlap','Bee_rich' = 'Bee richness', 'Year2022' = 'Year 2022', 'Morisita.hb.s' = 'HB resource overlap', 'Morisita.bl.s' = 'BL resource overlap',
                                   'SpeciesBombuspascuorum' = 'B. pascuorum', 'SpeciesBombusterrestris' = 'B. terrestris', "Connectance.s" = "Connectance", "total_bee_dens.s" = "Total bee density",
                                   "sum_nodes" = "Sum of nodes"),
         Virus = recode_factor(as.factor(Virus), 'DWVBabs' = 'DWV-B', 'BQCVabs' = 'BQCV', 'ABPVabs' = 'ABPV', 'SBVabs' = 'SBV')) %>% select(-id) %>% 
  mutate(Virus = factor(Virus, levels = c('DWV-B', 'BQCV', 'ABPV', 'SBV')), Bee_group = factor(Bee_group, levels = c('Honey bee', 'Bumble bee', 'Other wild bee')),
         Predictor = factor(Predictor, levels = c('Intercept', 'Viral exposure : HB resource overlap', 'Viral exposure : BL resource overlap', 'Viral exposure', 'HB resource overlap', 'BL resource overlap', 'Total bee density','Connectance','Bee richness', 
                                                  'Flower cover', "Sum of nodes",'Year 2022', 'B. pascuorum', 'B. terrestris'))) %>%
  relocate(Virus, Bee_group, term, Predictor) %>% arrange(Virus, Bee_group, term, Predictor)


# write.csv(t10, file = paste0('Data/Results/', date, '_t10.csv'), row.names = F)

# Table S11

t11 <- lapply(site_models_500, custom_summary) %>% bind_rows(.id = 'mod') %>%
  mutate(Predictor = recode_factor(as.factor(Predictor), 'FL.sum.s' = 'Flower density', "total_bee_dens.s" = "Bee density",
                                   "sum_nodes" = "Sum of nodes",
                                   'FL_per_agr:HB_per_agr' = 'HB density:Flower cover',  'Bee_rich' = 'Bee richness'),
         Response = recode_factor(as.factor(Response), 'logsumab' = 'Bee density', 'Morisitahbns' = 'HB niche overlap',
                                  'Morisitablns' = 'BL niche overlap', 'Beerich' = 'Bee richness')) %>% select(-mod)

# write.csv(t11, file = paste0('Data/Results/', date, '_t11.csv'), row.names = F)


# Table S12

t12 <- lapply(WB_models_combined_raw_noimputedint, custom_summary) %>% bind_rows(.id = 'id') 

t12 <- t12 %>% mutate(term = ifelse(grepl('hu_', Predictor), 'Probability of presence', 'Viral load')) %>% 
  mutate(Predictor = sub('hu_', '', Predictor), Bee_group = ifelse(grepl('.h$', id), 'Honey bee', 'Bumble bee')) %>% 
  mutate(Bee_group = ifelse(grepl('(.w)$', id), 'Other wild bee', Bee_group)) %>% rename(Virus = Response) %>%
  mutate(Predictor = recode_factor(as.factor(Predictor), 'dwvb.f.s' = 'Viral exposure', 'bqcv.f.s' = 'Viral exposure','bqcv.f.s:Morisita.hb.s' = 'Viral exposure : HB resource overlap', 'abpv.bl.f.s' = 'Viral exposure', 
                                   'abpv.bl.f.s:Morisita.bl.s' = 'Viral exposure : BL resource overlap','Bee_rich' = 'Bee richness', 'Year2022' = 'Year 2022', 'Morisita.hb.s' = 'HB resource overlap', 'Morisita.bl.s' = 'BL resource overlap',
                                   'SpeciesBombuspascuorum' = 'B. pascuorum', 'SpeciesBombusterrestris' = 'B. terrestris', "Connectance.s" = "Connectance", "total_bee_dens.s" = "Total bee density",
                                   "sum_nodes" = "Sum of nodes"),
         Virus = recode_factor(as.factor(Virus), 'DWVBabs' = 'DWV-B', 'BQCVabs' = 'BQCV', 'ABPVabs' = 'ABPV', 'SBVabs' = 'SBV')) %>% select(-id) %>% 
  mutate(Virus = factor(Virus, levels = c('DWV-B', 'BQCV', 'ABPV', 'SBV')), Bee_group = factor(Bee_group, levels = c('Honey bee', 'Bumble bee', 'Other wild bee')),
         Predictor = factor(Predictor, levels = c('Intercept', 'Viral exposure : HB resource overlap', 'Viral exposure : BL resource overlap', 'Viral exposure', 'HB resource overlap', 'BL resource overlap', 'Total bee density','Connectance','Bee richness', 
                                                  'Flower cover', "Sum of nodes",'Year 2022', 'B. pascuorum', 'B. terrestris'))) %>%
  relocate(Virus, Bee_group, term, Predictor) %>% arrange(Virus, Bee_group, term, Predictor)

# combining tables

# results of brms models might differ slightly between runs


tables_list <- list(
  "Content" = data.frame(), # list of content added manually
  "1. Location" = t1,
  "2. No. of ind. viral screen" = t2,
  "3. Primer info" = data.frame(), # molecular primers
  "4. Sensitivity R0" = t4,
  "5. Adjusted prevalence" = t5,
  "6. R0 results" = t6,
  "7. Prevalence simulation" = t7,
  "8. Bayes factor" = t8,
  "9. Hurdle models" = t9,
  "10. Sensitivity hurdle size" = t10,
  "11. Landscape models" = t11,
  "12. Sensitivity hurdle overlap" = t12
)

write_xlsx(tables_list, file = paste0("Data/Results/", date, "_supplementary_tables.xlsx"))

###--------------------------------------------------------------------------###

#######################
####### FIGURES #######
#######################

## first create individual plots, and then arrange them in the figures like in the ms

## PREVALENCE AND LOAD

p = data.pathogen.both  %>% 
  pivot_longer(cols = dwvb:abpv, names_to = 'virus', values_to = 'presence') %>% 
  filter(Species %in% c(unique(data.wb$Species), unique(data.bb$Species),'Apis mellifera')) %>%
  group_by(Group, virus, Year) %>%
  summarise(prevalence = mean(presence), n = n()) %>% mutate(Group = factor(Group, levels = c('hb', 'bb', 'wb'))) %>% 
  mutate(virus = case_match(virus, 'dwvb' ~ 'DWV-B', 'bqcv' ~ 'BQCV', 'abpv' ~ 'ABPV')) %>%
  mutate(virus = factor(virus, levels = c('DWV-B', 'BQCV', 'ABPV'))) %>% 
  ggplot(aes(virus, prevalence, fill = Group)) +
  geom_bar(stat = 'identity', col = 'black', position = position_dodge())+
  scale_fill_manual(values = c('#6a7f54', '#d55e00', '#f0e442'), labels = c('Honeybee', 'Bumblebee', 'Other wild bee'))+
  facet_wrap(~Year, nrow = 1)+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color = 'black'), legend.position = 'bottom',
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+  
  labs(x = NULL, y = 'Prevalence', fill = 'Bee group', title = expression(bold('a')))


l = data.pathogen.both %>% 
  mutate(across(DWVB.abs:ABPV.abs, function(x) x/BUFFER)) %>% 
  pivot_longer(cols = DWVB.abs:ABPV.abs, names_to = 'virus', values_to = 'load') %>% 
  filter(Species %in% c(unique(data.wb$Species), unique(data.bb$Species),'Apis mellifera')) %>%
  mutate(Group = factor(Group, levels = c('hb', 'bb', 'wb'))) %>%  
  mutate(virus = factor(virus, levels = c('DWVB.abs', 'BQCV.abs', 'ABPV.abs'))) %>% 
  filter(load > 0) %>% 
  ggplot(aes(virus, log(load), fill = Group))+
  geom_boxplot(col = 'black')+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color = 'black'), legend.position = 'none',
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+
  labs(x = NULL, y = 'Viral laod per uL (ln)', fill = 'Bee group', title = expression(bold('b')))+
  scale_fill_manual(values = c('#6a7f54', '#d55e00', '#f0e442'), labels = c('Honeybee', 'Bumblebee', 'Other wild bee'))+
  scale_x_discrete(labels = c('DWV-B', 'BQCV', 'ABPV'))+
  facet_wrap(~Year)

# create legend
leg.p = get_plot_component(p, "guide-box-bottom", return_all = T)

prevload.raw = plot_grid(p + theme(legend.position = 'none'), l,leg.p, rel_heights = c(0.4,0.4,0.2), nrow = 3)

# small prevalence plots for Fig. 4
dw = data.both  %>% 
  pivot_longer(cols = dwvb:abpv, names_to = 'virus', values_to = 'presence') %>% 
  group_by(Group, virus) %>%
  summarise(prevalence = mean(presence)) %>% mutate(Group = factor(Group, levels = c('hb', 'bb', 'wb'))) %>% 
  mutate(virus = factor(virus, levels = c('bqcv', 'dwvb', 'abpv'))) %>% filter(virus == 'dwvb') %>%
  ggplot(aes(virus, prevalence, fill = Group)) +
  geom_bar(stat = 'identity', col = 'black', position = position_dodge(), width = 1)+
  scale_fill_manual(values = c('#6a7f54', '#d55e00', '#f0e442'), labels = c('Honeybee', 'Bumblebee', 'Other wild bee'))+
  theme_bw(base_size = 16)+
  geom_hline(aes(yintercept = mean(prevalence)), linetype = 'dashed')+
  ylim(0,1)+
  labs(x = NULL, y = NULL, fill = NULL, title = expression(bold('d')))+
  theme(axis.text = element_text(color = 'black'), legend.position = 'left', panel.grid = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+
  scale_x_discrete(labels = c('DWV-B'))+
  add_phylopic(img = img_AM, x = 0.65, y = 0.9, height = 0.15)+
  add_phylopic(img = img_BT, x = 1, y = 0.9, height = 0.15)+
  add_phylopic(img = img_LS, x = 1.35, y = 0.9, height = 0.15)

bq = data.both  %>% 
  pivot_longer(cols = dwvb:abpv, names_to = 'virus', values_to = 'presence') %>% 
  group_by(Group, virus) %>%
  summarise(prevalence = mean(presence)) %>% mutate(Group = factor(Group, levels = c('hb', 'bb', 'wb'))) %>% 
  mutate(virus = factor(virus, levels = c('bqcv', 'dwvb', 'abpv'))) %>% filter(virus == 'bqcv') %>%
  ggplot(aes(virus, prevalence, fill = Group)) +
  geom_bar(stat = 'identity', col = 'black', position = position_dodge(), width = 1)+
  scale_fill_manual(values = c('#6a7f54', '#d55e00', '#f0e442'), labels = c('Honeybee', 'Bumblebee', 'Other wild bee'))+
  geom_hline(aes(yintercept = mean(prevalence)), linetype = 'dashed')+
  theme_bw(base_size = 16)+
  ylim(0,1)+
  labs(x = NULL, y = NULL, fill = NULL, title = expression(bold('e')))+
  theme(axis.text = element_text(color = 'black'), legend.position = 'none', panel.grid = element_blank(), 
        axis.ticks.x = element_blank(),
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+
  scale_x_discrete(labels = c('BQCV'))+
  add_phylopic(img = img_AM, x = 0.65, y = 0.9, height = 0.15)+
  add_phylopic(img = img_BT, x = 1, y = 0.9, height = 0.15)+
  add_phylopic(img = img_LS, x = 1.35, y = 0.9, height = 0.15)

leg.p = get_legend(dw)
leg.p = get_plot_component(dw, 'guide-box-left', return_all = T)
prev.sep = plot_grid(dw + theme(legend.position = 'none'),leg.p, bq, rel_heights = c(0.4,0.2,0.4), nrow = 3)
prev.sep = plot_grid(dw + theme(legend.position = 'none'), bq, nrow = 1)


########################
########## R0 ##########
########################

### plotting top 4 R0 species

top4.dwvb = t6 %>% 
  arrange(desc(r0_mean.dwvb)) %>% slice_max(r0_mean.dwvb, n = 4) 

dwvb.top.r0 <- r0.results.long$dwvb %>% 
  filter(Species %in% top4.dwvb$Species) %>%
  group_by(Species) %>%
  summarise(ymin  = min(r0),
    lower = quantile(r0, 0.25),
    middle = mean(r0, na.rm = T),      # mean instead of median
    upper = quantile(r0, 0.75),
    ymax  = max(r0)) %>%
  mutate(Species = sub(" ", "\n", Species)) %>%
  ggplot(aes(x = Species)) +
  geom_boxplot(aes(
      ymin   = ymin,
      lower  = lower,
      middle = middle,
      upper  = upper,
      ymax   = ymax,
      fill   = Species),
    stat = "identity", lwd = 0.8) +
  geom_abline(intercept = 1, slope = 0, 
              linetype = 21, linewidth = 1, col = 'black') +
  theme_classic(base_size = 15) +
  scale_fill_manual(values = c('darkgray', '#6a7f54', 'darkgray', 'darkgray')) +
  theme(legend.position = 'none',
    axis.text.y = element_text(color = 'black'),
    axis.text.x = element_text(color = 'black', size = 14, face = 'italic'),
    title = element_text(size = 20),
    plot.title.position = "plot",
    plot.caption.position = "plot") +
  labs(x = NULL,
    y = expression(paste('DWV-B ', R[0][","][i])),
    title = expression(bold('a'))) +
  ylim(0, 5)

top4.bqcv = t6 %>% 
  arrange(desc(r0_mean.bqcv)) %>% slice_max(r0_mean.bqcv, n = 4) 

bqcv.top.r0 <- r0.results.long$bqcv %>% 
  filter(Species %in% top4.bqcv$Species) %>%
  mutate(r0 = log10(r0)) %>%
  group_by(Species) %>%
  summarise(ymin  = min(r0),
            lower = quantile(r0, 0.25),
            middle = mean(r0, na.rm = T),      # mean instead of median
            upper = quantile(r0, 0.75),
            ymax  = max(r0)) %>%
  mutate(Species = sub(" ", "\n", Species)) %>%
  mutate(Species = factor(Species, levels = c("Bombus\nterrestris", "Apis\nmellifera", "Bombus\nlapidarius", "Bombus\npascuorum"))) %>%
  ggplot(aes(x = Species)) +
  geom_boxplot(aes(
    ymin   = ymin,
    lower  = lower,
    middle = middle,
    upper  = upper,
    ymax   = ymax,
    fill   = Species),
    stat = "identity", lwd = 0.8) +
  geom_abline(intercept = 1, slope = 0, 
              linetype = 21, linewidth = 1, col = 'black') +
  theme_classic(base_size = 15) +
  scale_fill_manual(values = c('darkgray','#6a7f54', 'darkgray',  'darkgray'))+
  theme(legend.position = 'none',
        axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', size = 14, face = 'italic'),
        title = element_text(size = 20),
        plot.title.position = "plot",
        plot.caption.position = "plot") +
  labs(x = NULL,
       y = expression(paste('BQCV ', R[0][","][i])),
       title = expression(bold('b'))) +
  scale_y_continuous(labels = c(1,10,100,1000), breaks = c(0,1,2,3))


top4.abpv = t6 %>% 
  arrange(desc(r0_mean.abpv)) %>% slice_max(r0_mean.abpv, n = 4) 

abpv.top.r0 <- r0.results.long$abpv %>% 
  filter(Species %in% top4.abpv$Species) %>%
  group_by(Species) %>%
  summarise(ymin  = min(r0),
            lower = quantile(r0, 0.25),
            middle = mean(r0, na.rm = T),      # mean instead of median
            upper = quantile(r0, 0.75),
            ymax  = max(r0)) %>%
  mutate(Species = sub(" ", "\n", Species)) %>%
  ggplot(aes(x = Species)) +
  geom_boxplot(aes(
    ymin   = ymin,
    lower  = lower,
    middle = middle,
    upper  = upper,
    ymax   = ymax,
    fill   = Species),
    stat = "identity", lwd = 0.8) +
  geom_abline(intercept = 1, slope = 0, 
              linetype = 21, linewidth = 1, col = 'black') +
  theme_classic(base_size = 15) +
  scale_fill_manual(values = c('darkgray','darkgray', '#d55e00',  'darkgray')) +
  theme(legend.position = 'none',
        axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', size = 14, face = 'italic'),
        title = element_text(size = 20),
        plot.title.position = "plot",
        plot.caption.position = "plot") +
  labs(x = NULL,
       y = expression(paste('ABPV ', R[0][","][i])),
       title = expression(bold('c'))) +
  ylim(0, 8)

r0_plots = plot_grid(dwvb.top.r0, bqcv.top.r0, abpv.top.r0, align = 'hv', nrow = 3)


## prevalence simulation plot (heatmap)

obs = t7 %>% select(Species, contains("obs_prev_mean")) %>%
  pivot_longer(cols = -Species, names_to = "Virus", values_to = "obs_prev") %>%
  mutate(Virus = sub(".*\\.", "", Virus))

sim = t7 %>% select(Species, contains(c("sim_prev"))) %>%
  select(Species, contains(c("mean"))) %>%
  pivot_longer(cols = -Species, names_to = "Type", values_to = "sim_prev") %>%
  mutate(Virus = sub(".*\\.", "", Type),
         Host = sub(".*_no", "", Type) %>% sub("_.*", "", .),
         Type = paste0(Virus, Host))

sim.logodds <- sim %>%
  left_join(obs, by = c("Species", "Virus")) %>%
  mutate(
    odds_obs = log(obs_prev / (1 - obs_prev)),
    odds_sim = log(sim_prev / (1 - sim_prev)),
    Fold.odds = odds_obs / odds_sim,
    fold_prev = sim_prev / obs_prev
  )


heat = ggplot(sim.logodds, aes(Species, Type, fill= fold_prev)) + 
  geom_tile()+
  labs(x = NULL, y = NULL, fill = 'Fold change from\nobserved prevalence')+
  theme_map()+
  scale_fill_distiller(palette = "Spectral", direction = 1)+
  scale_y_discrete(labels = c('A. minutula','B. lapidarius', 'A. mellifera','B. lapidarius', 'A. mellifera', 'L. pauxillum' ))+
  theme(axis.ticks = element_blank(), axis.text = element_text(color = 'black', face = 'italic'), axis.text.x = element_text(angle = 90, size = 12, hjust = 0.98))
  
## r0 vs resource overlap with key host (for supplement)

dwvb_r0_overlap <- r0.results.long$dwvb %>% 
  rename(Site_number = Site) %>%
  left_join(coord %>% arrange(Site) %>% 
              transmute(Site_number = row.names(.),
                        Site), by = "Site_number") %>%
  filter(main.host == "Lasioglossum pauxillum") %>%
  left_join(morisita.raw %>% 
              mutate(Species = sub(" agg.", "", Species)) %>%
                       group_by(Site, Year) %>% 
                       summarise(Morisita.lp = mean(`Lasioglossum pauxillum`, na.rm = T),
                                 Morisita.hb = mean(`Apis mellifera`, na.rm = T)), by = c("Year", "Site"))

dwvb.r0.overlap.lp <- dwvb_r0_overlap %>% filter(Species == "Lasioglossum pauxillum") %>%
  ggplot(aes(Morisita.lp, r0))+
  geom_point(size = 3)+
  theme_bw(base_size = 16)+
  labs(x = expression(paste('Resource overlap with ', italic('L. pauxillum'))),
       y = expression(paste(italic('L. pauxillum'), " DWV-B ", R[0])),
       title = expression(bold("b")))+
  theme(panel.grid = element_blank())

dwvb.r0.overlap.hb <- dwvb_r0_overlap %>% filter(Species == "Apis mellifera") %>%
  ggplot(aes(Morisita.hb, r0))+
  geom_point(size = 3)+
  theme_bw(base_size = 16)+
  labs(x = expression(paste('Resource overlap with ', italic('A. mellifera'))),
       y = expression(paste(italic('A. mellifera'),  " DWV-B ", R[0])),
       title = expression(bold("a")))+
  theme(panel.grid = element_blank())

##

bqcv_r0_overlap <- r0.results.long$bqcv %>% 
  rename(Site_number = Site) %>%
  left_join(coord %>% arrange(Site) %>% 
              transmute(Site_number = row.names(.),
                        Site), by = "Site_number") %>%
  filter(main.host == "Apis mellifera") %>%
  left_join(morisita.raw %>% 
              mutate(Species = sub(" agg.", "", Species)) %>%
              group_by(Site, Year) %>% 
              summarise(Morisita.bl = mean(`Bombus lapidarius`, na.rm = T),
                        Morisita.hb = mean(`Apis mellifera`, na.rm = T)), by = c("Year", "Site"))

bqcv.r0.overlap.bl <- bqcv_r0_overlap %>% filter(Species == "Bombus lapidarius") %>%
  ggplot(aes(Morisita.bl, r0))+
  geom_point(size = 3)+
  theme_bw(base_size = 16)+
  labs(x = expression(paste('Resource overlap with ', italic('B. lapidarius'))),
       y = expression(paste(italic('B. lapidarius'),  " BQCV ", R[0])),
       title = expression(bold("d")))+
  theme(panel.grid = element_blank())

bqcv.r0.overlap.hb <- bqcv_r0_overlap %>% filter(Species == "Apis mellifera") %>%
  ggplot(aes(Morisita.hb, r0))+
  geom_point(size = 3)+
  theme_bw(base_size = 16)+
  labs(x = expression(paste('Resource overlap with ', italic('A. mellifera'))),
       y = expression(paste(italic('A. mellifera'), " BQCV ", R[0])),
       title = expression(bold("c")))+
  theme(panel.grid = element_blank())

##

abpv_r0_overlap <- r0.results.long$abpv %>% 
  rename(Site_number = Site) %>%
  left_join(coord %>% arrange(Site) %>% 
              transmute(Site_number = row.names(.),
                        Site), by = "Site_number") %>%
  filter(main.host == "Bombus lapidarius") %>%
  left_join(morisita.raw %>% 
              mutate(Species = sub(" agg.", "", Species)) %>%
              group_by(Site, Year) %>% 
              summarise(Morisita.bl = mean(`Bombus lapidarius`, na.rm = T),
                        Morisita.am = mean(`Andrena minutula`, na.rm = T)), by = c("Year", "Site"))

abpv.r0.overlap.bl <- abpv_r0_overlap %>% filter(Species == "Bombus lapidarius") %>%
  ggplot(aes(Morisita.bl, r0))+
  geom_point(size = 3)+
  theme_bw(base_size = 16)+
  labs(x = expression(paste('Resource overlap with ', italic('B. lapidarius'))),
       y = expression(paste(italic('B. lapidarius'),  " ABPV ",R[0])),
       title = expression(bold("e")))+
  theme(panel.grid = element_blank())

abpv.r0.overlap.am <- abpv_r0_overlap %>% filter(Species == "Andrena minutula") %>%
  ggplot(aes(Morisita.am, r0))+
  geom_point(size = 3)+
  theme_bw(base_size = 16)+
  labs(x = expression(paste('Resource overlap with ', italic('A. minutula'))),
       y = expression(paste(italic('A. minutula'),  " ABPV ",R[0])),
       title = expression(bold("f")))+
  theme(panel.grid = element_blank())

r0_overlap_plot <- plot_grid(dwvb.r0.overlap.hb, dwvb.r0.overlap.lp, bqcv.r0.overlap.hb, bqcv.r0.overlap.bl,
          abpv.r0.overlap.bl, abpv.r0.overlap.am, ncol = 2, align = "hv")


########################################
########## STATISTICAL MODELS ##########
########################################

########################
###### COMMUNITY #######
########################
pd = p_direction(hb.dens.m) %>% filter(grepl("Density", Parameter))

hb.dens.colony <- ggplot(data.site, aes(x = as.factor(Density), y = hb_dens, fill = as.factor(Density))) +
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c('#ddcd4b', '#4b5bdd'))+
  labs(x = '\nHoneybee colony density', y = expression(paste('Honeybee density ', 'per ', m^2 ,'of flowers')))+
  theme(panel.grid = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none',
        axis.ticks.x = element_blank())+
  scale_x_discrete(labels = c('Low', 'Increased'))+
  annotate("text", label = paste0("pd = ", round(pd$pd, 2)), x=1.5,  y=max(data.site$hb_dens), hjust=0.5)
  
#######################
#### hurdle models ####
#######################

## VIRAL EXPOSURE
mod = WB_models_combined_raw[['dwvb.b']]
pd = p_direction(mod) %>% filter(grepl("hu_dwvb.f.s", Parameter))


new_dat = expand.grid(
  dwvb.f.s = seq(min(mod$data$dwvb.f.s), max(mod$data$dwvb.f.s), by = 0.1),
  Morisita.bl.s = 0,
  Morisita.hb.s = 0,
  Connectance.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  Site = "Goe1425",
  Species = 'Bombus lapidarius'
)
 
posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)
data.bb.raw = data.both %>% filter(Group == 'bb') %>% ungroup() 
m = mean(data.bb.raw$dwvb.f)
sd = sd(data.bb.raw$dwvb.f)
x = data.frame(z = new_dat$dwvb.f, r = new_dat$dwvb.f * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('dwvb.f.s' == 'z')) 

exp_dwvb_b = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#6a7f54') +
  labs(x = expression(paste("Viral exposure to ", italic("A. mellifera"))), y = "DWV-B prevalence", title = NULL)+
  annotate("text", label = paste0("pd = ", round(pd$pd, 2)), x=mean(range(posterior_preds2$r, na.rm = TRUE)),  y=0.92, hjust=0.5)+
  geom_jitter(data = mod$data %>% mutate(DWVB.abs = ifelse(DWVB.abs > 0, 1, 0)) %>%
                mutate(dwvb.f.r = dwvb.f.s * sd + m), aes(dwvb.f.r, DWVB.abs), shape = 21, height = 0.025, alpha = 0.4)+
  theme_bw(base_size = 14)+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))+
  add_phylopic(img = img_BT, x = min(posterior_preds2$r), y = 0.85, height = 0.2, hjust = 0.02)

##

mod = WB_models_combined_raw[['dwvb.w']]
pd = p_direction(mod) %>% filter(grepl("b_dwvb.f.s", Parameter))


new_dat = expand.grid(
  dwvb.f.s = seq(min(mod$data$dwvb.f.s), max(mod$data$dwvb.f.s), by = 0.1),
  Morisita.bl.s = 0,
  Morisita.hb.s = 0,
  Connectance.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  Site = "Goe1425",
  Species = 'Andrena flavipes'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'mu', re_formula = NA)
data.bb.raw = data.both %>% filter(Group == 'wb') %>% ungroup() 
m = mean(data.bb.raw$dwvb.f)
sd = sd(data.bb.raw$dwvb.f)
x = data.frame(z = new_dat$dwvb.f, r = new_dat$dwvb.f * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('dwvb.f.s' == 'z')) 


exp_dwvb_w = ggplot(posterior_preds2, aes(x = r, y = mu)) +
  stat_lineribbon(.width = c(0.95), fill = '#6a7f54') +
  labs(x = expression(paste("Viral exposure to ", italic("A. mellifera"))), y = "DWV-B load (ln)", title = NULL)+
  annotate("text", label = paste0("pd = ", round(pd$pd, 2)), x=mean(range(posterior_preds2$r, na.rm = TRUE)),  
           y=max(log(mod$data$DWVB.abs)), hjust=0.5)+
  geom_point(data = mod$data %>% filter(DWVB.abs> 0) %>%
                mutate(dwvb.f.r = dwvb.f.s * sd + m), aes(dwvb.f.r, log(DWVB.abs)), shape = 21,  alpha = 0.4)+
  theme_bw(base_size = 14)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))+
  add_phylopic(img = img_LS, x = min(posterior_preds2$r), y = max(log(mod$data$DWVB.abs)), height = 1.2, hjust = 0.02, vjust = 0.98)

##
mod = WB_models_combined_raw[['bqcv.b']]
pd = p_direction(mod) %>% filter(grepl("b_bqcv.f.s:", Parameter))

new_dat = expand.grid(
  bqcv.f.s = seq(min(mod$data$bqcv.f.s), max(mod$data$bqcv.f.s), by = 0.1),
  Morisita.hb.s = c(-1,1),
  Morisita.bl.s = 0,
  Connectance.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  Site = 'Goe1425',
  Species = 'Bombus pascuorum'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'mu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'bb') %>% ungroup()
m = mean(data.bb.raw$bqcv.f)
sd = sd(data.bb.raw$bqcv.f)
x = data.frame(z = new_dat$bqcv.f, r = new_dat$bqcv.f * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('bqcv.f.s' == 'z')) 


exp_bqcv_b = ggplot(posterior_preds2, aes(x = r, y = mu)) +
  labs(x = expression(paste("Viral exposure to ", italic("A. mellifera"))), 
       y = "BQCV load (ln)", title = NULL, fill = "Resource overlap")+
  annotate("text", label = paste0("pd = ", round(pd$pd[1],2)), x=mean(range(posterior_preds2$r, na.rm = TRUE)),  y=max(log(mod$data$BQCV.abs)), hjust=0.5)+
  geom_point(data = mod$data %>% filter(BQCV.abs > 0) %>%
                mutate(bqcv.f.r = bqcv.f.s * sd + m), aes(bqcv.f.r, log(BQCV.abs)), shape = 21, alpha = 0.3)+
  stat_lineribbon(.width = c(0.95), aes(fill = as.factor(Morisita.hb.s)), alpha = 0.7) +
  theme_bw(base_size = 14)+
  scale_fill_manual(values = c("#69547F","#6A7F54"), labels = c("0", "1"))+
  #scale_x_continuous(breaks = c(11,14,17,20))+
  # scale_y_continuous(breaks = c(0,0.5,1))+
  theme(legend.position = c(0.02, 0.07),
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        panel.grid = element_blank(), axis.text = element_text(color = 'black'))+
  guides(fill = guide_legend(nrow = 1, keywidth = 0.7, keyheight = 0.7))+
  add_phylopic(img = img_BT, x = min(posterior_preds2$r), y = max(log(mod$data$BQCV.abs))-0.5, 
               height = 1.5, hjust = 0.02)

##


mod = WB_models_combined_raw[['abpv.b']]
pd = p_direction(mod) %>% filter(grepl("hu_abpv.bl.f.s", Parameter))

new_dat = expand.grid(
  abpv.bl.f.s = seq(min(mod$data$abpv.bl.f.s), max(mod$data$abpv.bl.f.s), by = 0.1),
  Morisita.bl.s = 0,
  Morisita.hb.s = 0,
  Connectance.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  Site = "Goe1425",
  Species = 'Bombus pascuorum'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'bb' & Species != 'Bombus lapidarius') %>% ungroup()
m = mean(data.bb.raw$abpv.bl.f, na.rm = T)
sd = sd(data.bb.raw$abpv.bl.f, na.rm = T)
x = data.frame(z = new_dat$abpv.bl.f, r = new_dat$abpv.bl.f * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('abpv.bl.f.s' == 'z')) 

exp_abpv_b = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#d55e00') +
  labs(x = expression(paste("Viral exposure to ", italic("B. lapidarius"))), y = "ABPV prevalence", title = NULL)+
  annotate("text", label = paste0("pd = ", round(pd$pd, 2)), x=mean(range(posterior_preds2$r, na.rm = TRUE)),  y=0.92, hjust=0.5)+
  geom_jitter(data = mod$data %>% mutate(ABPV.abs = ifelse(ABPV.abs > 0, 1, 0)) %>%
                mutate(abpv.bl.f.r = abpv.bl.f.s * sd + m), aes(abpv.bl.f.r, ABPV.abs), shape = 21, height = 0.025, alpha = 0.4)+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme_bw(base_size = 14)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))+
  add_phylopic(img = img_BT, x = min(posterior_preds2$r), y = 0.85, height = 0.2, hjust = 0.02)

##

mod = WB_models_combined_raw[['bqcv.w']]
pd = p_direction(mod) %>% filter(grepl("b_bqcv.f.s:", Parameter))

new_dat = expand.grid(
  bqcv.f.s = seq(min(mod$data$bqcv.f.s), max(mod$data$bqcv.f.s), by = 0.1),
  Morisita.hb.s = c(-1,1),
  Morisita.bl.s = 0,
  Connectance.s = 0,
  total_bee_dens.s = 0,
  Year = 2022,
  Site = 'Goe1425',
  Species = 'Andrena flavipes'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'mu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'wb') %>% ungroup()
m = mean(data.bb.raw$bqcv.f)
sd = sd(data.bb.raw$bqcv.f)
x = data.frame(z = new_dat$bqcv.f, r = new_dat$bqcv.f * sd + m) %>% distinct(z,r)


posterior_preds2 = left_join(posterior_preds, x, by = join_by('bqcv.f.s' == 'z'))


exp_bqcv_w = ggplot(posterior_preds2, aes(x = r, y = mu)) +
  labs(x = expression(paste("Viral exposure to ", italic("A. mellifera"))),
       y = "BQCV prevalence", title = NULL, fill = "Resource overlap")+
  annotate("text", label = paste0("pd = ", round(pd$pd,2)), x=mean(range(posterior_preds2$r, na.rm = TRUE)),  y=max(log(mod$data$BQCV.abs)), hjust=0.5)+
  geom_point(data = mod$data %>% filter(BQCV.abs > 0) %>%
               mutate(bqcv.f.r = bqcv.f.s * sd + m), aes(bqcv.f.r, log(BQCV.abs)), shape = 21, alpha = 0.3)+
  stat_lineribbon(.width = c(0.95), aes(fill = as.factor(Morisita.hb.s)), alpha = 0.7) +
  theme_bw(base_size = 14)+
  scale_fill_manual(values = c("#69547F", "#6A7F54"), labels = c("0", "1"))+
  #scale_x_continuous(breaks = c(11,14,17,20))+
  # scale_y_continuous(breaks = c(0,0.5,1))+
  theme(legend.position = c(0.02, 0.02),
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill = alpha("white", 0)),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        panel.grid = element_blank(), axis.text = element_text(color = 'black'))+
  guides(fill = guide_legend(nrow = 1, keywidth = 0.7, keyheight = 0.7))+
  add_phylopic(img = img_LS, x = min(posterior_preds2$r), y = max(log(mod$data$BQCV.abs)),
               height = 1.5, hjust = 0.02, vjust = 0.98)

##

mod = WB_models_combined_raw[['abpv.h']]
pd = p_direction(mod) %>% filter(grepl("hu_abpv.bl.f.s", Parameter))

new_dat = expand.grid(
  abpv.bl.f.s = seq(min(mod$data$abpv.bl.f.s), max(mod$data$abpv.bl.f.s), by = 0.1),
  Morisita.bl.s = 0,
  Morisita.hb.s = 0,
  Connectance.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  Site = "Goe1425",
  Species = 'Bombus pascuorum'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'hb') %>% ungroup()
m = mean(data.bb.raw$abpv.bl.f, na.rm = T)
sd = sd(data.bb.raw$abpv.bl.f, na.rm = T)
x = data.frame(z = new_dat$abpv.bl.f, r = new_dat$abpv.bl.f * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('abpv.bl.f.s' == 'z')) 

exp_abpv_h = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#d55e00') +
  labs(x = expression(paste("Viral exposure to ", italic("B. lapidarius"))), y = "ABPV prevalence", title = NULL)+
  annotate("text", label = paste0("pd = ", round(pd$pd, 2)), x=mean(range(posterior_preds2$r, na.rm = TRUE)),  y=0.92, hjust=0.5)+
  geom_jitter(data = mod$data %>% mutate(ABPV.abs = ifelse(ABPV.abs > 0, 1, 0)) %>%
                mutate(abpv.bl.f.r = abpv.bl.f.s * sd + m), aes(abpv.bl.f.r, ABPV.abs), shape = 21, height = 0.025, alpha = 0.4)+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme_bw(base_size = 14)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))+
  add_phylopic(img = img_AM, x = min(posterior_preds2$r), y = 0.85, height = 0.2, hjust = 0.02)


exp_pred_6 = plot_grid(exp_dwvb_b, exp_dwvb_w,  exp_bqcv_b, exp_bqcv_w, exp_abpv_b, exp_abpv_h, align = 'hv', nrow = 3, ncol = 2)

## RESOURCE OVERLAP 

mod = WB_models_combined_raw[['dwvb.b']]
pd = p_direction(mod) %>% filter(grepl("hu_Morisita.hb.s", Parameter))

new_dat = expand.grid(
  Morisita.hb.s = seq(min(mod$data$Morisita.hb.s), max(mod$data$Morisita.hb.s), by = 0.1),
  Morisita.bl.s = 0,
  Connectance.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  dwvb.f.s = 0,
  bqcv.f.s = 0,
  abpv.bl.f.s = 0,
  Site = 'Goe1425',
  Species = 'Bombus lapidarius'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'bb') %>% ungroup()
m = mean(data.bb.raw$Morisita.hb, na.rm = T)
sd = sd(data.bb.raw$Morisita.hb, na.rm = T)
x = data.frame(z = new_dat$Morisita.hb.s, r = new_dat$Morisita.hb.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Morisita.hb.s' == 'z')) 


nich_dwvb_b = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#6a7f54') +
  labs(x = expression(paste("Resource overlap with ", italic("A. mellifera"))), y = "DWV-B prevalence", 
       title = expression(bold('a')))+
  annotate("text", label = paste0("pd = ", round(pd$pd, 2)), x=mean(range(posterior_preds2$r, na.rm = TRUE)),  y=0.92, hjust=0.5)+
  geom_jitter(data = mod$data %>% mutate(DWVB.abs = ifelse(DWVB.abs > 0, 1, 0)) %>%
                mutate(Morisita.hb.s.r = Morisita.hb.s * sd + m), aes(Morisita.hb.s.r, DWVB.abs),
              shape = 21, height = 0.025, alpha = 0.4)+
  theme_bw(base_size = 14)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'),
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+
  add_phylopic(img = img_BT, y = 0.85, x = min(posterior_preds2$r), height = 0.2, hjust = 0.02)

##

mod = WB_models_combined_raw[['bqcv.b']]
pd = p_direction(mod) %>% filter(grepl("hu_Morisita.hb.s", Parameter))

new_dat = expand.grid(
  Morisita.hb.s = seq(min(mod$data$Morisita.hb.s), max(mod$data$Morisita.hb.s), by = 0.1),
  Morisita.bl.s = 0,
  Connectance.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  dwvb.f.s = 0,
  bqcv.f.s = 0,
  abpv.bl.f.s = 0,
  Site = 'Goe1425',
  Species = 'Bombus pascuorum'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'bb') %>% ungroup()
m = mean(data.bb.raw$Morisita.hb, na.rm = T)
sd = sd(data.bb.raw$Morisita.hb, na.rm = T)
x = data.frame(z = new_dat$Morisita.hb.s, r = new_dat$Morisita.hb.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Morisita.hb.s' == 'z')) 

nich_bqcv_b = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#6a7f54') +
  labs(x = expression(paste("Resource overlap with ", italic("A. mellifera"))), y = "BQCV prevalence", 
       title = expression(bold('b')))+
  annotate("text", label = paste0("pd = ", round(pd$pd, 2)), x=mean(range(posterior_preds2$r, na.rm = TRUE)),  y=0.92, hjust=0.5)+
  geom_jitter(data = mod$data %>% mutate(BQCV.abs = ifelse(BQCV.abs > 0, 1, 0)) %>%
                mutate(Morisita.hb.s.r = Morisita.hb.s * sd + m), aes(Morisita.hb.s.r, BQCV.abs), 
              shape = 21, height = 0.025, alpha = 0.4)+ 
  theme_bw(base_size = 14)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'),
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+
  add_phylopic(img = img_BT, y = 0.85, x = min(posterior_preds2$r), height = 0.2, hjust = 0.02)

##

mod = WB_models_combined_raw[['abpv.w']]
pd = p_direction(mod) %>% filter(grepl("b_Morisita.bl.s", Parameter))

new_dat = expand.grid(
  Morisita.bl.s = seq(min(mod$data$Morisita.bl.s), max(mod$data$Morisita.bl.s), by = 0.1),
  Morisita.hb.s = 0,
  Connectance.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  dwvb.f.s = 0,
  bqcv.f.s = 0,
  abpv.bl.f.s = 0,
  Site = 'Goe1425',
  Species = 'Andrena flavipes'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'mu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'wb') %>% ungroup()
m = mean(data.bb.raw$Morisita.bl, na.rm = T)
sd = sd(data.bb.raw$Morisita.bl, na.rm = T)
x = data.frame(z = new_dat$Morisita.bl.s, r = new_dat$Morisita.bl.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Morisita.bl.s' == 'z')) 

nich_abpv_w = ggplot(posterior_preds2, aes(x = r, y = mu)) +
  stat_lineribbon(.width = c(0.95), fill = '#d55e00') +
  labs(x = expression(paste("Resource overlap with ", italic("B. lapidarius"))), y = "ABPV load (ln)", 
       title = expression(bold('c')))+
  annotate("text", label = paste0("pd = ", round(pd$pd,2)), x=mean(range(posterior_preds2$r, na.rm = T)),  y=12, hjust=0.5)+
  geom_jitter(data = mod$data[mod$data$ABPV.abs >0,] %>%
                mutate(Morisita.bl.s.r = Morisita.bl.s * sd + m), aes(Morisita.bl.s.r, log(ABPV.abs)), shape = 21, height = 0.025, 
              alpha = 0.4)+
  theme_bw(base_size = 14)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'),
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+
  add_phylopic(img = img_LS, y = max(log(mod$data$ABPV.abs)), x = min(posterior_preds2$r), 
               height = 1.2, hjust = 0.02, vjust = 0.98)

niche_pred = plot_grid(nich_dwvb_b , nich_bqcv_b,  nich_abpv_w, nrow = 1, align = 'hv')


### CONNECTANCE

mod = WB_models_combined_raw[['dwvb.b']]
pd = p_direction(mod) %>% filter(grepl("hu_Connectance.s", Parameter))

new_dat = expand.grid(
  Connectance.s = seq(min(mod$data$Connectance.s), max(mod$data$Connectance.s), by = 0.1),
  Morisita.bl.s = 0,
  Morisita.hb.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  dwvb.f.s = 0,
  bqcv.f.s = 0,
  abpv.bl.f.s = 0,
  Species = 'Bombus lapidarius'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'bb') %>% ungroup() 
m = mean(data.bb.raw$Connectance)
sd = sd(data.bb.raw$Connectance)
x = data.frame(z = new_dat$Connectance.s, r = new_dat$Connectance.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Connectance.s' == 'z')) 

con.a = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#547f7f') +
  labs(x = "Connectance", y = "DWV-B prevalence", title = expression(bold('a')))+
  geom_jitter(data = mod$data %>% mutate(Connectance.r = Connectance.s * sd + m) , aes(Connectance.r, ifelse(DWVB.abs > 0, 1, 0)), shape = 21, height = 0.025, alpha = 0.2)+
  annotate("text", label = paste0("pd = ", round(pd$pd,2)), x=mean(range(posterior_preds2$r), na.rm = T),  y=0.92, hjust=0.5)+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'),
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+
  add_phylopic(img = img_BT, x = min(posterior_preds2$r), y = 0.85, height = 0.2, hjust = 0.2)

##

mod = WB_models_combined_raw[['dwvb.w']]
pd = p_direction(mod) %>% filter(grepl("hu_Connectance.s", Parameter))

new_dat = expand.grid(
  Connectance.s = seq(min(mod$data$Connectance.s), max(mod$data$Connectance.s), by = 0.1),
  Morisita.bl.s = 0,
  Morisita.hb.s = 0,
  total_bee_dens.s = 0,
  Year = 2022,
  dwvb.f.s = 0,
  bqcv.f.s = 0,
  abpv.bl.f.s = 0,
  Species = 'Andrena flavipes'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'wb') %>% ungroup()
m = mean(data.bb.raw$Connectance)
sd = sd(data.bb.raw$Connectance)
x = data.frame(z = new_dat$Connectance.s, r = new_dat$Connectance.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Connectance.s' == 'z'))

con.b = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#547f7f') +
  labs(x = "Connectance", y = "DWV-B prevalence", title = expression(bold('b')))+
  geom_jitter(data = mod$data %>% mutate(Connectance.r = Connectance.s * sd + m) , aes(Connectance.r, ifelse(DWVB.abs > 0, 1, 0)), shape = 21, height = 0.025, alpha = 0.2)+
  annotate("text", label = paste0("pd = ", round(pd$pd,2)), x=mean(range(posterior_preds2$r), na.rm = T),  y=0.92, hjust=0.5)+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'),
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot",
        plot.caption.position =  "plot")+
  add_phylopic(img = img_LS, x = min(posterior_preds2$r), y = 0.85, height = 0.2, hjust = 0.2)

##

mod = WB_models_combined_raw[['bqcv.b']]
pd = p_direction(mod) %>% filter(grepl("hu_Connectance.s", Parameter))

new_dat = expand.grid(
  Connectance.s = seq(min(mod$data$Connectance.s), max(mod$data$Connectance.s), by = 0.1),
  Morisita.bl.s = 0,
  Morisita.hb.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  dwvb.f.s = 0,
  bqcv.f.s = 0,
  abpv.bl.f.s = 0,
  Species = 'Bombus pascuorum'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'bb') %>% ungroup()
m = mean(data.bb.raw$Connectance)
sd = sd(data.bb.raw$Connectance)
x = data.frame(z = new_dat$Connectance.s, r = new_dat$Connectance.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Connectance.s' == 'z')) 

con.c = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#547f7f') +
  labs(x = "Connectance", y = "BQCV prevalence", title = expression(bold('c')))+
  geom_jitter(data = mod$data %>% mutate(Connectance.r = Connectance.s * sd + m) , aes(Connectance.r, ifelse(BQCV.abs > 0, 1, 0)), shape = 21, height = 0.025, alpha = 0.2)+
  annotate("text", label = paste0("pd = ", round(pd$pd,2)), x=mean(range(posterior_preds2$r), na.rm = T),  y=0.92, hjust=0.5)+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'),
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+
  add_phylopic(img = img_BT, x = min(posterior_preds2$r), y = 0.85, height = 0.2, hjust = 0.2)


##

mod = WB_models_combined_raw[['bqcv.w']]
pd = p_direction(mod) %>% filter(grepl("b_Connectance.s", Parameter))

new_dat = expand.grid(
  Connectance.s = seq(min(mod$data$Connectance.s), max(mod$data$Connectance.s), by = 0.1),
  Morisita.bl.s = 0,
  Morisita.hb.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  dwvb.f.s = 0,
  bqcv.f.s = 0,
  abpv.bl.f.s = 0,
  Species = 'Andrena flavipes'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'mu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'wb') %>% ungroup()
m = mean(data.bb.raw$Connectance)
sd = sd(data.bb.raw$Connectance)
x = data.frame(z = new_dat$Connectance.s, r = new_dat$Connectance.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Connectance.s' == 'z')) 

con.b2 = ggplot(posterior_preds2, aes(x = r, y = mu)) +
  stat_lineribbon(.width = c(0.95), fill = '#547f7f') +
  labs(x = "Connectance", y = "BQCV load (ln)", title = expression(bold('b')))+
  geom_point(data = mod$data %>% filter(BQCV.abs > 0) %>% mutate(Connectance.r = Connectance.s * sd + m) , aes(Connectance.r, log(BQCV.abs)), shape = 21, alpha = 0.2)+
  annotate("text", label = paste0("pd = ", round(pd$pd,2)), x=mean(range(posterior_preds2$r), na.rm = T),  y=max(log(mod$data$BQCV.abs)), hjust=0.5)+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'),
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+
  add_phylopic(img = img_LS, x = min(posterior_preds2$r), y = max(log(mod$data$BQCV.abs)), height = 1.2, hjust = 0.2, vjust = 0.98)

##

mod = WB_models_combined_raw[['bqcv.w']]
pd = p_direction(mod) %>% filter(grepl("hu_Connectance.s", Parameter))

new_dat = expand.grid(
  Connectance.s = seq(min(mod$data$Connectance.s), max(mod$data$Connectance.s), by = 0.1),
  Morisita.bl.s = 0,
  Morisita.hb.s = 0,
  total_bee_dens.s = 0,
  Year = 2022, 
  dwvb.f.s = 0,
  bqcv.f.s = 0,
  abpv.bl.f.s = 0,
  Species = 'Andrena flavipes'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'wb') %>% ungroup()
m = mean(data.bb.raw$Connectance)
sd = sd(data.bb.raw$Connectance)
x = data.frame(z = new_dat$Connectance.s, r = new_dat$Connectance.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Connectance.s' == 'z')) 

con.d = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#547f7f') +
  labs(x = "Connectance", y = "ABPV prevalence", title = expression(bold('b')))+
  geom_jitter(data = mod$data %>% mutate(Connectance.r = Connectance.s * sd + m) , aes(Connectance.r, ifelse(BQCV.abs > 0, 1, 0)), shape = 21, height = 0.025, alpha = 0.2)+
  annotate("text", label = paste0("pd = ", round(pd$pd,2)), x=mean(range(posterior_preds2$r), na.rm = T),  y=0.92, hjust=0.5)+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'),
        plot.caption = element_text(hjust = 0),
        plot.title.position = "plot", 
        plot.caption.position =  "plot")+
  add_phylopic(img = img_LS, x = min(posterior_preds2$r), y = 0.85, height = 0.2, hjust = 0.2)


con.grid = plot_grid(con.a, con.b, con.c, ncol = 2, align = 'hv')

###############################
######### FIGURES MS ##########
###############################

## FIGURE 1 - made in powerpoint


## FIGURE 2
plot_grid(r0_plots, exp_pred_6, rel_widths = c(1.3/3, 1.6/3), ncol = 2, align = 'hv', scale = 0.95)

ggsave(file = 'Data/fig/f2.png', height = 10, width = 12)
ggsave(file = 'Data/fig/f2.pdf', height = 10, width = 12)

## FIGURE 3
heat

ggsave(file = 'Data/fig/f3.pdf', width = 9, height = 8)
ggsave(file = 'Data/fig/f3.png', width = 9, height = 8)

## FIGURE 4
niche_pred

ggsave('Data/fig/f4.pdf', width = 10, height = 4)
ggsave('Data/fig/f4.png', width = 10, height = 4, bg = 'transparent')


## FIGURE 5
plot_grid(con.a, con.b, con.c, prev.sep, nrow = 2,  align = 'hv')

ggsave('Data/fig/f5.pdf', width = 8, height = 8)
ggsave('Data/fig/f5.png', width = 8, height = 8)

## FIGURE S1 
# takes a long time to query osm
source("09_study_map.R")

plot_grid(map, empty, empty, hb.dens.colony,  ncol = 2, nrow = 2, rel_widths = c(2/3, 1/3), rel_heights = c(3/4, 1/4))

# Germany's outline added manually - osm query too big
ggsave(file = 'Data/fig/fs1.png', height = 8, width = 8)
ggsave(file = 'Data/fig/fs1.pdf', height = 8, width = 8)

## FIGURE S2
prevload.raw

ggsave('Data/fig/fs2.pdf', width = 8, height = 10)
ggsave('Data/fig/fs2.png', width = 8, height = 10)

## FIGURE S3
r0_overlap_plot

ggsave('Data/fig/fs3.pdf', width = 8, height = 12)
ggsave('Data/fig/fs3.png', width = 8, height = 12)


## FIGURE SA1 (appendix)

ggplot(network.metrics.both, aes(x = Connectance, y = NODF)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top") +
  theme_bw(base_size = 16)+
  labs(x = "Connectance", y = "NODF")

ggsave(file = "Data/fig/a1.png", height = 4, width = 4)

## FIGURE SA2 (appendix)

nodes_con = ggplot(network.metrics.both, aes(x = sum_nodes, y = Connectance)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top") +
  theme_bw(base_size = 16)+
  labs(x = "Sum of nodes", y = "Connectance", title = expression(bold("a")))

nodes_mo = morisita.raw %>%
  select(Species, Site, Year, `Apis mellifera`, `Bombus lapidarius`) %>%
  pivot_longer(cols = c(`Apis mellifera`, `Bombus lapidarius`)) %>% 
  group_by(Site, Year) %>%
  summarise(mean_morisita = mean(value)) %>%
  left_join(network.metrics.both %>% select(Site, Year, sum_nodes), by = c("Site", "Year")) %>%
  ggplot(aes(x = sum_nodes, y = mean_morisita)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top") +
  theme_bw(base_size = 16)+
  labs(x = "Sum of nodes", y = "Mean Morisita index", title = expression(bold("b")))

plot_grid(nodes_con, nodes_mo)

ggsave(file = "Data/fig/a2.png", width = 8, height = 4)

##

## Figure for R2:

gmean_con = ggplot(network.metrics.both, aes(x = sqrt(bee_rich*plant_rich), y = Connectance)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = "top") +
  theme_bw(base_size = 16)+
  labs(x = "Geometric mean of\nplants and pollinators", y = "Connectance", title = expression(bold("b")))

plot_grid(nodes_con, gmean_con)

ggsave(file = "Data/fig/r2.png", width = 8, height = 4)

