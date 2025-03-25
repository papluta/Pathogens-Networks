load('Data/Results/250260R0_results.RData')
load('Data/Results/250304_models_mor_con_sep.RData')
load('Data/Results/250304_models_connectance.RData')
load('Data/Results/250130_comm_models.RData')

source('05_main_data_file.R')
source('00_plotting_functions.R')

library(dplyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(bayestestR)

#######################################
####### RAW PREVALENCE AND LOAD #######
#######################################


p = data.pathogen.both  %>% pivot_longer(cols = dwvb:sbv, names_to = 'virus', values_to = 'presence') %>% 
  filter(Species %in% c(unique(data.wb$Species), unique(data.bb.forhb$Species),'Apis mellifera')) %>%
  group_by(Group, virus, Year) %>%
  summarise(prevalence = mean(presence), n = n()) %>% mutate(Group = factor(Group, levels = c('hb', 'bb', 'wb'))) %>% 
  mutate(virus = case_match(virus, 'dwvb' ~ 'DWV-B', 'bqcv' ~ 'BQCV', 'abpv' ~ 'ABPV')) %>%
  mutate(virus = factor(virus, levels = c('DWV-B', 'BQCV', 'ABPV'))) %>% filter(virus != 'sbv') %>%
  ggplot(aes(virus, prevalence, fill = Group)) +
  geom_bar(stat = 'identity', col = 'black', position = position_dodge())+
  scale_fill_manual(values = c('#6a7f54', '#d55e00', '#f0e442'), labels = c('Honey bee', 'Bumble bee', 'Other bee'))+
  facet_wrap(~Year, nrow = 1)+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color = 'black'), legend.position = 'bottom')+  
  labs(x = NULL, y = 'Prevalence', fill = 'Bee group', title = 'A')


l = data.pathogen.both %>% mutate(across(DWVB.abs:SBV.abs, function(x) x/BUFFER)) %>% pivot_longer(cols = DWVB.abs:SBV.abs, names_to = 'virus', values_to = 'load') %>% 
  filter(Species %in% c(unique(data.wb$Species), unique(data.bb.forhb$Species),'Apis mellifera')) %>%
  mutate(Group = factor(Group, levels = c('hb', 'bb', 'wb'))) %>% filter(virus != 'SBV.abs') %>% 
  mutate(virus = factor(virus, levels = c('DWVB.abs', 'BQCV.abs', 'ABPV.abs'))) %>% 
  filter(load > 0) %>% 
  ggplot(aes(virus, log10(load), fill = Group))+
  geom_boxplot(col = 'black')+
  theme_classic(base_size = 18)+
  theme(axis.text = element_text(color = 'black'), legend.position = 'none')+
  labs(x = NULL, y = 'Number of viral copies per uL (log10)', fill = 'Bee group', title = 'B')+
  scale_fill_manual(values = c('#6a7f54', '#d55e00', '#f0e442'), labels = c('Honey bee', 'Bumble bee', 'Other bee'))+
  scale_x_discrete(labels = c('DWV-B', 'BQCV', 'ABPV'))+
  facet_wrap(~Year)


leg.p = get_plot_component(p, "guide-box-bottom", return_all = T)

prevload.raw = plot_grid(p + theme(legend.position = 'none'), l,leg.p, rel_heights = c(0.4,0.4,0.2), nrow = 3)

ggsave('Data/fig/prev_raw.pdf', width = 8, height = 10)
ggsave('Data/fig/prev_raw.png', width = 8, height = 10)


dw = data.both  %>% pivot_longer(cols = dwvb:sbv, names_to = 'virus', values_to = 'presence') %>% group_by(Group, virus) %>%
  summarise(prevalence = mean(presence)) %>% mutate(Group = factor(Group, levels = c('hb', 'bb', 'wb'))) %>% 
  mutate(virus = factor(virus, levels = c('bqcv', 'dwvb', 'abpv', 'sbv'))) %>% filter(virus == 'dwvb') %>%
  ggplot(aes(virus, prevalence, fill = Group)) +
  geom_bar(stat = 'identity', col = 'black', position = position_dodge(), width = 1)+
  scale_fill_manual(values = c('#6a7f54', '#d55e00', '#f0e442'), labels = c('Honeybee', 'Bumblebee', 'Other bee'))+
  theme_bw(base_size = 16)+
  geom_hline(aes(yintercept = mean(prevalence)), linetype = 'dashed')+
  ylim(0,1)+
  labs(x = NULL, y = 'Prevalence', fill = NULL, title = 'C')+
  theme(axis.text = element_text(color = 'black'), legend.position = 'left', panel.grid = element_blank(), 
        axis.ticks.x = element_blank())+
  scale_x_discrete(labels = c('DWV-B'))

bq = data.both  %>% pivot_longer(cols = dwvb:sbv, names_to = 'virus', values_to = 'presence') %>% group_by(Group, virus) %>%
  summarise(prevalence = mean(presence)) %>% mutate(Group = factor(Group, levels = c('hb', 'bb', 'wb'))) %>% 
  mutate(virus = factor(virus, levels = c('bqcv', 'dwvb', 'abpv', 'sbv'))) %>% filter(virus == 'bqcv') %>%
  ggplot(aes(virus, prevalence, fill = Group)) +
  geom_bar(stat = 'identity', col = 'black', position = position_dodge(), width = 1)+
  scale_fill_manual(values = c('#6a7f54', '#d55e00', '#f0e442'), labels = c('Honeybee', 'Bumblebee', 'Other bee'))+
  geom_hline(aes(yintercept = mean(prevalence)), linetype = 'dashed')+
  theme_bw(base_size = 16)+
  ylim(0,1)+
  labs(x = NULL, y = 'Prevalence', fill = NULL, title = 'F')+
  theme(axis.text = element_text(color = 'black'), legend.position = 'none', panel.grid = element_blank(), 
        axis.ticks.x = element_blank())+
  scale_x_discrete(labels = c('BQCV'))

leg.p = get_legend(dw)
leg.p = get_plot_component(dw, 'guide-box-left', return_all = T)
prev.sep = plot_grid(dw + theme(legend.position = 'none'),leg.p, bq, rel_heights = c(0.4,0.2,0.4), nrow = 3)


########################
########## R0 ##########
########################

## R0 result table

r0.comb = r0.main.host.sp %>% bind_rows(.id = 'mod') %>% 
  mutate(across(r0:prev, function(x) round(x, digits = 2))) %>%
  mutate(Virus = sub('.main.202..sp',  '', sub('r0.', '', mod)), Year = ifelse(grepl('2021', mod), 2021, 2022)) %>%
  select(-mod)

ms.r0 = r0.comb %>% select(-starts_with('prev')) %>% pivot_wider(values_from = c(r0, sd), names_from = c(Virus, Year))

#write.csv(ms.r0, file = 'Data/Results/r0_results.csv')

## prevalence simulation result table

sim.ms = r0.main.host %>% bind_rows(.id = 'mod') %>% 
  mutate(across(mean.r0:mean.prev, function(x) round(x, digits = 2))) %>%
  mutate(Virus = sub('.main.202.',  '', sub('r0.', '', mod)), Year = ifelse(grepl('2021', mod), 2021, 2022)) %>%
  select(-mod) %>% mutate(Species = sub(' ', '\n', Species)) %>% group_by(Species, Virus) %>%
  summarise(prev.no.mh.m = mean(mean.prev.no.mh), prev.no.mh.sd = sd(mean.prev.no.mh), prev.m = mean(mean.prev), prev.sd = sd(mean.prev)) %>%
  rename(sim.m = prev.no.mh.m, sim.sd = prev.no.mh.sd) %>% pivot_wider(names_from = Virus, values_from = c(sim.m, sim.sd, prev.m, prev.sd)) %>% arrange(Species)

#write.csv(sim.ms, file = 'Data/Results/sim_results_r0.csv')

sim.ms.alt = r0.alt.host %>% bind_rows(.id = 'mod') %>% 
  mutate(across(mean.r0:mean.prev, function(x) round(x, digits = 2))) %>%
  mutate(Virus = sub('.alt.202.',  '', sub('r0.', '', mod)), Year = ifelse(grepl('2021', mod), 2021, 2022)) %>%
  select(-mod) %>% group_by(Species, Virus) %>%
  summarise(prev.no.mh.m = mean(mean.prev.no.mh), prev.no.mh.sd = sd(mean.prev.no.mh)) %>%
  rename(sim.m = prev.no.mh.m, sim.sd = prev.no.mh.sd) %>% pivot_wider(names_from = Virus, values_from = c(sim.m, sim.sd)) %>% arrange(Species)

#write.csv(sim.ms.alt, file = 'Data/Results/sim_results_r0_althost.csv')

### plotting top 4 R0 species

r0.plot = r0.main.host %>% bind_rows(.id = 'mod') %>% 
  mutate(across(mean.r0:mean.prev, function(x) round(x, digits = 2))) %>%
  mutate(Virus = sub('.main.202.',  '', sub('r0.', '', mod)), Year = ifelse(grepl('2021', mod), 2021, 2022)) %>%
  select(-mod) %>% mutate(Species = sub(' ', '\n', Species))

top4 = r0.plot %>% group_by(Virus, Species) %>% summarise(r0.m = mean(mean.r0, na.rm = T), sd = round(sd(mean.r0, na.rm = T), digits = 2)) %>%
  arrange(Virus, desc(r0.m)) %>% slice_max(r0.m, n = 4) 

dwvb.top.r0 = r0.plot %>% filter(Virus == 'dwvb') %>% filter(Species %in% top4[top4$Virus == 'dwvb',]$Species) %>%
  ggplot(aes(Species, mean.r0))+
  geom_boxplot(aes(fill = Species), lwd = 0.8)+
  geom_abline(aes(intercept = 1, slope = 0), linetype = 21, linewidth = 1, col = 'black')+
  theme_classic(base_size = 16)+
  scale_fill_manual(values = c('#E69F00', '#6a7f54', '#d55e00', '#56B4E9'))+
  theme(legend.position = 'none', axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', size = 14, face = 'italic'), title = element_text(size = 20))+
  labs(x = NULL, y = expression(paste('DWV-B ', R[0][","][i])), title = 'A')+
  ylim(0,3)

bqcv.top.r0 = r0.plot %>% filter(Virus == 'bqcv') %>% filter(Species %in% top4[top4$Virus == 'bqcv',]$Species) %>%
  ggplot(aes(Species, log10(mean.r0+1)))+
  geom_boxplot(aes(fill = Species), lwd = 0.8)+
  geom_abline(aes(intercept = 0, slope = 0), linetype = 21, linewidth = 1, col = 'black')+
  theme_classic(base_size = 16)+
  scale_fill_manual(values = c('#6a7f54', '#d55e00', '#56B4E9', '#CC79A7'))+
  theme(legend.position = 'none', axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', size = 14, face = 'italic'), title = element_text(size = 20))+
  labs(x = NULL, y = expression(paste('BQCV ', R[0][","][i])), title = 'B')+
  scale_y_continuous(labels = c(1,10,100,1000), breaks = c(0,1,2,3))


abpv.top.r0 = r0.plot %>% filter(Virus == 'abpv') %>% filter(Species %in% top4[top4$Virus == 'abpv',]$Species) %>%
  ggplot(aes(Species, mean.r0))+
  geom_boxplot(aes(fill = Species), lwd = 0.8)+
  geom_abline(aes(intercept = 1, slope = 0), linetype = 21, linewidth = 1, col = 'black')+
  theme_classic(base_size = 16)+
  scale_fill_manual(values = c('#6a7f54', '#d55e00', '#56B4E9', '#CC79A7'))+
  theme(legend.position = 'none', axis.text.y = element_text(color = 'black'),
        axis.text.x = element_text(color = 'black', size = 14, face = 'italic'), title = element_text(size = 20))+
  labs(x = NULL, y = expression(paste('ABPV ', R[0][","][i])), title = 'C')+
  ylim(0,6)

r0_plots = plot_grid(dwvb.top.r0, bqcv.top.r0, abpv.top.r0, align = 'hv', nrow = 3)


## simulation plot

obs = sim.ms %>% select(Species, starts_with('prev.m_')) %>% pivot_longer(cols = starts_with('prev.m_'), names_to = 'Virus', values_to = 'Prev') %>% mutate(host = 'obs', Virus = sub('prev.m_','', Virus))
sim1 = sim.ms %>% select(Species, starts_with('sim.m_')) %>% pivot_longer(cols = starts_with('sim.m_'), names_to = 'Virus', values_to = 'Prev') %>% mutate(host = 'main', Virus = sub('sim.m_','', Virus))
sim2 = sim.ms.alt %>% select(Species, starts_with('sim.m_')) %>% pivot_longer(cols = starts_with('sim.m_'), names_to = 'Virus', values_to = 'Prev') %>% mutate(host = 'alt', Virus = sub('sim.m_','', Virus)) %>% mutate(Species = sub(' ', '\n', Species))
sim = rbind(sim1, sim2) %>% left_join(obs, by = c('Species', 'Virus')) %>% mutate(fold = Prev.x/Prev.y) %>% 
  mutate(host.v = paste0(Virus, host.x)) %>% mutate(Species = sub('\n',' ', Species)) %>% mutate(host.v = factor(host.v, levels = c('abpvmain', 'abpvalt', 'bqcvalt', 'bqcvmain', 'dwvbalt', 'dwvbmain')))

ggplot(sim, aes(Species, host.v, fill= fold)) + 
  geom_tile()+
  labs(x = NULL, y = NULL, fill = 'Fold change from\nobserved prevalence')+
  theme_map()+
  scale_fill_viridis(option = 'magma')+
  scale_y_discrete(labels = c('B. lapidarius', 'A. mellifera','B. lapidarius', 'A. mellifera','B. lapidarius', 'A. mellifera' ))+
  theme(axis.ticks = element_blank(), axis.text = element_text(color = 'black', face = 'italic'), axis.text.x = element_text(angle = 90, hjust = 0.98))
  
ggsave(file = 'Data/fig/heatmap_simprev.pdf', width = 9, height = 8)
ggsave(file = 'Data/fig/heatmap_simprev.png', width = 9, height = 8)

########################################
########## STATISTICAL MODELS ##########
########################################

########################
###### COMMUNITY #######
########################
ggplot(data.site, aes(x = as.factor(Density), y = hb.dens, fill = as.factor(Density))) +
  geom_boxplot()+
  theme_bw()+
  scale_fill_manual(values = c('#d0d316', '#d32116'))+
  labs(x = '\nHoneybee colony density', y = expression(paste('Honeybee density ', 'per ', m^2 ,'of flowers')),
       title = 'C')+
  theme(panel.grid = element_blank(), axis.text = element_text(color = 'black'), legend.position = 'none',
        axis.ticks.x = element_blank())+
  scale_x_discrete(labels = c('Low', 'Increased'))

ggsave(file = 'Data/fig/colony_dens_pred.png', height = 3, width = 3)
ggsave(file = 'Data/fig/colony_dens_pred.pdf', height = 3, width = 3)

### JOINT SUMMARY OF THE MODELS
cs.site <- lapply(site_models_500, custom_summary) %>% bind_rows(.id = 'mod') %>%
  mutate(Predictor = recode_factor(as.factor(Predictor), 'FL_rich' = 'Flower richness', 
                                   'FL_per_agr:HB_per_agr' = 'HB density:Flower cover',  'Bee_rich' = 'Bee richness'),
         Response = recode_factor(as.factor(Response), 'logsumab' = 'Bee density', 'Morisitaz' = 'HB niche overlap',
                                  'Morisitablz' = 'BL niche overlap', 'Beerich' = 'Bee richness')) %>% select(-mod)

write.csv(cs.site, file = 'Data/Results/site_models.csv', row.names = F)


mod = site_models_500[['morisita.hb.m']]
p_direction(mod)
conditional_effects(mod)
new_dat = expand.grid(
  FL_rich = c(4,11),
  sum.ab.fl = seq(min(mod$data$sum.ab.fl), max(mod$data$sum.ab.fl), by = 0.05),
  Year = 2022
)

posterior_preds <- epred_draws(mod, newdata = new_dat, re_formula = NA)

mor.a = ggplot(posterior_preds, aes(x = sum.ab.fl, y = .epred)) +
  stat_lineribbon(.width = c(0.95), alpha = 0.5, aes(fill = as.factor(FL_rich)), show.legend = F) +
  labs(x = "Total bee density", y = "Resource overlap with\nhoneybee", title = 'A', col = 'Flower richness')+
  geom_jitter(data = mod$data, aes(sum.ab.fl, Morisita.z, col = cut(FL_rich, 2)), shape = 19, height = 0.025, alpha = 0.5)+ 
  #scale_y_continuous(limits = c(-5, 2.5))+
  theme_bw(base_size = 16)+
  scale_color_manual(values = c('#ddcd4b','#4b5bdd'))+
  scale_fill_manual(values = c('#ddcd4b','#4b5bdd'))+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'))


mod = site_models_500[['morisita.bl.m']]
p_direction(mod)
#conditional_effects(mod)
new_dat = expand.grid(
  FL_rich = c(4,11),
  sum.ab.fl = seq(min(mod$data$sum.ab.fl), max(mod$data$sum.ab.fl), by = 0.05),
  Year = 2022
)

posterior_preds <- epred_draws(mod, newdata = new_dat, re_formula = NA)


mor.b = ggplot(posterior_preds, aes(x = sum.ab.fl, y = .epred)) +
  stat_lineribbon(.width = c(0.95), alpha = 0.5, aes(fill = as.factor(FL_rich)), show.legend = F) +
  labs(x = "Total bee density", y = "Resource overlap with\nred-tailed bumblebee", title = 'B', col = 'Flower richness')+
  geom_jitter(data = mod$data, aes(sum.ab.fl, Morisita.bl.z, col = cut(FL_rich, 2)), shape = 19, height = 0.025, alpha = 0.5)+ 
  theme_bw(base_size = 14)+
  scale_y_continuous(limits = c(-5, 2.5))+
  scale_color_manual(values = c('#ddcd4b','#4b5bdd'))+
  scale_fill_manual(values = c('#ddcd4b','#4b5bdd'))+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))


legend <- get_legend(
  mor.a + theme(legend.box.margin = margin(-200, 0, 0, 0))
)

legend <- get_plot_component(mor.a,  'guide-box-bottom', return_all = T)

aa = plot_grid(mor.a + theme(legend.position="none"), mor.b, align = 'hv')

plot_grid(aa, legend, ncol = 1)

ggsave(file = 'Data/fig/pred_site_mor.pdf', height = 8, width = 8)
ggsave(file = 'Data/fig/pred_site_mor.png', height = 8, width = 8)

### Connectance


mod = site_models_500[['con.m']]
p_direction(mod)
#conditional_effects(mod)
#pairs(mod)

new_dat = expand.grid(
  FL_rich = mean(mod$data$FL_rich),
  sum.ab.fl.s = seq(min(mod$data$sum.ab.fl.s), max(mod$data$sum.ab.fl.s), by = 0.05),
  Year = 2021
)

posterior_preds <- epred_draws(mod, newdata = new_dat, re_formula = NA)

m = mean(data.site$sum.ab.fl)
sd = sd(data.site$sum.ab.fl)
x = data.frame(z = new_dat$sum.ab.fl.s, r = new_dat$sum.ab.fl.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('sum.ab.fl.s' == 'z')) 

con.a = ggplot(posterior_preds2, aes(x = r, y = .epred)) +
  stat_lineribbon(.width = c(0.95), fill = '#4b5bdd', show.legend = F) +
  labs(x = "Total bee density (ln)", y = "Network connectance", title = 'A')+
  geom_jitter(data = mod$data  %>%
                mutate(sum.ab.fl.r = sum.ab.fl.s * sd + m), aes(sum.ab.fl.r, Connectance), shape = 19, height = 0.025, alpha = 0.5)+ 
  #scale_y_continuous(limits = c(-5, 2.5))+
  theme_bw(base_size = 16)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))


new_dat = expand.grid(
  FL_rich = seq(min(mod$data$FL_rich), max(mod$data$FL_rich), by = 0.05),
  sum.ab.fl.s = 0,
  Year = 2022
)

posterior_preds <- epred_draws(mod, newdata = new_dat, re_formula = NA)


con.b = ggplot(posterior_preds, aes(x = FL_rich, y = .epred)) +
  stat_lineribbon(.width = c(0.95), fill = '#4b5bdd', show.legend = F) +
  labs(x = "Flower richness", y = "Network connectance", title = 'B')+
  geom_jitter(data = mod$data, aes(FL_rich, Connectance), shape = 19, height = 0.025, alpha = 0.5)+ 
  #scale_y_continuous(limits = c(-5, 2.5))+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'))

plot_grid(con.a, con.b + theme(legend.position = 'none'), ncol = 2, align = 'hv')

ggsave(file = 'Data/fig/connectance_pred.png', width = 8, height = 4)

########################
#### hurdle summary ####
########################

wb.cs <- lapply(WB_models_Connectance, custom_summary) %>% bind_rows(.id = 'id') 
wb.cs %>% filter(pd > 0.91) %>% 
  filter(Predictor != 'Intercept' & Predictor != 'hu_Intercept'& Predictor != 'SpeciesBombuspascuorum' & 
           Predictor != 'hu_SpeciesBombuspascuorum' & Predictor != 'SpeciesBombusterrestris' & Predictor != 'hu_SpeciesBombusterrestris') %>% 
  arrange(Response, Predictor) %>% rownames_to_column() %>% select(-rowname) %>% filter(Predictor != 'Year2022' & Predictor != 'hu_Year2022')

wb.cs2 <- wb.cs %>% mutate(term = ifelse(grepl('hu_', Predictor), 'Probability of presence', 'Viral load')) %>% 
  mutate(Predictor = sub('hu_', '', Predictor), Bee_group = ifelse(grepl('.h$', id), 'Honey bee', 'Bumble bee')) %>% 
  mutate(Bee_group = ifelse(grepl('(.w)$', id), 'Other wild bee', Bee_group)) %>% rename(Virus = Response) %>%
  mutate(Predictor = recode_factor(as.factor(Predictor), 'FL_per_agr' = 'Flower cover', 'dwvb.f' = 'Transmission potential', 'bqcv.f' = 'Transmission potential', 'abpv.bl.f' = 'Transmission potential', 'sbv.bl.f' = 'Transmission potential', 
                                   'Bee_rich' = 'Bee richness', 'Closeness.BL' = 'BL closeness', 'Closeness.AP' = 'HB closeness', 'Year2022' = 'Year 2022', 'Morisita.z' = 'HB niche overlap', 'Morisita.bl.z' = 'BL niche overlap',
                                   'SpeciesBombuspascuorum' = 'B. pascuorum', 'SpeciesBombusterrestris' = 'B. terrestris', 'dwvb.f:Morisita.z' = 'Transmission potential : niche overlap',
                                   'abpv.f:Morisita.bl.z' = 'Transmission potential : niche overlap', 'bqcv.f:Morisita.z' = 'Transmission potential : niche overlap', 'sum.ab.fl' = 'Total bee density'),
         Virus = recode_factor(as.factor(Virus), 'DWVBabs' = 'DWV-B', 'BQCVabs' = 'BQCV', 'ABPVabs' = 'ABPV', 'SBVabs' = 'SBV')) %>% select(-id) %>% 
  mutate(Virus = factor(Virus, levels = c('DWV-B', 'BQCV', 'ABPV', 'SBV')), Bee_group = factor(Bee_group, levels = c('Honey bee', 'Bumble bee', 'Other wild bee')),
         Predictor = factor(Predictor, levels = c('Intercept', 'Transmission potential', 'HB niche overlap', 'BL niche overlap', 'Total bee density','Connectance','HB closeness', 'BL closeness','Bee richness', 'Transmission potential : niche overlap', 'Transmission potential : closeness',
                                                  'Flower cover', 'Year 2022', 'B. pascuorum', 'B. terrestris'))) %>%
  relocate(Virus, Bee_group, term, Predictor) %>% arrange(Virus, Bee_group, term, Predictor)

write.csv(wb.cs2, file = 'Data/Results/250226bee_models_morisita.csv', row.names = F)
write.csv(wb.cs2, file = 'Data/Results/25026bee_models_connectance.csv', row.names = F)



### PLOTTING ###

## VIRAL EXPOSURE
mod = WB_models_Morisita[['dwvb.b']]
p_direction(mod)


new_dat = expand.grid(
  dwvb.f.s = seq(min(mod$data$dwvb.f.s), max(mod$data$dwvb.f.s), by = 0.1),
  Morisita.bl.z.s = 0,
  Morisita.z.s = 0,
  Year = 2022, 
  Site = 'Goe1425',
  Species = 'Bombus lapidarius'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)
data.bb.raw = data.both %>% filter(Group == 'bb') %>% ungroup() %>%
  filter(Morisita.z > -10)
m = mean(data.bb.raw$dwvb.f)
sd = sd(data.bb.raw$dwvb.f)
x = data.frame(z = new_dat$dwvb.f, r = new_dat$dwvb.f * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('dwvb.f.s' == 'z')) 


exp_dwvb_w = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#6a7f54') +
  labs(x = expression(paste("Viral exposure to ", italic("A. mellifera"))), y = "DWV-B prevalence", title = NULL)+
  annotate("text", label = "pd = 1.00", x=11,  y=0.92, hjust=0)+
  geom_jitter(data = mod$data %>% mutate(DWVB.abs = ifelse(DWVB.abs > 0, 1, 0)) %>%
                mutate(dwvb.f.r = dwvb.f.s * sd + m), aes(dwvb.f.r, DWVB.abs), shape = 21, height = 0.025, alpha = 0.4)+
  theme_bw(base_size = 14)+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))

##

mod = WB_models_Morisita[['bqcv.b']]
p_direction(mod)
new_dat = expand.grid(
  bqcv.f.s = seq(min(mod$data$bqcv.f.s), max(mod$data$bqcv.f.s), by = 0.1),
  Morisita.bl.z.s = 0,
  Morisita.z.s = 0,
  Year = 2022, 
  Site = 'Goe1425',
  Species = 'Bombus pascuorum'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'bb') %>% ungroup() %>%
  filter(Morisita.z > -10)
m = mean(data.bb.raw$bqcv.f)
sd = sd(data.bb.raw$bqcv.f)
x = data.frame(z = new_dat$bqcv.f, r = new_dat$bqcv.f * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('bqcv.f.s' == 'z')) 


exp_bqcv_w = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#6a7f54') +
  labs(x = expression(paste("Viral exposure to ", italic("A. mellifera"))), y = "BQCV prevalence", title = NULL)+
  annotate("text", label = "pd = 0.95", x=15,  y=0.92, hjust=0)+
  geom_jitter(data = mod$data %>% mutate(BQCV.abs = ifelse(BQCV.abs > 0, 1, 0)) %>%
                mutate(bqcv.f.r = bqcv.f.s * sd + m), aes(bqcv.f.r, BQCV.abs), shape = 21, height = 0.025, alpha = 0.4)+
  theme_bw(base_size = 14)+
  scale_x_continuous(breaks = c(11,14,17,20))+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))

##

mod = WB_models_Morisita[['abpv.b']]
p_direction(mod)
new_dat = expand.grid(
  abpv.bl.f.s = seq(min(mod$data$abpv.bl.f.s), max(mod$data$abpv.bl.f.s), by = 0.1),
  Morisita.bl.z.s = 0,
  Morisita.z.s = 0,
  Year = 2022, 
  Site = 'Goe1425',
  Species = 'Bombus pascuorum'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)

data.bb.raw = data.both %>% filter(Group == 'bb' & Species != 'Bombus lapidarius') %>% ungroup() %>%
  filter(Morisita.bl.z > -10)
m = mean(data.bb.raw$abpv.bl.f, na.rm = T)
sd = sd(data.bb.raw$abpv.bl.f, na.rm = T)
x = data.frame(z = new_dat$abpv.bl.f, r = new_dat$abpv.bl.f * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('abpv.bl.f.s' == 'z')) 

exp_abpv_w = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#d55e00') +
  labs(x = expression(paste("Viral exposure to ", italic("B. lapidarius"))), y = "ABPV prevalence", title = NULL)+
  annotate("text", label = "pd = 0.99", x=11,  y=0.92, hjust=0)+
  geom_jitter(data = mod$data %>% mutate(ABPV.abs = ifelse(ABPV.abs > 0, 1, 0)) %>%
                mutate(abpv.bl.f.r = abpv.bl.f.s * sd + m), aes(abpv.bl.f.r, ABPV.abs), shape = 21, height = 0.025, alpha = 0.4)+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme_bw(base_size = 14)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))


exp_pred = plot_grid(exp_dwvb_w, exp_bqcv_w, exp_abpv_w, align = 'hv', nrow = 3)


## RESOURCE OVERLAP 

mod = WB_models_Morisita[['bqcv.b']]
p_direction(mod)

new_dat = expand.grid(
  Morisita.z.s = seq(min(mod$data$Morisita.z.s), max(mod$data$Morisita.z.s), by = 0.1),
  Morisita.bl.z.s = 0,
  dwvb.f.s = 0,
  bqcv.f.s = 0,
  abpv.bl.f.s = 0,
  Year = 2022, 
  Site = 'Goe1425',
  Species = 'Bombus pascuorum'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)


nich_bqcv_b = ggplot(posterior_preds, aes(x = Morisita.z.s, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#6a7f54') +
  labs(x = expression(paste("Resource overlap with ", italic("A. mellifera"))), y = "BQCV prevalence", title = 'A')+
  annotate("text", label = "pd = 1.00", x=-1.8,  y=0.92, hjust=0)+
  geom_jitter(data = mod$data %>% mutate(BQCV.abs = ifelse(BQCV.abs > 0, 1, 0)), aes(Morisita.z.s, BQCV.abs), shape = 19, height = 0.025, alpha = 0.4)+ 
  theme_bw(base_size = 14)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))

##

mod = WB_models_Morisita[['bqcv.w']]
p_direction(mod)

new_dat = expand.grid(
  Morisita.z.s = seq(min(mod$data$Morisita.z.s), max(mod$data$Morisita.z.s), by = 0.1),
  Morisita.bl.z = 0,
  dwvb.f = 0,
  bqcv.f.s = 0,
  abpv.bl.f = 0,
  Year = 2022, 
  Site = 'Goe1425',
  Species = 'Bombus pascuorum'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'mu', re_formula = NA)


nich_bqcv_w = ggplot(posterior_preds, aes(x = Morisita.z.s, y = mu)) +
  stat_lineribbon(.width = c(0.95), fill = '#6a7f54') +
  labs(x = expression(paste("Resource overlap with ", italic("A. mellifera"))), y = "BQCV load (ln)", title = 'B')+
  annotate("text", label = "pd = 0.97", x=-0.4,  y=12, hjust=0)+
  geom_jitter(data = mod$data[mod$data$BQCV.abs >0,], aes(Morisita.z.s, log(BQCV.abs)), shape = 19, height = 0.025, alpha = 0.2)+
  theme_bw(base_size = 14)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))

##

mod = WB_models_Morisita[['abpv.w']]
p_direction(mod)

new_dat = expand.grid(
  Morisita.z = 0,
  Morisita.bl.z.s = seq(min(mod$data$Morisita.bl.z.s), max(mod$data$Morisita.bl.z.s), by = 0.1),
  abpv.bl.f.s = mean(mod$data$abpv.bl.f.s),
  Year = 2022, 
  Site = 'Goe1425',
  Species = 'Bombus lapidarius'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'mu', re_formula = NA)


nich_abpv_w = ggplot(posterior_preds, aes(x = Morisita.bl.z.s, y = mu)) +
  stat_lineribbon(.width = c(0.95), fill = '#d55e00') +
  labs(x = expression(paste("Resource overlap with ", italic("B. lapidarius"))), y = "ABPV load (ln)", title = 'C')+
  annotate("text", label = "pd = 0.98", x=0,  y=12, hjust=0)+
  geom_jitter(data = mod$data[mod$data$ABPV.abs >0,], aes(Morisita.bl.z.s, log(ABPV.abs)), shape = 19, height = 0.025, alpha = 0.2)+
  theme_bw(base_size = 14)+
  theme(legend.position = "none", panel.grid = element_blank(), axis.text = element_text(color = 'black'))

niche_pred = plot_grid(nich_bqcv_b, nich_bqcv_w, nich_abpv_w, nrow = 1, align = 'hv')


### CONNECTANCE

mod = WB_models_Connectance[['dwvb.b']]
p_direction(mod)

new_dat = expand.grid(
  Connectance.s = seq(min(mod$data$Connectance.s), max(mod$data$Connectance.s), by = 0.1),
  sum.ab.fl.s = 0,
  Year = 2022, 
  Site = 'Goe1425',
  Species = 'Bombus lapidarius'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)
data.bb.raw = data.both %>% filter(Group == 'bb') %>% ungroup() %>%
  filter(Morisita.z > -10)
m = mean(data.bb.raw$Connectance)
sd = sd(data.bb.raw$Connectance)
x = data.frame(z = new_dat$Connectance.s, r = new_dat$Connectance.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Connectance.s' == 'z')) 

con.a = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#547f7f') +
  labs(x = "Connectance", y = "DWV-B prevalence", title = 'A')+
  geom_jitter(data = mod$data %>% mutate(Connectance.r = Connectance.s * sd + m) , aes(Connectance.r, ifelse(DWVB.abs > 0, 1, 0)), shape = 21, height = 0.025, alpha = 0.2)+
  annotate("text", label = "pd = 0.99", x=0.38,  y=0.92, hjust=0)+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'))

##

mod = WB_models_Connectance[['dwvb.w']]
p_direction(mod)
new_dat = expand.grid(
  Connectance.s = seq(min(mod$data$Connectance.s), max(mod$data$Connectance.s), by = 0.1),
  sum.ab.fl.s = 0,
  Year = 2022, 
  Site = 'Goe1425',
  Species = 'Bombus lapidarius'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)
data.bb.raw = data.both %>% filter(Group == 'wb') %>% ungroup() %>%
  filter(Morisita.z > -10)
m = mean(data.bb.raw$Connectance)
sd = sd(data.bb.raw$Connectance)
x = data.frame(z = new_dat$Connectance.s, r = new_dat$Connectance.s * sd + m) %>% distinct(z,r)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Connectance.s' == 'z')) 

con.b = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#547f7f') +
  labs(x = "Connectance", y = "DWV-B prevalence", title = 'B')+
  geom_jitter(data = mod$data %>% mutate(Connectance.r = Connectance.s * sd + m) , aes(Connectance.r, ifelse(DWVB.abs > 0, 1, 0)), shape = 21, height = 0.025, alpha = 0.2)+
  annotate("text", label = "pd = 0.95", x=0.38,  y=0.92, hjust=0)+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'))

##

mod = WB_models_Connectance[['bqcv.b']]
p_direction(mod)

new_dat = expand.grid(
  Connectance.s = seq(min(mod$data$Connectance.s), max(mod$data$Connectance.s), by = 0.1),
  sum.ab.fl.s = 0,
  Year = 2022, 
  Site = 'Goe1425',
  Species = 'Bombus pascuorum'
)

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'hu', re_formula = NA)
data.bb.raw = data.both %>% filter(Group == 'bb') %>% ungroup() %>%
  filter(Morisita.z > -10)
m = mean(data.bb.raw$Connectance)
sd = sd(data.bb.raw$Connectance)
x = data.frame(z = new_dat$Connectance.s, r = new_dat$Connectance.s * sd + m) %>% distinct(z,r)
posterior_preds2 = left_join(posterior_preds, x, by = join_by('Connectance.s' == 'z')) 

con.c = ggplot(posterior_preds2, aes(x = r, y = 1-hu)) +
  stat_lineribbon(.width = c(0.95), fill = '#547f7f') +
  labs(x = "Connectance", y = "BQCV prevalence", title = 'D')+
  geom_jitter(data = mod$data %>% mutate(Connectance.r = Connectance.s * sd + m) , aes(Connectance.r, ifelse(BQCV.abs > 0, 1, 0)), shape = 21, height = 0.025, alpha = 0.2)+
  annotate("text", label = "pd = 0.98", x=0.38,  y=0.92, hjust=0)+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'))

posterior_preds <- epred_draws(mod, newdata = new_dat, dpar = 'mu', re_formula = NA)

posterior_preds2 = left_join(posterior_preds, x, by = join_by('Connectance.s' == 'z')) 

con.d = ggplot(posterior_preds2, aes(x = r, y = log(.epred))) +
  stat_lineribbon(.width = c(0.95), fill = '#547f7f') +
  labs(x = "Connectance", y = "BQCV load (ln)", title = 'E')+
  geom_point(data = mod$data %>% filter(BQCV.abs > 0) %>% mutate(Connectance.r = Connectance.s * sd + m) , aes(Connectance.r, log(BQCV.abs)), shape = 21,alpha = 0.2)+
  annotate("text", label = "pd = 0.94", x=0.38,  y=16, hjust=0)+
  theme_bw(base_size = 16)+
  theme(legend.position = "bottom", panel.grid = element_blank(), axis.text = element_text(color = 'black'))


con.grid = plot_grid(con.a, con.b, con.c, con.d, ncol = 2, align = 'hv')

###############################
######### FIGURES MS ##########
###############################

## FIGURE 2
r0_plots
exp_pred

plot_grid(r0_plots, exp_pred, rel_widths = c(1.9/3, 1.1/3), ncol = 2, align = 'hv', scale = 0.95)

ggsave(file = 'Data/fig/r0_pred.png', height = 10, width = 9)
ggsave(file = 'Data/fig/r0_pred.pdf', height = 10, width = 9)

## FIGURE 3
niche_pred

ggsave('Data/fig/pred_niche.pdf', width = 10, height = 3.5)
ggsave('Data/fig/pred_niche.png', width = 10, height = 3.5, bg = 'transparent')

## FIGURE 4
con.grid
prev.sep

plot_grid(con.grid, prev.sep, nrow = 1, rel_widths = c(3/4, 1/4), align = 'v')

ggsave('Data/fig/pred_conn_prev.pdf', width = 10, height = 8)
ggsave('Data/fig/pred_conn.png', width = 10, height = 8)
