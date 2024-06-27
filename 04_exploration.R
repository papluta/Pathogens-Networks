source('03_data_processing.R')
source('05_plotting_functions.R')


# correlation between extrapolated honeybees and flower cover standardized to agricultural area
plot(data$FL_per_agr ~ data$HB_per_agr)
plot(data$connectance ~ data$niche.overlap.HL)
plot(data$connectance ~ data$weighted.nested)
plot(data$weighted.nested ~ data$niche.overlap.HL)
dotchart(data$HB_per_agr)
dotchart(data$FL_per_agr)
dotchart(log(data$FL_per_agr+1))
dotchart(log(data$sum.fl+1))
dotchart(data$sum.fl)


## separate data frames for each bee group
data.hb <- data %>% filter(Group == 'hb') %>% ungroup() %>% 
  mutate(across(SH:sbv.f, ~ scale(.)[,1]))

data.bb <- data %>% filter(Group == 'bb' | Group == 'bp') %>% ungroup() %>% #filter(Species == 'Bombus lapidarius') %>%
  mutate(across(SH:sbv.f, ~ scale(.)[,1]))


data.wb <- data %>% filter(Group == 'wb') %>% 
  mutate(Genus = gsub(' .*', '', Species)) %>%
  # filtering out erroneous barcodes
  filter(Genus != '') %>% filter(Genus != ('Sipha')) %>% filter(Genus != ('Trypoxylon')) %>% filter(Genus != ('Orisarma')) %>% 
  filter(Genus != ('Lindenius')) %>% filter(Genus != ('Oedogonium')) %>% filter(Genus != ('Orasema')) %>% filter(Genus != ('Megalocoleus')) %>% 
  filter(Genus != ('Bombus')) %>% ungroup() %>%
  mutate(across(SH:sbv.f, ~ scale(.)[,1]))

b <- data.wb %>% group_by(Species) %>% summarise(S_n = n()) %>% filter(S_n == 1) # 236 Andrena, 55 Lasioglossum
print(data.wb %>% group_by(Site) %>% summarise(n = n()), n = Inf) # Nor 174 only 6 samples, Wm1316 7

data.wb <- data.wb %>% left_join(b, by = 'Species') %>% filter(is.na(S_n))

data.site <- data %>% filter(Species == 'Apis mellifera') %>% distinct(Site, sum.hb, HB_per_agr, sum.fl, FL_per_agr, Density, 
                                                                       dwvb.hb.agr, bqcv.hb.agr, abpv.hb.agr, sbv.hb.agr) %>%
  ungroup() %>% mutate(across(4:5, ~scale(.)[,1]))

# separate models for each virus and bee group
# virus.hb.agr = honey bee prevalence * extrapolated honey bee abundance
# not running honey bee models, because the infection force is based on their viral prevalence

library(glmmTMB)
library(car)
library(performance)
library(effects)
library(ggeffects)
library(DHARMa)
library(MuMIn)
options(na.action = 'na.fail')

# hb

data.hb.dwvb <- data.hb %>% filter(DWVB.abs > 0)
data.hb.bqcv <- data.hb %>% filter(BQCV.abs > 0)
data.hb.abpv <- data.hb %>% filter(ABPV.abs > 0)
data.hb.sbv <- data.hb %>% filter(SBV.abs > 0)
dotchart(log10(data.hb.abpv$ABPV.abs), labels = data.hb.dwvb$Site)
dotchart(data.hb.dwvb$sbv.hb.agr, labels = data.hb.dwvb$Site)
dotchart(log10(data.hb.abpv$DWVB.abs), labels = data.hb.abpv$Site)

mod <- glmmTMB(dwvb ~ HB_per_agr + FL_per_agr + (1 |Site), family = 'binomial', data = data.hb)
#mod <- glmmTMB(bqcv ~ HB_per_agr + FL_per_agr + (1 |Site), family = 'binomial', data = data.hb) # 100% prevalence
mod <- glmmTMB(abpv ~ HB_per_agr + FL_per_agr + (1 |Site), family = 'binomial', data = data.hb) 
mod <- glmmTMB(sbv ~ HB_per_agr + FL_per_agr + (1 |Site), family = 'binomial',data = data.hb)
summary(mod)
check_collinearity(mod)
plot(simulateResiduals(mod))
plot(allEffects(mod))
plot(ggpredict(mod, terms = c('HB_per_agr', 'FL_per_agr[-1,1,3]'), back_transform = F))

mod <- glmmTMB(log10(DWVB.abs) ~ HB_per_agr +  FL_per_agr  +  (1 |Site), data = data.hb.dwvb) # 
mod <- glmmTMB(log10(BQCV.abs) ~ HB_per_agr +  FL_per_agr  +   (1 |Site), data = data.hb.bqcv) # !
mod <- glmmTMB(log10(ABPV.abs) ~ HB_per_agr + FL_per_agr +   (1|Site), data = data.hb.abpv) # !
mod <- glmmTMB(log10(SBV.abs) ~ HB_per_agr  + FL_per_agr +  (1 |Site), data = data.hb.sbv)
summary(mod)
check_collinearity(mod)
plot(simulateResiduals(mod))
plot(allEffects(mod))

plot(ggpredict(mod, terms = c('HB_per_agr'), back_transform = F))+
  geom_point(data = data.hb.bqcv, aes(HB_per_agr, log10(BQCV.abs)), alpha = 0.3)
plot(ggpredict(mod, terms = c('FL_per_agr'), back_transform = F))+
  geom_point(data = data.hb.abpv, aes(FL_per_agr, log10(ABPV.abs)), alpha = 0.3)

## bb
data.bb.dwvb <- data.bb %>% filter(DWVB.abs > 0)
data.bb.bqcv <- data.bb %>% filter(BQCV.abs > 0)
data.bb.abpv <- data.bb %>% filter(ABPV.abs > 0) %>% filter(ABPV.abs < 10^10)
data.bb.sbv <- data.bb %>% filter(SBV.abs > 0)
dotchart(log10(data.bb.abpv$ABPV.abs+1), labels = data.bb.abpv$Site)
dotchart(data.bb.abpv$abpv.hb.agr, labels = data.bb.abpv$Site)

mod <- glmmTMB(dwvb ~ HB_per_agr * dwvb.hb + FL_per_agr + Species + (1 |Site), family = 'binomial', data = data.bb)
mod <- glmmTMB(bqcv ~ HB_per_agr + FL_per_agr + Species + (1 |Site), family = 'binomial', data = data.bb)
mod <- glmmTMB(abpv ~ HB_per_agr * abpv.hb + FL_per_agr + Species + (1 |Site), family = 'binomial', data = data.bb) #!
mod <- glmmTMB(sbv ~ HB_per_agr * sbv.hb + FL_per_agr + Species + (1 |Site), family = 'binomial',data = data.bb)
summary(mod)
check_collinearity(mod)
plot(simulateResiduals(mod))
plot(allEffects(mod))
plot(ggpredict(mod, terms = 'FL_per_agr', back_transform = F))+
  geom_point(data = data.bb %>% group_by(FL_per_agr, abpv) %>% summarise(s = length(FL_per_agr)), aes(FL_per_agr, abpv, size = s), alpha = 0.5)


mod <- glmmTMB(log10(DWVB.abs) ~ HB_per_agr * dwvb.hb +  FL_per_agr  + Species + (1 |Site), data = data.bb.dwvb) # residuals
mod <- glmmTMB(log10(BQCV.abs) ~ HB_per_agr + FL_per_agr  +  Species + (1 |Site), data = data.bb.bqcv)
mod <- glmmTMB(log10(ABPV.abs) ~ HB_per_agr + abpv.hb +  Species + (1|Site), data = data.bb.abpv) 
mod <- glmmTMB(log10(SBV.abs) ~ HB_per_agr * sbv.hb  + FL_per_agr + Species + (1 |Site), data = data.bb.sbv) #!

summary(mod)
check_collinearity(mod)
plot(simulateResiduals(mod))
testDispersion(mod)
plot(allEffects(mod))
plot(ggpredict(mod, terms = c( 'HB_per_agr', 'dwvb.hb'), back_transform = F))
p <- ggpredict(mod, terms = c('HB_per_agr', 'dwvb.hb'), back_transform = F)

p %>% ggplot(aes(x, predicted))+
  geom_line(aes(col = group))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.3)+
  geom_point(data = data.bb.abpv, aes(abpv.hb, log10(ABPV.abs), col = cut(FL_per_agr, breaks = 2)), alpha = 0.3)

dotchart(data.bb.abpv$FL_per_agr)
# wb

data.wb.dwvb <- data.wb %>% filter(DWVB.abs > 0)
data.wb.bqcv <- data.wb %>% filter(BQCV.abs > 0)
data.wb.abpv <- data.wb %>% filter(ABPV.abs > 0)
data.wb.sbv <- data.wb %>% filter(SBV.abs > 0)
dotchart(log10(data.wb.abpv$ABPV.abs), labels = data.wb.abpv$Site)
plot(data.wb.abpv$FL_per_agr ~ log10(data.wb.abpv$ABPV.abs))
plot(density(log10(data.wb.abpv$ABPV.abs)))
dotchart(data.wb.abpv$FL_per_agr, labels = data.wb.abpv$Site)

mod <- glmmTMB(dwvb ~ HB_per_agr * dwvb.hb   + FL_per_agr + (1|Species) + (1 |Site), family = 'binomial', data = data.wb) # 
mod <- glmmTMB(bqcv ~ bqcv.hb.agr  + FL_per_agr+  (1|Species) + (1 |Site), family = 'binomial', data = data.wb)
mod <- glmmTMB(abpv ~ HB_per_agr * abpv.hb  + FL_per_agr +  (1|Species) + (1 |Site), family = 'binomial', data = data.wb) 
mod <- glmmTMB(sbv ~ HB_per_agr * sbv.hb  + FL_per_agr + (1|Species) + (1 |Site), family = 'binomial', data = data.wb)
summary(mod)
check_collinearity(mod)
testDispersion(mod)
plot(simulateResiduals(mod))
plot(allEffects(mod))
plot(ggpredict(mod, terms = c('abpv.f', 'sum.fl'), back_transform = F))


mod <- glmmTMB(log10(DWVB.abs) ~ dwvb.hb.agr + FL_per_agr + (1|Species) + (1 |Site), data = data.wb.dwvb)
mod <- glmmTMB(log10(BQCV.abs) ~ bqcv.hb.agr + FL_per_agr + (1|Species) + (1 |Site), data = data.wb.bqcv)
mod <- glmmTMB(log10(ABPV.abs) ~ abpv.hb.agr + FL_per_agr + (1|Species) + (1 |Site), data = data.wb.abpv) ## residuals!
mod <- glmmTMB(log10(SBV.abs) ~ sbv.hb.agr + FL_per_agr + (1|Species) + (1 |Site), data = data.wb.sbv)
summary(mod)
check_collinearity(mod)
DHARMa::testDispersion(mod)
plot(simulateResiduals(mod))
plot(allEffects(mod))
plot(ggpredict(mod, terms = c('FL_per_agr', 'HB_per_agr'), back_transform = F))

mod <- glmmTMB(abpv ~ abpv.f + FL_per_agr + Species + (1 |Site), family = 'binomial', data = data.bb)




### bayesian
library(brms)
library(tidybayes)
library(tidyr)
library(ggplot2)
library(rstan)
library(bayesplot)
library(bayestestR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


dwvb.b <- brm(bf(DWVB.abs ~ dwvb.f + FL_per_agr + Species + (1 |Site),
                  hu ~ dwvb.f + FL_per_agr + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.w <- brm(bf(DWVB.abs ~ dwvb.f + FL_per_agr + (1 |Site) + (1|Species),
                  hu ~ dwvb.f + FL_per_agr +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

bqcv.b <- brm(bf(BQCV.abs ~ bqcv.f + FL_per_agr +  Species + (1 |Site),
                  hu ~ bqcv.f + FL_per_agr + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.w <- brm(bf(BQCV.abs ~ bqcv.f + FL_per_agr + (1 |Site) + (1|Species),
                  hu ~ bqcv.f + FL_per_agr +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

abpv.b <- brm(bf(ABPV.abs ~ abpv.f * sum.fl +  Species +(1 |Site),
                  hu ~ abpv.f * sum.fl + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.w <- brm(bf(ABPV.abs ~ abpv.f + FL_per_agr + (1 |Site) + (1|Species),
                  hu ~ abpv.f + FL_per_agr +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

sbv.b <- brm(bf(SBV.abs ~ sbv.f + FL_per_agr +  Species + (1 |Site),
                  hu ~ sbv.f + FL_per_agr + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)
sbv.w <- brm(bf(SBV.abs ~ sbv.f + FL_per_agr + (1 |Site) + (1|Species),
                  hu ~ sbv.f + FL_per_agr +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)

FL_FORCE_models <- list(dwvb.b, dwvb.w, 
                            bqcv.b, bqcv.w, 
                            abpv.b, abpv.w, 
                            sbv.b, sbv.w)

names(FL_FORCE_models) <- c('dwvb.b', 'dwvb.w',
                                'bqcv.b', 'bqcv.w', 
                                'abpv.b', 'abpv.w', 
                                'sbv.b', 'sbv.w')

save(FL_FORCE_models, file = 'Data/Models/240627FL_FORCE_models.RData')

## plotting and diagnosing

abpv.b2 <- brm(log10(ABPV.abs) ~ abpv.f * sum.fl +Species + (1|Site), family = gaussian, data = data.bb.abpv, prior = 
                c(prior(normal(0,5), class = 'b')), sample_prior = TRUE)

pp_hurdle(abpv.b)
pp_check(abpv.b2)
summary(abpv.b)
p_direction(abpv.b2)
conditional_effects(abpv.b2)
#load('FLOWER_FORCE_models.RData')

cs.all <- lapply(FL_FORCE_models, custom_summary) %>% bind_rows(.id = 'id')
rownames(cs.all) <- NULL

# to see which relationships were likely:
cs.all %>% filter(pd > 0.9) %>% filter(Predictor != 'Intercept' & Predictor != 'hu_Intercept' & Predictor != 'SpeciesBombuspascuorum' & 
                                         Predictor != 'hu_SpeciesBombuspascuorum') %>% arrange(Predictor)

# POSTERIOR CHECK
pp.all <- lapply(FL_FORCE_models, pp_hurdle) #check for each model
pp.all[['abpv.b']]
pp.all[['abpv.w']]

## PLOTTING EFFECTS
conditional_effects(FL_FORCE_models[['abpv.w']], effects = 'abpv.f:FL_per_agr', dpar = 'mu')

mu_force <- lapply(CON_FORCE_models, function(x) conditional_effects(x, dpar = 'mu')[[1]])
mu_flower <- lapply(CON_FORCE_models, function(x) conditional_effects(x, dpar = 'mu')[[2]])
hu_force <- lapply(CON_FORCE_models, function(x) conditional_effects(x, dpar = 'hu')[[1]])
hu_flower <- lapply(CON_FORCE_models, function(x) conditional_effects(x, dpar = 'hu')[[2]])

conditional_effects(abpv.w, dpar = 'mu')


plots <- lapply(mu_force, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Honey bee force', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plots <- lapply(mu_flower, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Flower cover', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plots <- lapply(hu_force, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Honey bee force', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plots <- lapply(hu_flower, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Flower cover', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)


plot(density(data.site$HB_per_agr))

bee.ab <- brm(log(HB_per_agr+1) ~ FL_per_agr * Density, family = gaussian, data = data.site, sample_prior = T)
summary(bee.ab)
pp_check(bee.ab, ndraws = 100)
p_direction(bee.ab)
conditional_effects(bee.ab)

dwvb.force <- brm(log(dwvb.f+1) ~ FL_per_agr * Density, family = gaussian, data = data.site, sample_prior = T)
summary(dwvb.force)
pp_check(dwvb.force, ndraws = 100)
p_direction(dwvb.force)
conditional_effects(dwvb.force)

library(DHARMa)

m <- lm(log(HB_per_agr+1) ~ FL_per_agr * Density, data.site)
summary(m)
plot(simulateResiduals(m))


dwvb.h <- brm(bf(DWVB.abs ~ HB_per_agr + FL_per_agr + (1 |Site),
                 hu ~ HB_per_agr + FL_per_agr +  (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

bqcv.h <- brm(bf(BQCV.abs ~ HB_per_agr + FL_per_agr +  (1 |Site),
                 hu ~ HB_per_agr + FL_per_agr + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.h <- brm(bf(ABPV.abs ~ HB_per_agr + FL_per_agr + (1 |Site),
                 hu ~ HB_per_agr + FL_per_agr +  (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

sbv.h <- brm(bf(SBV.abs ~ HB_per_agr + FL_per_agr +  (1 |Site),
                hu ~ HB_per_agr + FL_per_agr + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)

hb.mods <- list(dwvb.h = dwvb.h, bqcv.h = bqcv.h, abpv.h = abpv.h, sbv.h = sbv.h)
pairs(bqcv.h)

cs.all <- lapply(hb.mods, custom_summary) %>% bind_rows(.id = 'id')
rownames(cs.all) <- NULL

# to see which relationships were likely:
cs.all %>% filter(pd > 0.9) %>% filter(Predictor != 'Intercept' & Predictor != 'hu_Intercept') %>% arrange(Predictor)

# POSTERIOR CHECK
pp.all <- lapply(hb.mods, pp_hurdle)
pp.all[['dwvb.h']]


## PLOTTING EFFECTS
mu_sum <- lapply(hb.mods, function(x) conditional_effects(x, dpar = 'mu')[[1]])
mu_flower <- lapply(hb.mods, function(x) conditional_effects(x, dpar = 'mu')[[2]])
hu_sum <- lapply(hb.mods, function(x) conditional_effects(x, dpar = 'hu')[[1]])
hu_flower <- lapply(hb.mods, function(x) conditional_effects(x, dpar = 'hu')[[2]])


plots <- lapply(mu_sum, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Honey bee abundance', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 1, nrow = 4)

plots <- lapply(mu_flower, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Flower cover', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 1, nrow = 4)


######################### OTHER STUFF, NOT SORTED OUT ############################

library(bipartite)
ex <- data %>% group_by(Species) %>% summarise(n = n()) %>% filter(n >2)

data1 <- data[c(1:7, 13)] %>% group_by(Density, Species) %>% summarise(dwvb = sum(dwvb), bqcv = sum(bqcv), sbv = sum(sbv), abpv = sum(abpv)) %>%
  mutate(Density = factor(Density, levels = c('0', '1'))) %>% filter(Species %in% ex$Species)



net.data <- split(data1, data1$Density)
net.data[[1]] <- net.data[[1]][-1]
net.data[[2]] <- net.data[[2]][-1]
low <- net.data$`0`
high <- net.data$`1`
low <- low[-1]
high <- high[-1]
low2 <- as.matrix(low)
high2 <- as.matrix(high)
rownames(low2) <- net.data[[1]]$Species
rownames(high2) <- net.data[[2]]$Species
a <- betalinkr(webs2array(low2, high2), partitioning="commondenom")
a

low3 <- t(low2)
plotweb(low2)
#low <- low %>% mutate(across(everything(), function(x) ifelse(x > 0, 1, 0)))


### WORKS
col <- adjustcolor(c('#7f5454','#6a7f54', '#547f7f', '#6a547f'), alpha = 0.6)
low$Species <- net.data[[1]]$Species
g.data.low <- low %>% pivot_longer(cols = dwvb:abpv, names_to = 'Virus', values_to = 'Yes') %>% 
  mutate(e.col = rep(col, times = 17)) %>% filter(Yes != 0) %>% rename(weight = Yes) %>%
  #mutate(Species = sub('(.).* ', '. ', Species)) 
  mutate(Species = sub('Andrena', "A.", Species)) %>% mutate(Species = sub('Apis', "A.", Species)) %>% 
  mutate(Species = sub('Bombus', "B.", Species)) %>% 
  mutate(Species = sub('Lasioglossum', "L.", Species)) %>% 
  mutate(Species = sub('Melitta', "M.", Species)) %>% 
  mutate(Species = sub('Eucera', "E.", Species)) %>% mutate(Species = sub('Colletes', "C.", Species)) %>%
  mutate(Virus = recode_factor(as.factor(Virus), 'dwvb' = 'DWV-B', bqcv = 'BQCV', sbv = 'SBV', abpv = 'ABPV'))
#deg <- low %>% pivot_longer(cols = dwvb:abpv, names_to = 'Virus', values_to = 'Yes') %>% filter(Yes != 0) %>% select(Yes)
g.data.low <- g.data.low[,c(2,1,3,4)]
#colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
g <- graph.data.frame(g.data.low, directed = T)

V(g)$type <- bipartite.mapping(g)$type
V(g)$shape <- ifelse(V(g)$type == TRUE, "circle", "circle")
V(g)$color <- c('#547f7f', '#7f5454','#6a547f','#6a7f54',  rep('#ecde51', times = 17))
V(g)$frame.color <- c('#547f7f', '#7f5454','#6a547f','#6a7f54', rep('#ecde51', times = 17))
V(g)$size <- degree(g)+25
E(g)$color <- E(g)$e.col

E(g)$width <- E(g)$weight*0.1
plot(g, layout = layout_in_circle(g), edge.arrow.size = 0.5, vertex.label.color = 'black', vertex.label.cex = 0.8, vertex.frame.color = NULL)

## high dens
high$Species <- net.data[[2]]$Species
g.data.high <- high %>% pivot_longer(cols = dwvb:abpv, names_to = 'Virus', values_to = 'Yes') %>% 
  mutate(e.col = rep(col, times = 15)) %>% filter(Yes != 0) %>% rename(weight = Yes) %>%
  #mutate(Species = sub('(.).* ', '. ', Species)) 
  mutate(Species = sub('Andrena', "A.", Species)) %>% mutate(Species = sub('Apis', "A.", Species)) %>% 
  mutate(Species = sub('Bombus', "B.", Species)) %>% 
  mutate(Species = sub('Lasioglossum', "L.", Species)) %>% 
  mutate(Species = sub('Melitta', "M.", Species)) %>% mutate(Species = sub('Seladonia', "S.", Species)) %>% 
  mutate(Species = sub('Eucera', "E.", Species)) %>% mutate(Species = sub('Colletes', "C.", Species)) %>%
  mutate(Virus = recode_factor(as.factor(Virus), 'dwvb' = 'DWV-B', bqcv = 'BQCV', sbv = 'SBV', abpv = 'ABPV'))
#deg <- high %>% pivot_longer(cols = dwvb:abpv, names_to = 'Virus', values_to = 'Yes') %>% filter(Yes != 0) %>% select(Yes)
g.data.high <- g.data.high[,c(2,1,3,4)]
#colrs <- adjustcolor( c("gray50", "tomato", "gold", "yelhighgreen"), alpha=.6)
g <- graph.data.frame(g.data.high, directed = T)

V(g)$type <- bipartite.mapping(g)$type
#V(g)$shape <- ifelse(V(g)$type == TRUE, "circle", "circle")
V(g)$color <- ifelse(V(g)$type  == TRUE, "#ecde51", "#ec5a51")
V(g)$color <- c('#547f7f', '#7f5454','#6a547f','#6a7f54',  rep('#ecde51', times = 15))
V(g)$frame.color <- c('#547f7f', '#7f5454','#6a547f','#6a7f54', rep('#ecde51', times = 15))
V(g)$size <- degree(g)+25
E(g)$color <- E(g)$e.col

E(g)$width <- E(g)$weight*0.1
plot(g, layout = layout_in_circle(g), edge.arrow.size = 0.5, vertex.label.color = 'black', vertex.label.cex = 0.8, vertex.frame.color = NULL)


dwvb.h <- brm(bf(DWVB.abs ~ totalFC_h + totalHB_h + (1 |Site),
                 hu ~ totalFC_h + totalHB_h + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.b <- brm(bf(DWVB.abs ~ totalFC_h + dwvb.hb.K + Species + (1 |Site),
                 hu ~ totalFC_h + dwvb.hb.K + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.w <- brm(bf(DWVB.abs ~ totalFC_h + dwvb.hb.K + (1 |Site) + (1|Species),
                 hu ~ totalFC_h + dwvb.hb.K +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.h <- brm(bf(BQCV.abs ~ totalFC_h + totalHB_h +  (1 |Site),
                 hu ~ totalFC_h + totalHB_h + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.b <- brm(bf(BQCV.abs ~ totalFC_h + bqcv.hb.K +  Species + (1 |Site),
                 hu ~ totalFC_h + bqcv.hb.K + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
bqcv.w <- brm(bf(BQCV.abs ~ totalFC_h + bqcv.hb.K + (1 |Site) + (1|Species),
                 hu ~ totalFC_h + bqcv.hb.K +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

abpv.h <- brm(bf(ABPV.abs ~ totalFC_h + totalHB_h +  (1 |Site),
                 hu ~ totalFC_h + totalHB_h + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.b <- brm(bf(ABPV.abs ~ totalFC_h + abpv.hb.K +  Species + (1 |Site),
                 hu ~ totalFC_h + abpv.hb.K + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
abpv.w <- brm(bf(ABPV.abs ~ totalFC_h + abpv.hb.K +(1 |Site) + (1|Species),
                 hu ~ totalFC_h + abpv.hb.K +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)

sbv.h <- brm(bf(SBV.abs ~ totalFC_h + totalHB_h +  (1 |Site),
                hu ~ totalFC_h + totalHB_h + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)
sbv.b <- brm(bf(SBV.abs ~ totalFC_h + sbv.hb.K +  Species + (1 |Site),
                hu ~ totalFC_h + sbv.hb.K + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)
sbv.w <- brm(bf(SBV.abs ~ totalFC_h + sbv.hb.K + (1 |Site) + (1|Species),
                hu ~ totalFC_h + sbv.hb.K +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
               c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                 prior(normal(0,5), class = 'b')), sample_prior = TRUE)

Kat_fl_force_models <- list(dwvb.h, dwvb.b, dwvb.w, 
                        bqcv.h, bqcv.b, bqcv.w, 
                        abpv.h, abpv.b, abpv.w, 
                        sbv.h, sbv.b, sbv.w)

names(Kat_fl_force_models) <- c('dwvb.h', 'dwvb.b', 'dwvb.w', 
                            'bqcv.h', 'bqcv.b', 'bqcv.w', 
                            'abpv.h', 'abpv.b', 'abpv.w', 
                            'sbv.h', 'sbv.b', 'sbv.w')

pairs(bqcv.b)

save(FL_dens_models, file = 'C:/Users/patry/OneDrive/PhD/tralala/Pathogens-Networks/Data/Models/FL_dens_mod.RData')
save(FL_abu_models, file = 'C:/Users/patry/OneDrive/PhD/tralala/Pathogens-Networks/Data/Models/FL_abu_mod.RData')
save(FL_force_models, file = 'C:/Users/patry/OneDrive/PhD/tralala/Pathogens-Networks/Data/Models/FL_force_mod.RData')
save(FL_forceNI_models, file = 'C:/Users/patry/OneDrive/PhD/tralala/Pathogens-Networks/Data/Models/FL_force_mod.RData')
save(nest_forceNI_models, file = 'C:/Users/patry/OneDrive/PhD/tralala/Pathogens-Networks/Data/Models/nest_forceNI_models.RData')
save(con_forceNI_models, file = 'C:/Users/patry/OneDrive/PhD/tralala/Pathogens-Networks/Data/Models/con_forceNI_models.RData')
save(FL_forceNI_models2, file = 'C:/Users/patry/OneDrive/PhD/tralala/Pathogens-Networks/Data/Models/FL_force_mod2.RData')

cs.all <- lapply(Kat_fl_force_models, custom_summary) %>% bind_rows(.id = 'id')
rownames(cs.all) <- NULL
cs.all %>% filter(pd > .9) %>% filter(Predictor != 'Intercept' & Predictor != 'hu_Intercept' & 
                                        Predictor != 'SpeciesBombuspascuorum' & Predictor != 'hu_SpeciesBombuspascuorum') %>% arrange(pd)

cs.all_Pat <- lapply(FL_forceNI_models2, custom_summary) %>% bind_rows(.id = 'id')
rownames(cs.all_Pat) <- NULL
hm <- cs.all_Pat %>% filter(pd > .9) %>% filter(Predictor != 'Intercept' & Predictor != 'hu_Intercept' & Predictor != 'SpeciesBombuspascuorum' & 
                                            Predictor != 'hu_SpeciesBombuspascuorum') %>% arrange(pd)

conditional_effects(FL_forceNI_models2[['sbv.b']], dpar = 'hu')

SpeciesBombuspascuorum
plot(data$totalFC_h ~ data$dwvb.hb.K)

Kat_only_wb <- list(dwvb.b, dwvb.w, 
                            bqcv.b, bqcv.w, 
                            abpv.b, abpv.w, 
                            sbv.b, sbv.w)
e <- conditional_effects(dwvb.w, dpar = 'mu')

mu_force <- lapply(Kat_only_wb, function(x) conditional_effects(x, dpar = 'mu')[[2]])
mu_flower <- lapply(Kat_only_wb, function(x) conditional_effects(x, dpar = 'mu')[[1]])


plots <- lapply(mu_force, function(x) x %>% 
               ggplot(aes(effect1__, estimate__))+
               #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
               geom_line(linewidth = 1)+
               geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
               #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
               #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
               labs(x = 'Honey bee force', y = 'Viral load')+
               theme_bw(base_size = 16)+
               #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
               #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
               theme(panel.grid = element_blank(),  legend.position = 'none',
                     plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                     plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plots <- lapply(mu_flower, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Flower cover', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)


hu_force <- lapply(Kat_only_wb, function(x) conditional_effects(x, dpar = 'hu')[[2]])
hu_flower <- lapply(Kat_only_wb, function(x) conditional_effects(x, dpar = 'hu')[[1]])


plots <- lapply(hu_force, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Honey bee force', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

plots <- lapply(hu_flower, function(x) x %>% 
                  ggplot(aes(effect1__, estimate__))+
                  #geom_point(data = dat, aes(x = x, y = log(y), col = Flower2), alpha = 0.2, shape = 21) +
                  geom_line(linewidth = 1)+
                  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`, alpha = 0.3))+
                  #annotation_custom(grob = grobTree(textGrob(paste0("pd = ", pd), x=0.65,  y=0.94, hjust=0,
                  #                                           gp=gpar(col="black", fontsize=10, fontface=NULL))))+
                  labs(x = 'Flower cover', y = 'Viral load')+
                  theme_bw(base_size = 16)+
                  #scale_fill_manual(values = c('#f1c40f','#00806A', '#0a81f1'))+
                  #scale_color_manual(values = c('#f1c40f','#00806A', '#0a81f1', '#f1c40f','#00806A', '#0a81f1'))+
                  theme(panel.grid = element_blank(),  legend.position = 'none',
                        plot.margin = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
                        plot.title = element_text(size = 16, hjust = 0.5)))
plot_grid(plotlist = plots, ncol = 2, nrow = 4)

pp.all <- lapply(FL_forceNI_models2, pp_hurdle)
pp.all[[12]] #+ coord_cartesian(xlim = c(0, 6))

conditional_effects(FL_forceNI_models2[['sbv.b']], effects = 'sum.fl', dpar = 'mu')
conditional_effects(FL_forceNI_models2[['dwvb.b']], effects = 'dwvb.hb.f', dpar = 'mu')







dens_plot_mu(SH_dens_models[[4]], effects = 'sum.fl:Density') 

hm <- lapply(FL_abu_models, function(x) force_plot(x, effects = 'sum.hb:sum.fl'))
plot_grid(plotlist = hm, align = 'vh',nrow = 4, ncol = 3)

conditional_effects(FL_forceNI_models[[2]], effects = 'dwvb.hb.f', dpar = 'mu')
conditional_effects(FL_forceNI_models[[6]], effects = 'bqcv.hb.f', dpar = 'mu')
conditional_effects(FL_forceNI_models[[9]], effects = 'abpv.hb.f', dpar = 'mu')
conditional_effects(FL_forceNI_models[[12]], effects = 'sbv.hb.f', dpar = 'mu')
conditional_effects(FL_forceNI_models[[12]], effects = 'weighted', dpar = 'mu')
force_plot_hu(FL_force_models[[6]], effects = 'bqcv.hb.f:sum.fl')
force_plot_hu(FL_force_models[[9]], effects = 'abpv.hb.f:sum.fl')
force_plot_hu(FL_force_models[[11]], effects = 'sbv.hb.f:sum.fl')
conditional_effects(nest_forceNI_models[[12]], effects = 'weighted.nested', dpar = 'mu')
conditional_effects(nest_forceNI_models[[1]], effects = 'weighted.nested', dpar = 'hu')

hm_hu <- lapply(FL_abu_models, function(x) force_plot_hu(x, effects = 'sum.hb:sum.fl'))
plot_grid(plotlist = hm_hu, align = 'vh',nrow = 4, ncol = 3)

dens_plot_hu(SH_dens_models[[11]], effects = 'SH') 

conditional_effects(SH_dens_models[[11]], dpar = 'hu')

nm <- posterior_density(SH_dens_models, 'b_Density1')
plot(nm)

plot_forest(SH_dens_models[[1]], 'SH:Density', dpar = 'hu')

plot_forest(SH_dens_models, effect = 'Density')

plot <- lapply(SH_dens_models, custom_summary) 
plot2 <- lapply(plot, function(x) filter(x, Predictor == 'SH')) %>% bind_rows(.id = 'mod') %>%
  ggplot(aes(x = Estimate, y = mod))+
  geom_vline(xintercept = 0, color = 'grey', linewidth = 0.8)+
  geom_linerange(aes(xmin = CI_low_90, xmax = CI_high_90), linewidth = 1, color = 'black')+
  #geom_linerange(aes(xmin = CI_l90, xmax = CI_h90), linewidth = 2.1, color = '#a40b0b')+
  geom_point(shape = 21, size = 4, color = 'black', fill = '#dc3a3a')+
  theme_bw(base_size = 16, base_line_size = 16/44)+
  #theme(axis.text.y = element_blank())+
  labs(x = NULL, y = NULL)
plot2

## net

library(lme4)
library(glmmTMB)
library(DHARMa)
library(car)
library(performance)
library(effects)
library(splines)

net <- data %>% distinct(Site, connectance, niche.overlap.HL, `weighted.nested`, Density, SH, sum.fl, sum.hb, totalFC_h, FL_per_agr, totalHB, HB_per_agr)

m <- lm(connectance ~ Density + FL_per_agr, net)
m <- lm(weighted.nested ~ Density + FL_per_agr, net)
m <- lm(log(niche.overlap.HL+1) ~ Density + FL_per_agr, net)
m <- lm(log(sum.hb+1) ~ Density + FL_per_agr, net)
m <- lm(log(totalHB+1) ~ Density + FL_per_agr, net)

plot(net$totalHB ~ net$FL_per_agr)
plot(net$sum.hb ~ net$sum.fl)
summary(m)
DHARMa::testDispersion(m)
plot(so <- simulateResiduals(m))

plot(allEffects(m))
hist(log(net$connectance))

dwvb. <- glmer(dwvb ~ sum.fl + sum.hb + (1 |Site), family = 'binomial', data = data.hb)
dwvb. <- glmer(dwvb ~ sum.fl + dwvb.hb.f + Species + (1 |Site), family = 'binomial', data = data.bb)
dwvb. <- glmer(abpv ~ sum.fl + abpv.hb.f + Species + (1 |Site), family = 'binomial', data = data.bb)
summary(dwvb.)
vif(dwvb.)
DHARMa::testDispersion(dwvb.)
plot(so <- simulateResiduals(dwvb.))
plot(allEffects(dwvb.))


dwvb.b.f <- glmmTMB(log(DWVB.abs) ~ sum.fl + dwvb.hb.f + Species + (1 |Site), data = data.bb[data.bb$DWVB.abs > 0,])
dwvb.b.f <- glmmTMB(log(BQCV.abs) ~ sum.fl  + bqcv.hb.f + Species + (1 |Site), data = data.bb[data.bb$BQCV.abs > 0,])
dwvb.b.f <- glmmTMB(log(ABPV.abs) ~ sum.fl + abpv.hb.f + Species + (1 |Site), data = data.bb[data.bb$ABPV.abs > 0,])
dwvb.b.f <- glmmTMB(log(SBV.abs) ~ sum.fl + Species + (1 |Site), data = data.bb[data.bb$SBV.abs > 0,])
summary(dwvb.b.f)
DHARMa::testDispersion(dwvb.b.f)
plot(so <- simulateResiduals(dwvb.b.f))            
check_collinearity(dwvb.b.f)
plot(allEffects(dwvb.b.f))

dwvb.h <- brm(bf(DWVB.abs ~ weighted.nested + Density + (1 |Site),
                 hu ~ weighted.nested + Density + (1 |Site)), family = hurdle_lognormal(), data = data.hb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.b <- brm(bf(DWVB.abs ~ weighted.nested + Density + Species + (1 |Site),
                 hu ~ weighted.nested + Density + Species + (1 |Site)), family = hurdle_lognormal(), data = data.bb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
dwvb.w <- brm(bf(DWVB.abs ~ weighted.nested + Density + (1 |Site) + (1|Species),
                 hu ~ weighted.nested + Density +  (1 |Site) + (1|Species)), family = hurdle_lognormal(), data = data.wb, prior = 
                c(prior(normal(0,5), class = 'b', dpar = 'hu'),
                  prior(normal(0,5), class = 'b')), sample_prior = TRUE)
