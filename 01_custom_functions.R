date <- format(Sys.Date(), "%y%m%d") # for saving files

## FLOWER AND BEE DENSITY EXTRAPOLATION FUNCTION

extrapolation_function <- function(relative.data, landuse, normalization) {
  extrapolated <- relative.data %>% ungroup() %>%
    # taking mean across all sites in case of NAs (here even if the habitat was absent from a landscape, the mean will be calculated, but then it will be multiplied by zero anyway)
    mutate(across(any_of(c("Flower_fieBS2", "Flower_fieBS12", 
                           "Flower_fieBS11", "Other_AUM", 
                           "Fallow", "CropBV1", 
                           "semi_natur", "Grassy_str")), ~ ifelse(is.na(.x), mean(.x, na.rm = T), .x))) %>%
    left_join(landuse, by = 'Site') %>% 
    # extrapolating (.x are the flower estimates, .y is the area in hectares)
    mutate(Org.ex = CropBV1.x * CropBV1.y/normalization,
           Fl_BS11_ex = Flower_fieBS11.x * Flower_fieBS11.y/normalization,
           Fl_BS12_ex = Flower_fieBS12.x * Flower_fieBS12.y/normalization,
           Fl_BS2_ex = Flower_fieBS2.x * Flower_fieBS2.y/normalization,
           Fallow_ex = Fallow.x * Fallow.y/normalization,
           Other_AUM_ex = Other_AUM.x * Other_AUM.y/normalization,
           Grassy_str_ex = Grassy_str.x * Grassy_str.y/normalization,
           semi_natur_ex = semi_natur.x * semi_natur.y/normalization) %>%
    # summing areas
    mutate(sum = Org.ex + Fl_BS11_ex + Fl_BS12_ex + Fl_BS2_ex + Fallow_ex + Other_AUM_ex + semi_natur_ex + Grassy_str_ex) %>%
    # mutate(EX_per_agr = sum/(Crop + CropBV1.y + Flower_fieBS2.y + Flower_fieBS12.y + Flower_fieBS11.y + Fallow.y +
    #                                 Other_AUM.y + semi_natur.y + Grassy_str.y)) %>%
    select(Site, sum)
  return(extrapolated)
}

### PATHOGEN DATA FILTERING
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

## MODEL SUMMARY TABLE

custom_summary <- function(mod) {
  sum <- summary(mod)[['fixed']] %>% select(Estimate, Est.Error, Bulk_ESS, Tail_ESS) %>% 
    cbind(ci(as_draws_array(mod), ci = 0.95, method = 'ETI') %>% 
            filter(grepl('b_', Parameter)) %>% filter(!grepl('prior', Parameter)) %>% left_join(p_direction(mod)[,1:2], by = 'Parameter') %>% 
            select(-c(CI, Parameter)) %>% rename(CI_low_95 = CI_low, CI_high_95 = CI_high)) %>% 
    mutate(across(everything(), ~ round(., 2))) %>%
    relocate(c(CI_low_95, CI_high_95), .after = Est.Error) %>% mutate(Response = mod[["formula"]][["resp"]], Predictor = rownames(.)) %>% 
    relocate(c(Response, Predictor), .before = Estimate)
  rownames(sum) <- NULL
  return(sum)
}

## MODEL FIT CHECKS

pp_hurdle <- function(mod) {
  pred <- posterior_predict(mod)
  y1 <- mod[["data"]][[1]]
  ppc_dens_overlay(y = log1p(y1), 
                   yrep = log1p(pred[1:length(y1),]))
}

posterior_correlation <- function(fit = fit, cor = 0.5) {
  post <- posterior_samples(fit) %>% 
    select(starts_with("b_"))
  
  cor_mat <- cor(post)
  
  cor_df <- cor_mat %>%
    as.data.frame() %>%
    mutate(var1 = row.names(.)) %>%
    pivot_longer(
      cols = -var1,
      names_to = "var2",
      values_to = "correlation"
    )
  
  cor_df_clean <- cor_df %>%
    filter(var1 < var2) 
  
  
  cor_df_clean %>% filter(correlation > cor |
                            correlation < -1*cor)
  
}

## REVERSING SCALING FOR PREDICTIORS

posterior_pred_unscaled <- function(mod, newdata, dpar, variable, bee_group) {
  unscaled_var <- sub(".s", "",variable)
  posterior_preds <- epred_draws(mod, newdata = newdata, dpar = dpar, re_formula = NA)
  data.raw = data.both %>% filter(Group == bee_group) %>% ungroup()
  m = mean(data.raw[[unscaled_var]])
  sd = sd(data.raw[[unscaled_var]])
  x = data.frame(z = new_dat[, variable], r = new_dat[, variable] * sd + m) %>% distinct(z,r)
  
  posterior_preds2 = left_join(posterior_preds, x, by = join_by(!!variable == 'z'))
  posterior_preds2$m = m
  posterior_preds2$sd = sd
  return(posterior_preds2)
}

######################
#### R0 FUNCTIONS ####
######################

dat_into_matrix = function(abundance, true.prev, data, virus, sociality = 2) {
  names(data) = tolower(names(data))
  names(true.prev) = tolower(names(true.prev))
  names(abundance) = tolower(names(abundance))
  t.prev <- true.prev %>% 
    group_by(site, species) %>% 
    rename(vir = paste0(virus, '.true')) %>% 
    summarise(prev = mean(vir), n = n()) %>% 
    filter(n > 2) %>% 
    select(-n)
  load <- data %>% 
    rename(vir = virus, vir.load = paste0(virus, '.abs')) %>% 
    filter(vir > 0) %>% 
    group_by(site, species) %>% 
    summarise(load = mean(vir.load))
  
  mat.dat <- full_join(abundance %>% 
                         select(site, species, bee_abundance), t.prev, by = c('site','species')) %>% 
    full_join(load, by = c('site','species')) %>% 
    group_by(species) %>%
    mutate(prev = ifelse(is.na(prev), mean(prev, na.rm = T), prev),
           load = ifelse(is.na(load), mean(load, na.rm = T), load)) %>%
    ungroup() %>% 
    filter(bee_abundance > 0) %>%
    mutate(load = ifelse(is.na(load), 100, load)) %>% 
    drop_na(prev) %>% 
    left_join(data %>% distinct(species, social), by = 'species') %>%
    mutate(social = as.numeric(ifelse(social == 0, 1, sociality)))
  return(mat.dat)
}

niche_fun = function(mat.dat.list, morisita) {
  inter.niche <- list()
  inter.niche2 <- list()
  for (i in 1:length(mat.dat.list)) {
    inter.niche[[i]] <- as.data.frame(morisita[[i]]) %>% rownames_to_column() %>% mutate(rowname = sub(' agg.', '', rowname)) %>% 
      column_to_rownames()
    diag(inter.niche[[i]]) = 1
    colnames(inter.niche[[i]]) <- rownames(inter.niche[[i]])
    inter.niche2[[i]] <- inter.niche[[i]] %>% filter(row.names(inter.niche[[i]]) %in% mat.dat.list[[i]]$species) 
    inter.niche2[[i]] <- inter.niche2[[i]] %>% select(rownames(inter.niche2[[i]]))
    excluded_species <- setdiff(colnames(inter.niche[[i]]), colnames(inter.niche2[[i]])) 
    message(glue("Species removed from site {i}: {paste0(excluded_species, collapse = ',')}"))
  }
  return(inter.niche2)
}

r0_func_fac <- function(data.matrix.raw, inter.niche.ex, species_out = 'Apis mellifera') {
  result <- matrix(NA, nrow = nrow(data.matrix.raw), ncol = 1)
  r0 <- matrix(NA, nrow = nrow(data.matrix.raw), ncol = nrow(data.matrix.raw))
  transmission <- matrix(NA, nrow = nrow(data.matrix.raw), ncol = nrow(data.matrix.raw))
  data.matrix <- as.matrix(data.matrix.raw[,3:6])
  
  # Calculate niche overlap for each pair of species
  for (j in seq_len(nrow(r0))) {
    for (i in seq_len(nrow(r0))) {
      transmission[j,i] <- inter.niche.ex[j,i]
      transmission[j,j] <- data.matrix[j, 4]
      colnames(transmission) = colnames(inter.niche.ex)
      rownames(transmission) = rownames(inter.niche.ex)
      r0[j,i] <- log10(data.matrix[i, 3])/log10(data.matrix[j, 3]) * data.matrix[i, 1]/data.matrix[j, 1] * data.matrix[i, 2] * 
        (transmission[j,i]/data.matrix[j, 4])
      #r0[j,j] <- 1
    }
    result[j,1] <- 1 / (((1 - data.matrix[j, 2])/data.matrix[j, 2]) * sum(r0[j,], na.rm = T))
    rownames(result) <- rownames(inter.niche.ex)
    colnames(result) <- 'r0'
    #result_com = list(r0 = result, transmission = transmission)
  }
  # remove main host
  data.matrix2 = data.matrix.raw %>% filter(!species %in% species_out)
  transmission2 = as.data.frame(transmission) %>% filter(row.names(transmission) %in% data.matrix2$species)
  transmission2 = as.data.frame(transmission2) %>% select(rownames(transmission2))
  r02 = as.data.frame(result)  %>% filter(row.names(result) %in% data.matrix2$species)
  pred.prev <- matrix(NA, nrow = nrow(data.matrix2), ncol = 1)
  mat <- matrix(NA, nrow = nrow(data.matrix2), ncol = nrow(data.matrix2))
  data.matrix2 <- as.matrix(data.matrix2[,3:6])
  ## run the code solving for prevalence
  for (o in seq_len(nrow(mat))) {
    for (p in seq_len(nrow(mat))) {
      mat[o,p] <- log10(data.matrix2[p, 3])/log10(data.matrix2[o, 3]) * data.matrix2[p, 1]/data.matrix2[o, 1] * data.matrix2[p, 2] * 
        (transmission2[o,p]/data.matrix2[o, 4])
    }
    pred.prev[o,1] <- r02[o,1] * sum(mat[o,], na.rm = T) / (1 + r02[o,1] * sum(mat[o,], na.rm = T))
    rownames(pred.prev) <- rownames(transmission2)
    #colnames(pred.prev) <- 'prev_no_hb'
  }
  result_com = as.data.frame(result) %>% rownames_to_column() %>% left_join(as.data.frame(pred.prev) %>% rownames_to_column(), by = 'rowname') %>% rename(Prev.no.mh = V1) %>%
    left_join(data.matrix.raw %>% select(species, prev), by = join_by('rowname' == 'species'))
  return(result_com)
}

r0_wrapper <- function(abundance, true.prev, data, virus, is_subset, morisita, species_out, sociality = 2) {
  
  mat.dat = dat_into_matrix(abundance, true.prev, data, virus, sociality)
  
  if (is_subset) {
    
    mat.dat = mat.dat %>% filter(site %in% subset)
    
  }
  
  mat.dat.list <- split(mat.dat, mat.dat$site)
  
  inter.niche = niche_fun(mat.dat.list, morisita)
  
  for (i in 1:length(mat.dat.list)) {
    mat.dat.list[[i]] <- mat.dat.list[[i]] %>% filter(species %in% rownames(inter.niche[[i]]))
  }
  
  r0 <- list()
  for (s in 1:length(mat.dat.list)){
    r0[[s]]  <- r0_func_fac(mat.dat.list[[s]], inter.niche[[s]], species_out = species_out)
  }
  
  return(r0)
}