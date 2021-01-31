library(lubridate)
library(sf)
library(dggridR)
library(raster)
library(ebirdst)
library(fields)
library(tidyverse)
library(rstan)

# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

ebird <- read_csv("data/ebd_bluetit_junejuly_zf.csv") %>% 
  mutate(year = year(observation_date),
         # occupancy modeling requires an integer response
         species_observed = as.integer(species_observed))
# modis land cover covariates
habitat <- read_csv("data/modis_pland_location-year.csv") %>% 
  mutate(year = as.integer(year))

# combine ebird and habitat data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

# filter prior to creating occupancy model data
ebird_filtered <- filter(ebird_habitat, 
                         number_observers <= 5,
                         year == max(year))

time_freq <- ebird_filtered %>%
  mutate(tod_bins = cut(time_observations_started, 
                        breaks = 0:24, 
                        labels = 0:23,
                        include.lowest = TRUE),
         tod_bins = as.numeric(as.character(tod_bins))) %>% 
  group_by(tod_bins) %>% 
  summarise(n_checklists = n(),
            n_detected = sum(species_observed),
            det_freq = mean(species_observed))
plot(time_freq$tod_bins, time_freq$det_freq )

# load gis data for making maps
map_proj <- st_crs(4326)
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

N_localities <- n_distinct(ebird_filtered$locality_id)
loc_ids <- tibble(locality_id = unique(ebird_filtered$locality_id), locality_simpleid = 1:N_localities)

ebird_filtered <- inner_join(ebird_filtered, loc_ids, by = c("locality_id"))
ebird_filtered$locality_simpleid

habitats <-
  c("pland_00_water", 
    "pland_01_evergreen_needleleaf", 
    "pland_03_deciduous_needleleaf", 
    "pland_04_deciduous_broadleaf", 
    "pland_05_mixed_forest", 
    "pland_07_open_shrubland", 
    "pland_08_woody_savanna", 
    "pland_09_savanna", 
    "pland_10_grassland", 
    "pland_11_wetland", 
    "pland_12_cropland", 
    "pland_13_urban",
    "pland_15_barren")

habitat_prop <- distinct(inner_join(loc_ids, ebird_filtered[,append(c("locality_id"),habitats)], by=c("locality_id")))
habitat_prop

forest_prop <- ebird_filtered %>%
  select(c("pland_01_evergreen_needleleaf", 
           "pland_03_deciduous_needleleaf", 
           "pland_04_deciduous_broadleaf", 
           "pland_05_mixed_forest")) %>%
  mutate(sum = rowSums(.[1:4])) %>%
  select("sum")

ppcheck_indices <- 1:dim(ebird_filtered)[1]#sample(1:dim(ebird_filtered)[1],300)
ppcheck_indices_habprop <- ebird_filtered$locality_simpleid#ebird_filtered[ppcheck_indices,]$locality_simpleid

stan_fit_1 <- stan("model.stan", data=list(N=dim(ebird_filtered)[1],
                                         observed=ebird_filtered$species_observed,
                                         N_habitats=length(habitats),
                                         N_localities=dim(loc_ids)[1],
                                         habitat_prop=habitat_prop[,habitats],
                                         localities=ebird_filtered$locality_simpleid,
                                         duration_minutes=ebird_filtered$duration_minutes,
                                         distance_traveled=ebird_filtered$effort_distance_km,
                                         forest_prop=forest_prop$sum,
                                         check_N=length(ppcheck_indices),
                                         check_indices_habprop=ppcheck_indices_habprop,
                                         check_indices_others=ppcheck_indices),
                 iter=2000, chains=4)
save("stan_fit_1",file="data/model1fit")
load("data/model1fit")

stan_fit_3 <- stan("model3.stan", data=list(N=dim(ebird_filtered)[1],
                                           observed=ebird_filtered$species_observed,
                                           N_habitats=length(habitats),
                                           N_localities=dim(loc_ids)[1],
                                           habitat_prop=habitat_prop[,habitats],
                                           localities=ebird_filtered$locality_simpleid,
                                           duration_minutes=ebird_filtered$duration_minutes,
                                           distance_traveled=ebird_filtered$effort_distance_km,
                                           forest_prop=forest_prop$sum,
                                           start_time=ebird_filtered$time_observations_started,
                                           check_N=length(ppcheck_indices),
                                           check_indices_habprop=ppcheck_indices_habprop,
                                           check_indices_others=ppcheck_indices),
                   iter=2000, chains=4)
save("stan_fit_3",file="data/model3fit")
load("data/model3fit")

stan_fit_4 <- stan("model4.stan", data=list(N=dim(ebird_filtered)[1],
                                            observed=ebird_filtered$species_observed,
                                            N_habitats=length(habitats),
                                            N_localities=dim(loc_ids)[1],
                                            habitat_prop=habitat_prop[,habitats],
                                            localities=ebird_filtered$locality_simpleid,
                                            duration_minutes=ebird_filtered$duration_minutes,
                                            distance_traveled=ebird_filtered$effort_distance_km,
                                            forest_prop=forest_prop$sum,
                                            start_time=ebird_filtered$time_observations_started,
                                            check_N=length(ppcheck_indices),
                                            check_indices_habprop=ppcheck_indices_habprop,
                                            check_indices_others=ppcheck_indices),
                   iter=2000, chains=4)
save("stan_fit_4",file="data/model4fit")
load("data/model4fit")

draws_model1 <- rstan::extract(stan_fit_1, permuted = T)
draws_model3 <- rstan::extract(stan_fit_3, permuted = T)
draws_model4 <- rstan::extract(stan_fit_4, permuted = T)


#----------------------DO THE LOGREG POSTERIOR PREDICTIVE THING FROM THE LECTURE SLIDES---------------

plot_ppcheck <- function(draws_model,modeltitle){
  ebird_ppcheck <- tibble(ebird_filtered[ppcheck_indices,], estimated_prob = rowMeans(t(draws_model$observed_check)))
  breaks <- seq(0,0.8,0.1)
  labels <- (breaks[-1] + breaks[-(length(breaks))])/2
  binomial_ppcheck <- ebird_ppcheck %>%
    mutate(bin=cut(ebird_ppcheck$estimated_prob, breaks = breaks, 
                   labels = labels,
                   include.lowest = TRUE),
           bin = as.numeric(as.character(bin))) %>% 
    group_by(bin) %>% 
    summarise(n_checklists = n(),
              n_detected = sum(species_observed),
              det_freq = mean(species_observed),
              sd = sqrt(det_freq*(1-det_freq)/n_checklists))
  
  plot <- ggplot(binomial_ppcheck) + 
    geom_pointrange(aes(x=bin, y=det_freq, ymin=det_freq-2*sd,ymax=det_freq+2*sd)) + 
    geom_point(data=ebird_ppcheck,aes(x=estimated_prob, y=species_observed)) + 
    geom_line(data=tibble(prob=seq(0,1,0.01),prop=seq(0,1,0.01)),aes(x=prob,y=prop)) + 
    xlab("Estimated probability") + ylab("") + labs(title=modeltitle)
  return(plot)
}
library(patchwork)
g <- plot_ppcheck(draws_model1, "Model 1")+ ylab("Bird observed/not observed and binned frequency") + 
  plot_ppcheck(draws_model3, "Model 2") + plot_ppcheck(draws_model4, "Model 3")
g
ggsave(g, file="posterior_predictive_1.png" , width=10, height=4)


#---------------------DISTANCE TRAVELED-----------------------
library(gridExtra)
plot_ppcheck_dist <- function(draws_model,modeltitle){
# summarize data by 500m bins
  breaks <- seq(0, 5, by = 0.5)
  labels <- breaks[-length(breaks)] + diff(breaks) / 2
  ebird_dist <- ebird_filtered[ppcheck_indices,] %>% 
    mutate(dist_bins = cut(effort_distance_km, 
                           breaks = breaks, 
                           labels = labels,
                           include.lowest = TRUE),
           dist_bins = as.numeric(as.character(dist_bins))) %>% 
    group_by(dist_bins) %>% 
    summarise(n_checklists = n(),
              n_detected = sum(species_observed),
              det_freq = mean(species_observed),
              sd = sqrt(det_freq*(1-det_freq)/n_checklists))
  
  #summarize posterior predictive draws
  pp_dist <- tibble(effort_distance_km=ebird_filtered[ppcheck_indices,]$effort_distance_km,t(draws_model$observed_check)) %>%
    mutate(mean_probs = rowMeans(.),
           dist_bins = cut(effort_distance_km, 
                           breaks = breaks, 
                           labels = labels,
                           include.lowest = TRUE),
           dist_bins = as.numeric(as.character(dist_bins))) %>% 
    group_by(dist_bins) %>% 
    summarise(det_freq = mean(mean_probs))
  
  #frequency of detection for posterior predictive draws
  g_dist_freq <- ggplot(pp_dist) +
    aes(x = dist_bins, y = det_freq) +
    geom_line(col="red") +
    geom_point(col="red", size=2.4) +
    scale_x_continuous(breaks = 0:5) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Distance travelled (km)",
         y = "% checklists with detections",
         title = modeltitle) + 
    geom_line(data=ebird_dist,aes(x = dist_bins, y = det_freq)) +
    geom_point(data=ebird_dist,aes(x = dist_bins, y = det_freq)) + 
    geom_pointrange(data=ebird_dist, aes(x=dist_bins, y=det_freq, ymin=det_freq-2*sd,ymax=det_freq+2*sd))
  return(g_dist_freq)
}

g <-plot_ppcheck_dist(draws_model1,"Model 1") + plot_ppcheck_dist(draws_model3,"Model 2") + plot_ppcheck_dist(draws_model4,"Model 3")
g
ggsave(g, file="posterior_predictive_dist.png" , width=10, height=4)

#-----------------------START TIME PLOTS------------------------
plot_ppcheck_startime <- function(draws_model,modeltitle){
  ebird_ppcheck <- tibble(ebird_filtered[ppcheck_indices,], estimated_prob = rowMeans(t(draws_model$observed_check)))
  breaks <- 0:24
  labels <- breaks[-length(breaks)] + diff(breaks) / 2
  starttime_ppcheck <- ebird_ppcheck %>%
    mutate(tod_bins = cut(time_observations_started, 
                          breaks = breaks, 
                          labels = labels,
                          include.lowest = TRUE),
           tod_bins = as.numeric(as.character(tod_bins))) %>%
    group_by(tod_bins) %>% 
    summarise(n_checklists = n(),
              n_detected = sum(species_observed),
              det_freq = mean(species_observed),
              estimated_prob = mean(estimated_prob))
  g_return <- ggplot(starttime_ppcheck) + 
    geom_line(aes(x=tod_bins,y=det_freq)) + geom_point(aes(x=tod_bins,y=det_freq)) + 
    geom_line(aes(x=tod_bins,y=estimated_prob),color='red') + geom_point(aes(x=tod_bins,y=estimated_prob),color='red') + 
    labs(x = "Start time (h)",
         y = "% checklists with detections",
         title = modeltitle) 
  return(g_return)
}
g <-plot_ppcheck_startime(draws_model1,"Model 1") + plot_ppcheck_startime(draws_model3,"Model 2") + plot_ppcheck_startime(draws_model4,"Model 3")
g
ggsave(g, file="posterior_predictive_time.png" , width=10, height=4)

#PSIS-LOO CHECKS
probs_model1 <- draws_model1$occupied_prob*draws_model1$detected_prob
probs_model3 <- draws_model3$occupied_prob*draws_model3$detected_prob
probs_model4 <- draws_model4$occupied_prob*draws_model4$detected_prob
loo(log(probs_model1))
loo(log(probs_model3))
loo(log(probs_model4))

#ACCURACIES
mean(ebird_filtered$species_observed == (colMeans(probs_model1) > 0.5))
mean(ebird_filtered$species_observed == (colMeans(probs_model3) > 0.5))
mean(ebird_filtered$species_observed == (colMeans(probs_model4) > 0.5))

#----------------------------PRIOR SENSITIVITY CHECKS---------------------------------

stan_fit_4_2 <- stan("model4.stan", data=list(N=dim(ebird_filtered)[1],
                                            observed=ebird_filtered$species_observed,
                                            N_habitats=length(habitats),
                                            N_localities=dim(loc_ids)[1],
                                            habitat_prop=habitat_prop[,habitats],
                                            localities=ebird_filtered$locality_simpleid,
                                            duration_minutes=ebird_filtered$duration_minutes,
                                            distance_traveled=ebird_filtered$effort_distance_km,
                                            forest_prop=forest_prop$sum,
                                            start_time=ebird_filtered$time_observations_started,
                                            check_N=length(ppcheck_indices),
                                            check_indices_habprop=ppcheck_indices_habprop,
                                            check_indices_others=ppcheck_indices),
                   iter=2000, chains=4)
save("stan_fit_4_2",file="data/model4fit_2")
load("data/model4fit_2")

stan_fit_4_3 <- stan("model4.stan", data=list(N=dim(ebird_filtered)[1],
                                              observed=ebird_filtered$species_observed,
                                              N_habitats=length(habitats),
                                              N_localities=dim(loc_ids)[1],
                                              habitat_prop=habitat_prop[,habitats],
                                              localities=ebird_filtered$locality_simpleid,
                                              duration_minutes=ebird_filtered$duration_minutes,
                                              distance_traveled=ebird_filtered$effort_distance_km,
                                              forest_prop=forest_prop$sum,
                                              start_time=ebird_filtered$time_observations_started,
                                              check_N=length(ppcheck_indices),
                                              check_indices_habprop=ppcheck_indices_habprop,
                                              check_indices_others=ppcheck_indices),
                     iter=2000, chains=4)
save("stan_fit_4_3",file="data/model4fit_3")
load("data/model4fit_3")

draws_model4_2 <- rstan::extract(stan_fit_4_2, permuted = T)
draws_model4_3 <- rstan::extract(stan_fit_4_3, permuted = T)

g <- plot_ppcheck(draws_model4, "Model 3, normal")+ ylab("Bird observed/not observed and binned frequency") + 
  plot_ppcheck(draws_model4_2, "Model 3, tight prior") + plot_ppcheck(draws_model4_3, "Model 3, wide prior")
g
ggsave(g, file="posterior_predictive_sensitivity.png" , width=10, height=4)


#----------------------------RESULT ANALYSIS-----------------------------

habitat_names <- c("Water","Evergreen needleleaf","Deciduous needleleaf", "Deciduous broadleaf", "Mixed forest", "Open shrubland", "Woody savanna", 
                   "Savanna", "Grassland", "Wetland", "Cropland", "Urban", "Barren")

habitat_coef_draws <- data.frame(draws_model4$habitat_coef)
names(habitat_coef_draws) <- habitat_names
occ_plot <- ggplot(stack(habitat_coef_draws), aes(x=ind,y=values)) +
  geom_boxplot()+ggtitle("Occupancy coefficients") + 
  theme(axis.text.x = element_text(face = "bold", color = "#993333", size = 8, angle = 45,hjust = 1)) + 
  labs(x="",y="")

detection_coef_draws <- tibble(duration_coef = draws_model$duration_coef, 
                               distance_coef=draws_model$distance_coef, 
                               forest_prop_coef=draws_model$forest_prop_coef) %>%
                        mutate(duration_coef = duration_coef * sd(ebird_filtered$duration_minutes),
                               distance_coef = distance_coef * sd(ebird_filtered$effort_distance_km),
                               forest_prop_coef = forest_prop_coef * sd(forest_prop$sum))
det_plot <- ggplot(pivot_longer(detection_coef_draws,cols=c("duration_coef", "distance_coef", "forest_prop_coef")), aes(x=name,y=value)) + 
  geom_boxplot()+ggtitle("Normalized detectability coefficients") + 
  theme(axis.text.x = element_text(face = "bold", color = "#993333", size = 8, angle = 45,hjust = 1)) + 
  labs(x="",y="")

g <- occ_plot+det_plot
ggsave(g, file="coefs.png" , width=10, height=4)
