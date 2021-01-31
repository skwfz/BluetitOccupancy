library(auk)
library(lubridate)
library(sf)
library(dggridR)
library(unmarked)
library(raster)
library(ebirdst)
library(MuMIn)
library(AICcmodavg)
library(fields)
library(tidyverse)

# resolve namespace conflicts
select <- dplyr::select
projection <- raster::projection

# setup output directory for saved results
if (!dir.exists("output")) {
  dir.create("output")
}

ebird <- read_csv("data/ebd_greattit_junejuly_zf.csv") %>% 
  mutate(year = year(observation_date),
         # occupancy modeling requires an integer response
         species_observed = as.integer(species_observed))

# modis land cover covariates
habitat <- read_csv("data/modis_pland_location-year.csv") %>% 
  mutate(year = as.integer(year))

# combine ebird and habitat data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

# load gis data for making maps
map_proj <- st_crs(4326)
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()

# filter prior to creating occupancy model data
ebird_filtered <- filter(ebird_habitat, 
                         number_observers <= 5,
                         year == max(year))

occ <- filter_repeat_visits(ebird_filtered, 
                            min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id", "observer_id"))

nrow(ebird)
nrow(occ)
n_distinct(occ$site)

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

# format for unmarked
occ_wide <- format_unmarked_occu(occ, 
                                 site_id = "site", 
                                 response = "species_observed",
                                 site_covs = append(c("n_observations", 
                                               "latitude", "longitude"), habitats),
                                 obs_covs = append(c("time_observations_started", 
                                              "duration_minutes", 
                                              "effort_distance_km", 
                                              "number_observers", 
                                              "protocol_type"), habitats))
