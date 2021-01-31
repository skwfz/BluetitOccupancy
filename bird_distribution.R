library(sf)
library(raster)
library(MODIS)
library(exactextractr)
library(viridis)
library(tidyverse)
library(boot)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

#-------------LOAD MAP DATA---------------
ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>%
  #Project to the native MODIS projection
  st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                           "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))


#-------------LOAD MODIS DATA---------------
landcover <- list.files("data/modis", "^modis_mcd12q1_umd", 
                        full.names = TRUE) %>% 
  stack()
landcover <- names(landcover) %>% 
  str_extract("(?<=modis_mcd12q1_umd_)[0-9]{4}") %>% 
  paste0("y", .) %>% 
  setNames(landcover, .)
landcover

max_lc_year <- names(landcover) %>% 
  str_extract("[0-9]{4}") %>% 
  as.integer() %>% 
  max()

#-----------GET PREDICTION SURFACE----------
pland <-read_csv("data/modis_pland_location-year.csv") %>% 
  mutate(year = as.integer(year))

agg_factor <- round(2 * neighborhood_radius / res(landcover))
r <- raster(landcover) %>% 
  aggregate(agg_factor) 
r <- ne_land %>% 
  st_transform(crs = projection(r)) %>% 
  rasterize(r, field = 1) %>% 
  # remove any empty cells at edges
  trim()
r <- writeRaster(r, filename = "data/prediction-surface.tif", overwrite = TRUE)

# get cell centers and create neighborhoods
r_centers <- rasterToPoints(r, spatial = TRUE) %>% 
  st_as_sf() %>% 
  transmute(id = row_number())
r_cells <- st_buffer(r_centers, dist = neighborhood_radius)

# extract landcover values within neighborhoods, only needed most recent year
lc_extract_pred <- landcover[[paste0("y", max_lc_year)]] %>% 
  exact_extract(r_cells, progress = FALSE) %>% 
  map(~ count(., landcover = value)) %>% 
  tibble(id = r_cells$id, data = .) %>% 
  unnest(data)

# calculate the percent for each landcover class
pland_pred <- lc_extract_pred %>%
  group_by(id) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

# convert names to be more descriptive
pland_pred <- pland_pred %>% 
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>%
  select(-landcover)

# tranform to wide format, filling in implicit missing values with 0s
pland_pred <- pland_pred %>% 
  pivot_wider(names_from = lc_name, 
              values_from = pland, 
              values_fill = list(pland = 0)) %>% 
  mutate(year = max_lc_year) %>% 
  select(id, year, everything())

# join in coordinates
pland_coords <- st_transform(r_centers, crs = 4326) %>% 
  st_coordinates() %>% 
  as.data.frame() %>% 
  cbind(id = r_centers$id, .) %>% 
  rename(longitude = X, latitude = Y) %>% 
  inner_join(pland_pred, by = "id")

#---------------LOAD MODEL-------------------
load("data/model4fit")
draws_model4 <- rstan::extract(stan_fit_4, permuted = T)

pland_coords
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

pland_coords[,5:17]

draws_model4$habitat_coef
occupancy_means <- 0*(1:dim(pland_coords)[1])
occupancy_sds <- 0*(1:dim(pland_coords)[1])

n_samples <- dim(draws_model4$habitat_coef)[1]
dim(pland_coords)

for(i in 1:dim(pland_coords)[1]){
  probs <- inv.logit(rowMeans(draws_model4$habitat_coef * pland_coords[i,5:17][rep.int(1, n_samples),]))
  occupancy_means[i] <- mean(probs)
  occupancy_sds[i] <- sd(probs)
  if(i %% 1000 == 0){
    print(i)
  }
}

occupancy_tibble <- tibble(id=1:dim(pland_coords)[1], occupancy_means=occupancy_means, occupancy_sds=occupancy_sds)

expected_occupancy_cover <- pland_coords %>% inner_join(occupancy_tibble, by="id") %>%
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r, field = occupancy_means) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs(4326)$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

uncertainty_occupancy_cover <- pland_coords %>% inner_join(occupancy_tibble, by="id") %>%
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r, field = occupancy_sds) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs(4326)$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

# make a map
par(mfrow=c(1,2),mar = c(0.25, 0.25, 2, 0.25))

t <- str_glue("Estimated occupancy of\n",
              "the Eurasian Blue Tit in 2020, June-July")
plot(expected_occupancy_cover, axes = FALSE, box = FALSE, col = viridis(10), main = t)


t <- str_glue("Occupancy uncertainty of\n",
              "the Eurasian Blue Tit in 2020, June-July")
p2 <- plot(uncertainty_occupancy_cover, axes = FALSE, box = FALSE, col = viridis(10), main = t)

