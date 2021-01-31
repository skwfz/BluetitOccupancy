library(sf)
library(raster)
library(MODIS)
library(exactextractr)
library(viridis)
library(tidyverse)
# resolve namespace conflicts
select <- dplyr::select
map <- purrr::map
projection <- raster::projection

ne_land <- read_sf("data/gis-data.gpkg", "ne_land") %>%
  #Project to the native MODIS projection
  st_transform(crs = paste("+proj=sinu +lon_0=0 +x_0=0 +y_0=0",
                           "+a=6371007.181 +b=6371007.181 +units=m +no_defs"))
plot(ne_land)

ebird <- read_csv("data/ebd_greattit_junejuly_zf.csv")
tiles <- getTile(ne_land)
tiles@tile

MODIS::EarthdataLogin(usr = "severirissanen", pwd = "Hn4VqQ4hc68Sh99P")

begin_year <- format(min(ebird$observation_date), "%Y.01.01")
# end date for ebird data
end_year <- format(max(ebird$observation_date), "%Y.12.31")

tifs <- runGdal(product = "MCD12Q1", collection = "006", SDSstring = "01", 
                extent = ne_land %>% st_buffer(dist = 10000), 
                begin = begin_year, end = end_year, 
                outDirPath = "data", job = "modis",
                MODISserverOrder = "LPDAAC") %>% 
  pluck("MCD12Q1.006") %>% 
  unlist()

new_names <- format(as.Date(names(tifs)), "%Y") %>% 
  sprintf("modis_mcd12q1_umd_%s.tif", .) %>% 
  file.path(dirname(tifs), .)
file.rename(tifs, new_names)

#-------------LOAD DATA---------------
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

#------------GET LANDCOVER---------------

neighborhood_radius <- 5 * ceiling(max(res(landcover))) / 2
ebird_buff <- ebird %>% 
  distinct(year = format(observation_date, "%Y"),
           locality_id, latitude, longitude) %>%
  mutate(year_lc = if_else(as.integer(year) > max_lc_year, 
                           as.character(max_lc_year), year),
         year_lc = paste0("y", year_lc)) %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  # transform to modis projection
  st_transform(crs = projection(landcover)) %>% 
  # buffer to create neighborhood around each point
  st_buffer(dist = neighborhood_radius) %>% 
  # nest by year
  nest(data = c(year, locality_id, geometry))

# function to summarize landcover data for all checklists in a given year
calculate_pland <- function(yr, regions, lc) {
  locs <- st_set_geometry(regions, NULL)
  exact_extract(lc[[yr]], regions, progress = FALSE) %>% 
    map(~ count(., landcover = value)) %>% 
    tibble(locs, data = .) %>% 
    unnest(data)
}

# iterate over all years extracting landcover for all checklists in each
lc_extract <- ebird_buff %>% 
  mutate(pland = map2(year_lc, data, calculate_pland, lc = landcover)) %>% 
  select(pland) %>% 
  unnest(cols = pland)

pland <- lc_extract %>% 
  # calculate proporiton
  group_by(locality_id, year) %>% 
  mutate(pland = n / sum(n)) %>% 
  ungroup() %>% 
  select(-n) %>% 
  # remove NAs after tallying so pland is relative to total number of cells
  filter(!is.na(landcover))

# convert names to be more descriptive
lc_names <- tibble(landcover = 0:15,
                   lc_name = c("pland_00_water", 
                               "pland_01_evergreen_needleleaf", 
                               "pland_02_evergreen_broadleaf", 
                               "pland_03_deciduous_needleleaf", 
                               "pland_04_deciduous_broadleaf", 
                               "pland_05_mixed_forest",
                               "pland_06_closed_shrubland", 
                               "pland_07_open_shrubland", 
                               "pland_08_woody_savanna", 
                               "pland_09_savanna", 
                               "pland_10_grassland", 
                               "pland_11_wetland", 
                               "pland_12_cropland", 
                               "pland_13_urban", 
                               "pland_14_mosiac", 
                               "pland_15_barren"))
pland <- pland %>% 
  inner_join(lc_names, by = "landcover") %>% 
  arrange(landcover) %>% 
  select(-landcover)

# tranform to wide format, filling in implicit missing values with 0s%>% 
pland <- pland %>% 
  pivot_wider(names_from = lc_name, 
              values_from = pland, 
              values_fill = list(pland = 0))

# save
write_csv(pland, "data/modis_pland_location-year.csv")

pland <-read_csv("data/modis_pland_location-year.csv") %>% 
  mutate(year = as.integer(year))

#---------------VISUALIZE HABITATS--------------

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

#--------------PLOT----------------

habitats <-
  c("pland_00_water", 
    "pland_01_evergreen_needleleaf", 
    "pland_02_evergreen_broadleaf", 
    "pland_03_deciduous_needleleaf", 
    "pland_04_deciduous_broadleaf", 
    "pland_05_mixed_forest",
    "pland_06_closed_shrubland", 
    "pland_07_open_shrubland", 
    "pland_08_woody_savanna", 
    "pland_09_savanna", 
    "pland_10_grassland", 
    "pland_11_wetland", 
    "pland_12_cropland", 
    "pland_13_urban", 
    "pland_14_mosiac", 
    "pland_15_barren")

habitat = habitats[10]

habitat_cover <- pland_coords %>% 
  # convert to spatial features
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = projection(r)) %>% 
  # rasterize points
  rasterize(r, field = habitat) %>% 
  # project to albers equal-area for mapping
  projectRaster(crs = st_crs(4326)$proj4string, method = "ngb") %>% 
  # trim off empty edges of raster
  trim()

# make a map
par(mar = c(0.25, 0.25, 2, 0.25))
t <- str_glue("Proportion of {habitat}\n",
              "{max_lc_year} MODIS Landcover")
plot(habitat_cover, axes = FALSE, box = FALSE, col = viridis(10), main = t)
