###########################
### 5. Submarine Cables ###
###########################

# Clear environment
rm(list = ls())

# Calculate start time of code (determine how long it takes to complete all code)
start <- Sys.time()

#####################################
#####################################

# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(docxtractr,
               dplyr,
               elsa,
               fasterize,
               fs,
               ggplot2,
               janitor,
               ncf,
               paletteer,
               pdftools,
               plyr,
               purrr,
               raster,
               RColorBrewer,
               reshape2,
               rgdal,
               rgeoda,
               rgeos,
               rmapshaper,
               rnaturalearth, # use devtools::install_github("ropenscilabs/rnaturalearth") if packages does not install properly
               sf,
               sp,
               stringr,
               terra, # is replacing the raster package
               tidyr)

# Inspect packages and system
sessionInfo()

#####################################
#####################################

# Set directories
## Define data directory (as this is an R Project, pathnames are simplified)
### Input directories
submarine_cable_dir <- "data/a_raw_data/SubmarineCable/NOAAChartedSubmarineCables.gdb"
submarine_cable_area_dir <- "data/a_raw_data/SubmarineCableArea/SubmarineCableArea.gpkg"

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Analysis directories
industry_operations_submodel <- "data/c_submodel_data/oregon_industry_operations_submodel.gpkg"

#### Intermediate directories
submarine_cable_gpkg <- "data/b_intermediate_data/oregon_submarine_cable.gpkg"

#####################################
#####################################

# Inspect geodatabase and geopackage to find layer names
## Submarine Cable
sf::st_layers(dsn = submarine_cable_dir,
              do_count = TRUE)

## Submarine Cable Area
sf::st_layers(dsn = submarine_cable_area_dir,
              do_count = TRUE)

## Study Area
sf::st_layers(dsn = study_area_gpkg,
              do_count = TRUE)

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

## layer names
layer_500 <- "submarine_cable500m"
layer_1000 <- "submarine_cable1000m"

## designate date
date <- format(Sys.time(), "%Y%m%d")

#####################################
#####################################

# Load data
## Oregon call areas
oregon_call_areas <- sf::st_read(dsn = wind_area_gpkg,
                                 layer = paste(sf::st_layers(dsn = wind_area_gpkg,
                                                             do_count = TRUE)))

## Oregon hex areas (original data)
oregon_hex <- sf::st_read(dsn = study_area_gpkg,
                          layer = paste(sf::st_layers(dsn = study_area_gpkg,
                                                      do_count = TRUE)[[1]][4]))

#####################################

## NOAA Charted submarine cable data (source: https://marinecadastre.gov/downloads/data/mc/SubmarineCable.zip)
### Metadata: https://www.fisheries.noaa.gov/inport/item/57238
submarine_cables_noaa <- sf::st_read(dsn = submarine_cable_dir,
                                     layer = paste(sf::st_layers(dsn = submarine_cable_dir,
                                                                 do_count = TRUE)))

### ***Note: due to a change in how sf casts data, instead of being able to cast the entire dataset, the
###          multicurve is no longer supported by GEOS, and breaks the st_is_empty() [as of 14 July 2023]
class(st_geometry(submarine_cables_noaa[1171,]))
submarine_cables_noaa$Shape[1171] <- st_cast(submarine_cables_noaa$Shape[1171], "MULTILINESTRING")
class(st_geometry(submarine_cables_noaa[1171,]))

submarine_cables_noaa <- submarine_cables_noaa %>%
  # make sure all geometries are valid
  sf::st_make_valid() %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910") # EPSG 26910 (https://epsg.io/26910)

#####################################
#####################################

# Submarine cables in Oregon call areas
oregon_submarine_cables <- submarine_cables_noaa %>%
  # obtain only submarine cables in the study area
  rmapshaper::ms_clip(target = .,
                      clip = oregon_call_areas) %>%
  # create field called "layer" and fill with "submarine cables" for summary
  dplyr::mutate(layer = "submarine cables")
  
#####################################

# Add setback distances
## 500-meter setback
oregon_submarine_cable500 <- oregon_submarine_cables %>%
  #  add a setback (buffer) distance of 500 meters
  sf::st_buffer(dist = 500) %>%
  # group all features by the "layer"
  dplyr::group_by(layer) %>%
  # summarise data to obtain single feature
  dplyr::summarise()

## 1000-meter setback (501 - 1000m setback)
oregon_submarine_cable1000 <- oregon_submarine_cables %>%
  #  add a setback (buffer) distance of 1000 meters
  sf::st_buffer(dist = 1000) %>%
  # remove areas between 0 - 500 meters
  rmapshaper::ms_erase(oregon_submarine_cable500) %>%
  # group all features by the "layer"
  dplyr::group_by(layer) %>%
  # summarise data to obtain single feature
  dplyr::summarise()

#####################################
#####################################

# Oregon hex
## Submarine Cable with 500-meter setback
oregon_hex_submarine_cable500 <- oregon_hex[oregon_submarine_cable500, ] %>%
  sf::st_join(x = .,
              y = oregon_submarine_cable500,
              join = st_intersects)

## Submarine Cable with 501-1000-meter setback
oregon_hex_submarine_cable1000 <- oregon_hex[oregon_submarine_cable1000, ] %>%
  sf::st_join(x = .,
              y = oregon_submarine_cable1000,
              join = st_intersects)


#####################################
#####################################

# Export data
## Submodel geopackage
sf::st_write(obj = oregon_hex_submarine_cable500, dsn = industry_operations_submodel, layer = paste0(region, "_hex_", layer_500), append = F)
sf::st_write(obj = oregon_hex_submarine_cable1000, dsn = industry_operations_submodel, layer = paste0(region, "_hex_", layer_1000), append = F)

## Submarine Cable geopackage
sf::st_write(obj = oregon_hex_submarine_cable500, dsn = submarine_cable_gpkg, layer = paste0(region, "_hex_", layer_500), append = F)
sf::st_write(obj = oregon_hex_submarine_cable1000, dsn = submarine_cable_gpkg, layer = paste0(region, "_hex_", layer_1000), append = F)

sf::st_write(obj = submarine_cables_noaa, dsn = submarine_cable_gpkg, layer = "noaa_submarine_cable", append = F)
sf::st_write(obj = oregon_submarine_cables, dsn = submarine_cable_gpkg, layer = "oregon_submarine_cables", append = F)
sf::st_write(obj = oregon_submarine_cable500, dsn = submarine_cable_gpkg, layer = "oregon_submarine_cable500", append = F)
sf::st_write(obj = oregon_submarine_cable1000, dsn = submarine_cable_gpkg, layer = "oregon_submarine_cable1000", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
