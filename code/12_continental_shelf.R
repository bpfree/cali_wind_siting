###################################
### 12. Continental Shelf Break ###
###################################

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

#####################################
#####################################

# Set directories
## Define data directory (as this is an R Project, pathnames are simplified)
### Input directories
shelf_gpkg <- "data/a_raw_data/NCCOS_Share.gdb"

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

## Output directories
### Analysis directory
natural_resources_submodel <- "data/c_submodel_data/oregon_natural_resources_submodel.gpkg"

### Intermediate directory
continental_shelf_geopackage <- "data/b_intermediate_data/oregon_continental_shelf.gpkg"

#####################################

# Inspect available layers and names
sf::st_layers(dsn = shelf_gpkg,
              do_count = T)

#####################################
#####################################

# Set parameters
## setback (buffer) distance
buffer <- 10000

## designate region name
region <- "california"

## layer names
layer <- "continental_shelf"

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

## ***Note: Curt Whitmire (curt.whitmire@noaa.gov) generated the continental shelf break data by selecting the
##          shelf polygon from the surficial geological habitat version 4 data, converting to polyline features
##          and then cut by the northern and southern extents.
### These data were manually created

continental_shelf <- sf::st_read(dsn = shelf_gpkg, layer = "V4_0_SGH_WA_OR_NCA_Shelf_line_diss") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910") # EPSG 26910 (https://epsg.io/26910)

#####################################
#####################################

# 10-km buffer around continental shelf
continental_shelf_10km <- continental_shelf %>%
  # create field called "layer" and fill with "continental shelf" for summary
  dplyr::mutate(layer = "continental shelf") %>%
  # group by layer and summarise to have single feature
  dplyr::group_by(layer) %>%
  dplyr::summarise() %>%
  # apply a buffer of 10km (10000m)
  sf::st_buffer(dist = buffer)

plot(continental_shelf_10km)

# Buffered continental shelf within Oregon call areas
oregon_continental_shelf <- continental_shelf_10km %>%
  # limit to Oregon call areas
  rmapshaper::ms_clip(target = .,
                      clip = oregon_call_areas)

plot(oregon_continental_shelf)

#####################################
#####################################

# Continental shelf hex grids
oregon_hex_continental_shelf <- oregon_hex[oregon_continental_shelf, ]  %>%
  # spatially join continental shelf values to Oregon hex cells
  sf::st_join(x = .,
              y = oregon_continental_shelf,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(index, layer)

#####################################
#####################################

# Export data
## Natural resources submodel
sf::st_write(obj = oregon_hex_continental_shelf, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", layer, "_10km"), append = F)

## Continental shelf geopackage
sf::st_write(obj = continental_shelf, dsn = natural_resources_submodel, layer = "continental_shelf", append = F)
sf::st_write(obj = continental_shelf_10km, dsn = natural_resources_submodel, layer = "continental_shelf_10km", append = F)
sf::st_write(obj = oregon_continental_shelf, dsn = natural_resources_submodel, layer = "oregon_continental_shelf_10km", append = F)
sf::st_write(obj = oregon_hex_continental_shelf, dsn = natural_resources_submodel, layer = "oregon_hex_continental_shelf_10km", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
