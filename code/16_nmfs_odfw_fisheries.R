######################################
### 16. NOAA NMFS / ODFW Fisheries ###
######################################

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
fisheries_dir <- "data/a_raw_data/nmfs_odfw_fisheries/NMFS_Fisheries_Final"

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Submodel directory
fisheries_submodel <- "data/c_submodel_data/oregon_fisheries_submodel.gpkg"

#### Intermediate directories
##### Fisheries geopackage
fisheries_gpkg <- "data/b_intermediate_data/oregon_fisheries.gpkg"

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

## layer names
layer <- "fisheries"

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

nmfs_odfw_fisheries <- sf::read_sf(dsn = fisheries_dir,
                                   layer = "Static_data_for_each_scenario")

#####################################
#####################################

oregon_fisheries <- nmfs_odfw_fisheries %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910") %>%
  # limit to only Oregon call areas
  rmapshaper::ms_clip(target = .,
                      clip = oregon_call_areas) %>%
  # select fields of importance
  dplyr::select(CellID,
                SCEN_4) %>%
  # rename fields
  dplyr::rename(cell_id = CellID,
                scenario4 = SCEN_4) %>%
  # create field called "layer" and fill with "nmfs / odfw fisheries" for summary
  dplyr::mutate(layer = "nmfs / odfw fisheries")

#####################################
#####################################

# NMFS / ODFW fisheries hex grid
oregon_hex_fisheries <- oregon_hex[oregon_fisheries, ] %>%
  # spatially join continental shelf values to Oregon hex cells
  sf::st_join(x = .,
              y = oregon_fisheries,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(index, layer,
                scenario4) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(index) %>%
  # summarise the fisheries score values
  ## take the minimum value of the fisheries score for any that overlap
  ## ***Note: this will provide the most conservation given that low
  ##          values are less desirable
  dplyr::summarise(fisheries_index = min(scenario4))

#####################################
#####################################

# Export data
## Fisheries submodel
sf::st_write(obj = oregon_hex_fisheries, dsn = fisheries_submodel, layer = paste0(region, "_hex_", layer), append = F)

## Fisheries geopackage
sf::st_write(obj = nmfs_odfw_fisheries, dsn = fisheries_gpkg, layer = "nmfs_odfw_fisheries", append = F)
sf::st_write(obj = oregon_fisheries, dsn = fisheries_gpkg, layer = "oregon_fisheries", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
