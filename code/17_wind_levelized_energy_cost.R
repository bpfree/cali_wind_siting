####################################
### 17. Levelized Cost of Energy ###
####################################

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
lcoe_dir <- "data/a_raw_data/LCOE_2027"

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Submodel directory
wind_submodel <- "data/c_submodel_data/oregon_wind_submodel.gpkg"

#### Intermediate directories
wind_gpkg <- "data/b_intermediate_data/wind.gpkg"

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

## layer names
layer <- "lcoe_2027"

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

## Levelized cost of energy (2027)
### ***Note: These data are not presently publicly available, yet a presentation of the data
###          are viewable here: https://www.nrel.gov/docs/fy22osti/80908.pdf (Musial et al. 2021).
###          Additional information can be located here: https://www.osti.gov/biblio/1825061
lcoe_2027 <- sf::st_read(dsn = lcoe_dir, "LCOE_2027") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

#####################################
#####################################

# Energy cost in Oregon call areas
oregon_lcoe_2027 <- lcoe_2027 %>%
  # limit to only Oregon call areas
  rmapshaper::ms_clip(target = .,
                      clip = oregon_call_areas) %>%
  # make sure all features are polygons
  sf::st_cast("MULTIPOLYGON")
  

#####################################

# Normalize the levelized cost data
## The analysis seeks to have values get normalized with a linear function
## and place them between 0.8 and 1.0
### Linear function = (value - min) / (max - min)
### Linear function between values = ((value - min) / (max - min)) * (target_max - target_min) + target_min

## Targets
### target minimum
tmin <- 0.8

### target maximum
tmax <- 1.0

## Values
### value minimum
min <- min(oregon_lcoe_2027$lcoe_mid)

### value maximum
max <- max(oregon_lcoe_2027$lcoe_mid)

## Oregon levelized cost of energy normalized
oregon_lcoe_2027_norm <- oregon_lcoe_2027 %>%
  # normalize data between 0.8 and 1
  dplyr::mutate(lcoe_norm = 
                  # flip the normalization so high costs get low
                  # scores (0.8) and low costs get higher scores (1.0)
                  (tmax + tmin) - 
                  # normalize the data (will become between 0.8 and 1.0)
                  (((lcoe_mid - min) / (max - min)) *
                  # and then rescale to be between 0.8 and 1.0
                    (tmax - tmin) + tmin)) %>%
  # create field called "layer" and fill with "levelized cost of energy (2027)" for summary
  dplyr::mutate(layer = "levelized cost of energy (2027)") %>%
  # select fields of interest
  dplyr::select(Id,
                layer,
                lcoe_mid,
                lcoe_norm)

#####################################
#####################################

# Levelized cost of energy (2027) hex grid
oregon_hex_lcoe_2027 <- oregon_hex[oregon_lcoe_2027_norm, ] %>%
  # spatially join continental shelf values to Oregon hex cells
  sf::st_join(x = .,
              y = oregon_lcoe_2027_norm,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(index, layer,
                lcoe_norm) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(index) %>%
  # summarise the normalized values for the levelized cost of energy in 2027
  ## take the minimum value of the index
  dplyr::summarise(lcoe_norm_index = min(lcoe_norm))

#####################################
#####################################

# Export data
## Wind submodel
sf::st_write(obj = oregon_hex_lcoe_2027, dsn = wind_submodel, layer = paste0(region, "_hex_", layer), append = F)

## Wind geopackage
sf::st_write(obj = oregon_lcoe_2027, dsn = wind_gpkg, layer = "oregon_lcoe_2027", append = F)
sf::st_write(obj = oregon_lcoe_2027_norm, dsn = wind_gpkg, layer = "oregon_lcoe_2027_norm", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
