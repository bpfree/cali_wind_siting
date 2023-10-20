############################
### 2. Define Study Area ###
############################

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

# Commentary on R and code formulation:
## ***Note: If not familiar with dplyr notation
##          dplyr is within the tidyverse and can use %>%
##          to "pipe" a process, allowing for fluidity
## Can learn more here: https://style.tidyverse.org/pipes.html

## Another  common coding notation used is "::"
## For instance, you may encounter it as dplyr::filter()
## This means use the filter function from the dplyr package
## Notation is used given sometimes different packages have
## the same function name, so it helps code to tell which
## package to use for that particular function.
## The notation is continued even when a function name is
## unique to a particular package so it is obvious which
## package is used

#####################################
#####################################

# Set directories
## Define data directory (as this is an R Project, pathnames are simplified)
### Input directories
wind_area_dir <- "data/a_raw_data/BOEM-Renewable-Energy-Geodatabase/BOEMWindLayers_4Download.gdb"
original_hex_gdb <- "data/a_raw_data/Draft_5_final.gdb"

### Output directories
#### Intermediate directories
study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

#####################################
#####################################

# Inspect available layers and names within BOEM geodatabase
## Note: as of October 6th 2023, there are 5 elements within the BOEM geodatabase;
##       previously there were 4 elements and the order was different
sf::st_layers(dsn = wind_area_dir,
              do_count = TRUE)

# Inspect available layers and names within original final Oregon study geodatabase
sf::st_layers(dsn = original_hex_gdb,
              do_count = TRUE)

#####################################
#####################################

# Load original data
original_hex_grid <- sf::st_read(dsn = original_hex_gdb, layer = "most_conservative_model") %>%
  # add field "index" that will be populated with the row_number
  dplyr::mutate(index = row_number()) %>%
  dplyr::select(index)

#####################################
#####################################

# Load call areas
wind_areas <- sf::st_read(dsn = wind_area_dir,
                          # select the layer for call areas and use to open desired layer (planning area outlines)
                          ## ***Note: the planning area outlines is the 2nd dataset
                          ## Using this method will avoid having to update the layer name each time BOEM updates the dataset
                          ## ***Note: If BOEM changes which layers are available it is possible in the future that the planning
                          ## area outlines will no longer be the 2nd dataset
                          layer = paste(sf::st_layers(dsn = wind_area_dir,
                                                # [[1]] --> first component, which is the column "layer_name"
                                                # [2] --> 2nd element of that list, which is the planning area outlines
                                                do_count = TRUE)[[1]][2])) %>%
  # filter for only Oregon call areas
  dplyr::filter(grepl(pattern = "california",
                      x = ADDITIONAL_INFORMATION)) %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

## verify that units are in meters and no longer degrees
st_crs(wind_areas, parameters = TRUE)$units_gdal

#####################################
#####################################

# Hexagon area = ((3 * sqrt(3))/2) * side length ^ 2 (https://pro.arcgis.com/en/pro-app/latest/tool-reference/data-management/generatetesellation.htm)
# 1 acre equals approximately 4064.86 square meters
# 10 acres = 40468.6 square meters
# 40468.6 = ((3 * sqrt(3))/2) * side length ^ 2
# 40468.6 * 2 = 3 * sqrt(3) * side length ^ 2 --> 80937.2
# 80937.2 / 3 = sqrt(3) * side length ^ 2 --> 26979.07
# 26979.07 ^ 2 = 3 * side length ^ 4 --> 727870218
# 727870218 / 3 = side length ^ 4 --> 242623406
# 242623406 ^ (1/4) = side length --> 124.8053

# Create 10-acre grid around call areas
wind_area_grid <- sf::st_make_grid(x = wind_areas,
                                   ## see documentation on what cellsize means when relating to hexagons: https://github.com/r-spatial/sf/issues/1505
                                   ## cellsize is the distance between two vertices (short diagonal --> d = square root of 3 * side length)
                                   ### So in this case, square-root of 3 * 124.8053 = 1.73205080757 * 124.8053 = 216.1691
                                   cellsize = 216.1691,
                                   # make hexagon (TRUE will generate squares)
                                   square = FALSE,
                                   # make hexagons orientation with a flat topped (FALSE = pointy top)
                                   flat_topped = TRUE) %>%
  # convert back as sf
  sf::st_as_sf()

#####################################
#####################################

# Subset by location: hexagonal grids that intersect with wind areas
wind_area_hex <- wind_area_grid[wind_areas, ] %>%
  # add field "index" that will be populated with the row_number
  dplyr::mutate(index = row_number())

# Oregon call area hexes as single feature
## ***Note: This dataset will be used to extract any data from datasets
##          within the model that will impact wind siting suitability
oregon_call_area_hex_dissolve <- wind_area_hex %>%
  # create field called "call area"
  dplyr::mutate(call_area = "call_area") %>%
  # group all rows by the different elements with "call area" field -- this will create a row for the grouped data
  dplyr::group_by(call_area) %>%
  # summarise all those grouped elements together -- in effect this will create a single feature
  dplyr::summarise()

#####################################
#####################################

# Export data
## Study Area
sf::st_write(obj = wind_area_grid, dsn = study_area_gpkg, layer = paste(region, "call_area_grid", sep = "_"), append = F)
sf::st_write(obj = wind_area_hex, dsn = study_area_gpkg, layer = paste(region, "call_area_hex", sep = "_"), append = F)
sf::st_write(obj = oregon_call_area_hex_dissolve, dsn = study_area_gpkg, layer = paste(region, "call_areas_dissolve", sep = "_"), append = F)

sf::st_write(obj = original_hex_grid, dsn = study_area_gpkg, layer = paste(region, "original_hex_grid", sep = "_"), append = F)

## Wind Call Areas
sf::st_write(obj = wind_areas, dsn = wind_area_gpkg, layer = "oregon_wind_call_areas", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
