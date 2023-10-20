#######################################
### 3. DoD Offshore Wind Assessment ###
#######################################

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
dod_opnav_dir <- "data/a_raw_data/opnav_combined_assessment_May2022"

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Constraint submodel
constraint_gpkg <- "data/c_submodel_data/constraints.gpkg"

#### Intermediate directories
intermediate_dir <- "data/b_intermediate_data"
dod_gpkg <- "data/b_intermediate_data/oregon_dod_opnav.gpkg"

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

## layer name
layer <- "dod_opnav"

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
#####################################

## DoD combined Oregon offshore wind assessment (OPNAV May 2022)
### ***Note: currently these data are not available for download 
### A data portal contains the data for viewing: https://offshorewind.westcoastoceans.org/
### To view the data: Human > Military > Combined Oregon Offshore Wind Assessment, OPNAV, May 2022
### Metadata: https://www.coastalatlas.net/waf/boem/OPNAV_CombinedAssesment_May2022.xml
### Dataset Source Provider: https://gis.lcd.state.or.us/server/rest/services/Projects/OCMP_OceanPlanning_Human/MapServer/21
### Alternative link on dataset: https://portal.westcoastoceans.org/geoportal/rest/metadata/item/45b6aa29abe7427a91d8f430eac0ab75/html
### InPort: https://www.fisheries.noaa.gov/inport/item/48875

dod_opnav <- sf::st_read(dsn = dod_opnav_dir, layer = "OPNAV_CombinedAssesment_May2022") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

#####################################
#####################################

# Exclusion areas
## Explanation of data
### Green = no mitigation necessary -- open areas
### Red = mitigation not feasible -- closed areas

dod_opnav_exclusion <- dod_opnav %>%
  # limit areas to only ones that will be constraints (no mitigation feasible)
  dplyr::filter(Color == "Red") %>%
  # obtain only DoD exclusion areas in the study area
  rmapshaper::ms_clip(target = .,
                      clip = oregon_call_areas) %>%
  # create field called "layer" and fill with "DoD exclusion" for summary
  dplyr::mutate(layer = "DoD exclusion")

#####################################

oregon_hex_dod_opnav <- oregon_hex[dod_opnav_exclusion, ] %>%
  # spatially join DoD OPNAV data to Oregon hex cells
  sf::st_join(x = .,
              y = dod_opnav_exclusion,
              join = st_intersects) %>%
  # select fields of interest
  dplyr::select(index, layer)

#####################################
#####################################

# Export data
## Constraints
sf::st_write(obj = oregon_hex_dod_opnav, dsn = constraint_gpkg, layer = paste0(region, "_hex_", layer), append = F)

## DoD OPNAV geopackage
sf::st_write(obj = dod_opnav, dsn = dod_gpkg, layer = "oregon_call_area_dod_opnav", append = F)
sf::st_write(obj = dod_opnav_exclusion, dsn = dod_gpkg, layer = "oregon_call_area_dod_opnav_exclusion", append = F)
sf::st_write(obj = oregon_hex_dod_opnav, dsn = dod_gpkg, layer = "oregon_hex_dod_opnav", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
