##############################
### 22. Wind Submodel ########
##############################

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
## Input directories
study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_submodel <- "data/c_submodel_data/oregon_wind_submodel.gpkg"

## Output directories
### Oregon suitability geopackage
oregon_suitability_gpkg <- "data/d_suitability_data/oregon_suitability_model.gpkg"

### Wind directory
suitability_dir <- "data/d_suitability_data"
dir.create(paste0(suitability_dir, "/",
                  "wind_suitability"))

oregon_wind_dir <- "data/d_suitability_data/wind_suitability"
oregon_wind_suitability <- "data/d_suitability_data/wind_suitability/oregon_wind_suitability.gpkg"

#####################################

# Inspect available layers and names within wind submodel geopackage
sf::st_layers(dsn = wind_submodel,
              do_count = T)

#####################################
#####################################

clean_function <- function(data){
  
  data <- data %>%
    as.data.frame() %>%
    dplyr::select(-geom)
  
  return(data)
}

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

## submodel
submodel <- "wind"

## designate date
date <- format(Sys.time(), "%Y%m%d")

## geometric mean weight
wind_wt <- 1/1

#####################################
#####################################

# Load data
## Oregon hex areas (original data)
oregon_hex <- sf::st_read(dsn = study_area_gpkg,
                          layer = paste(sf::st_layers(dsn = study_area_gpkg,
                                                      do_count = TRUE)[[1]][4]))

## Fisheries
oregon_hex_wind <- sf::st_read(dsn = wind_submodel,
                               layer = "oregon_hex_lcoe_2027")

#####################################
#####################################

# Prepare the fisheries data
oregon_wind_values <- oregon_hex_wind %>%
  clean_function()

#####################################
#####################################

# Calculate geometric mean for industry and operations submodel
oregon_wind <- oregon_hex %>%
  # join the wind values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = oregon_wind_values,
                   by = "index") %>%
  
  # add value of 1 for datasets when hex cell has value of NA
  ## for hex cells not impacted by a particular dataset, that cell gets a value of 1
  ### this indicates  suitability with wind energy development
  dplyr::mutate(across(.cols = lcoe_norm_index,
                       .fns = ~replace(x = .,
                                       list = is.na(.),
                                       # replacement values
                                       values = 1))) %>%
  
  # calculate the geometric mean
  ## geometric mean = nth root of the product of the variable values
  dplyr::mutate(wind_geom_mean = lcoe_norm_index ^ wind_wt) %>%
  
  # select the fields of interest
  dplyr::select(index,
                lcoe_norm_index,
                wind_geom_mean) %>%
  dplyr::rename(wind_value = lcoe_norm_index)

### Check to see if there are any duplicates of the indices
### There are none
wind_duplicates <- oregon_wind %>%
  # create frequency field based on index
  dplyr::add_count(index) %>%
  # see which ones are duplicates
  dplyr::filter(n>1) %>%
  # show distinct options
  dplyr::distinct()

#####################################
#####################################

# Export data
## Suitability
sf::st_write(obj = oregon_wind, dsn = oregon_suitability_gpkg, layer = paste0(region, "_", submodel, "_suitability"), append = F)

## Submodel
### Wind
base::saveRDS(obj = oregon_wind_values, file = paste(oregon_wind_dir, "oregon_wind_values.rds", sep = "/"))

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
