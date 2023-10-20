#############################################
### 11. Deep-sea Coral and Sponge Habitat ###
#############################################

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
coral_sponge_dir <- "data/a_raw_data/deep_sea_coral_sponge_habitat/0276883/1.1/data/0-data/NCCOS_USWestCoast_DSC_Models_2016_2020"
intermediate_dir <- "data/b_intermediate_data"

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Submodel directory
natural_resources_submodel <- "data/c_submodel_data/oregon_natural_resources_submodel.gpkg"

#### Intermediate directories
##### Deep sea coral and sponge habitat geopackage (vector data)
coral_sponge_habitat_gpkg <- "data/b_intermediate_data/oregon_coral_sponge_habitat.gpkg"

##### Deep sea coral and sponge habitat subdirectory (raster data)
dir.create(paste0(intermediate_dir, "/",
                  "deep_sea_coral_sponge_habitat"))

deep_coral_sponge_dir <- "data/b_intermediate_data/deep_sea_coral_sponge_habitat"

#####################################
#####################################

## Create z-shape membership function
### Adapted from https://www.mathworks.com/help/fuzzy/zmf.html
zmf_function <- function(richness){
  
  # calculate minimum value
  min <- min(richness$richness_index)
  
  # calculate maximum value
  max <- max(richness$richness_index)
  
  # create a field and populate with the value determined by the z-shape membership scalar
  richness <- richness %>%
    # calculate the z-shape membership value (more desired values get a score of 1 and less desired values will decrease till 0.01)
    ## ***Note: in other words, habitats with higher richness values will be closer to 0
    dplyr::mutate(z_value = ifelse(richness_index == min, 1, # if value is equal to minimum, score as 1
                                   # if value is larger than minimum but lower than mid-value, calculate based on scalar equation
                                   ifelse(richness_index > min & richness_index < (min + max) / 2, 1 - 2 * ((richness_index - min) / (max - min)) ** 2,
                                          # if value is lower than maximum but larger than than mid-value, calculate based on scalar equation
                                          ifelse(richness_index >= (min + max) / 2 & richness_index < max, 2 * ((richness_index - max) / (max - min)) ** 2,
                                                 # if value is equal to maximum, value is equal to 0.01 [all other values should get an NA]
                                                 ifelse(richness_index == max, 0.01, NA)))))
  
  # return the layer
  return(richness)
}

### ***Note: see NOAA's comments for why the maximum species richness value receives a z-value of 0.01 instead of 0.0

#####################################
#####################################

# Unzip the data for the coral richness (source: https://www.ncei.noaa.gov/archive/archive-management-system/OAS/bin/prd/jquery/download/276883.1.1.tar.gz)
### NCEI: https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0276883
### Metadata: https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0276883;view=iso
#### Link for habitat suitability data only: https://www.nodc.noaa.gov/archive/arc0211/0276883/1.1/data/0-data/NCCOS_USWestCoast_DSC_Models_2016_2020/Coral_Richness.zip

## Examine if the filename contains the pattern "Coral_Richness"
### grepl returns a logic statement when pattern "Coral_Richness" is met in the file
file <- basename(list.files(coral_sponge_dir)[3])

unzip_function <- function(data_dir, file){
  
  if (grepl("Coral_Richness", file)){
  
  # grab text before ".zip" and keep only text before that
  new_dir_name <- sub(".zip", "", file)
  
  # create new directory for data
  new_dir <- file.path(data_dir, new_dir_name)
  
  # unzip the file
  unzip(zipfile = file.path(paste0(coral_sponge_dir, "/", file)),
        # export file to the new data directory
        exdir = new_dir)
  
  # remove original zipped file
  file.remove(file.path(data_dir, file))
  }
}

unzip_function(coral_sponge_dir, file)

#####################################
#####################################

# Set parameters
## setback (buffer) distance
buffer <- 500

## designate region name
region <- "california"

## layer names
layer_high <- "high_habitat_coral_sponge"
layer_robust <- "robust_habitat_coral_sponge"

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

# Extend Oregon call areas
oregon_call_areas_500 <- oregon_call_areas %>%
  # extend call areas by 500 meters
  ## this will help extract the deep-sea coral and sponge richness areas
  sf::st_buffer(dist = buffer)

#####################################

# Create deep-sea coral and sponge species habitat richness data
## set directory
richness_dir <- "data/a_raw_data/deep_sea_coral_sponge_habitat/0276883/1.1/data/0-data/NCCOS_USWestCoast_DSC_Models_2016_2020/Coral_Richness/Coral_Richness"

## Open high and robust high datasets
high_habitat <- terra::rast(paste(richness_dir, "Num_Taxa_Corals_HardSubstrate_with_High_Habitat_Suitability.tif", sep = "/"))
robust_habitat <- terra::rast(paste(richness_dir, "Num_Taxa_Corals_HardSubstrate_with_Robust_High_Habitat_Suitability.tif", sep = "/"))

## Convert raster to points
high_habitat_point <- terra::as.points(x = high_habitat) %>%
  # change to simple feature (sf)
  sf::st_as_sf() %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910") %>%
  # simplify column name to "richness" (this is the first column of the object, thus the colnames(.)[1] means take the first column name from the high_habitat object)
  dplyr::rename(richness = colnames(.)[1]) %>%
  # obtain areas with species (richness >= 1)
  dplyr::filter(richness >= 1) %>%
  # add field "layer" and populate with "deep-sea coral and sponge (high habitat)"
  dplyr::mutate(layer = "deep-sea coral and sponge (high habitat)") %>%
  # limit to the extended call areas (with 500m setback)
  rmapshaper::ms_clip(clip = oregon_call_areas_500)

## Create point setbacks and limit to original call areas
oregon_high_habitat_point <- high_habitat_point %>%
  # create setback
  sf::st_buffer(dist = buffer) %>%
  # obtain data within Oregon call areas
  rmapshaper::ms_clip(clip = oregon_call_areas)

#####################################

## Robust high habitat suitability
robust_habitat_point <- terra::as.points(x = robust_habitat) %>%
  sf::st_as_sf() %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910") %>%
  # simplify column name to "richness"
  dplyr::rename(richness = colnames(.)[1]) %>%
  # obtain areas with species (richness >= 1)
  dplyr::filter(richness >= 1) %>%
  # add field "layer" and populate with "deep-sea coral and sponge (robust high habitat)"
  dplyr::mutate(layer = "deep-sea coral and sponge (robust high habitat)") %>%
  # limit to the extended call areas (with 500m setback)
  rmapshaper::ms_clip(clip = oregon_call_areas_500)

## Create point setbacks and limit to original call areas
oregon_robust_habitat_point <- robust_habitat_point %>%
  # create setback
  sf::st_buffer(dist = buffer) %>%
  # obtain data within Oregon call areas
  rmapshaper::ms_clip(clip = oregon_call_areas)

#####################################
#####################################

# Deep-sea coral and sponge habitat in call area
## High habitat suitability
oregon_hex_coral_sponge_high <- oregon_hex[oregon_high_habitat_point, ] %>%
  # spatially join NMFS EFHCA values to Oregon hex cells
  sf::st_join(x = .,
              y = oregon_high_habitat_point,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(index, layer,
                richness) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(index) %>%
  # summarise the richness values
  ## take the minimum value of the richness for any that overlap
  ## ***Note: this will provide the most conservation given that low
  ##          values are less desirable
  dplyr::summarise(richness_index = min(richness)) %>%
  # z-membership function to add the z-value based on the richness values
  zmf_function(.)
  
#####################################

## Robust high habitat suitability
oregon_hex_coral_sponge_robust <- oregon_hex[oregon_robust_habitat_point, ] %>%
  # spatially join NMFS EFHCA values to Oregon hex cells
  sf::st_join(x = .,
              y = oregon_robust_habitat_point,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(index, layer,
                richness) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(index) %>%
  # summarise the richness values
  ## take the minimum value of the richness for any that overlap
  ## ***Note: this will provide the most conservation given that low
  ##          values are less desirable
  dplyr::summarise(richness_index = min(richness)) %>%
  # z-membership function to add the z-value based on the richness values
  zmf_function(.)

#####################################
#####################################

# Export data
## Natural resources submodel
sf::st_write(obj = oregon_hex_coral_sponge_high, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", layer_high, "_", buffer, "m"), append = F)
sf::st_write(obj = oregon_hex_coral_sponge_robust, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", layer_robust, "_", buffer, "m"), append = F)

## Deep-sea coral and sponge habitat
### Richness point
sf::st_write(obj = high_habitat_point, dsn = coral_sponge_habitat_gpkg, layer = "coral_sponge_high_habitat_point", append = F)
sf::st_write(obj = robust_habitat_point, dsn = coral_sponge_habitat_gpkg, layer = "coral_sponge_robust_habitat_point", append = F)

### Oregon
sf::st_write(obj = oregon_high_habitat_point, dsn = coral_sponge_habitat_gpkg, layer = "coral_sponge_high_habitat_point", append = F)
sf::st_write(obj = oregon_robust_habitat_point, dsn = coral_sponge_habitat_gpkg, layer = "coral_sponge_robust_habitat_point", append = F)

## Oregon call area hex grid
sf::st_write(obj = oregon_hex_coral_sponge_high, dsn = coral_sponge_habitat_gpkg, layer = "oregon_hex_high_habitat_coral_sponge", append = F)
sf::st_write(obj = oregon_hex_coral_sponge_robust, dsn = coral_sponge_habitat_gpkg, layer = "oregon_hex_robust_high_habitat_coral_sponge", append = F)

terra::writeRaster(high_habitat, filename = file.path(deep_coral_sponge_dir, "coral_sponge_high_habitat.grd"), overwrite = T)
terra::writeRaster(robust_habitat, filename = file.path(deep_coral_sponge_dir, "coral_sponge_robust_high_habitat.grd"), overwrite = T)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
