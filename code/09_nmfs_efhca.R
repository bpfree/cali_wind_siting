###################################
### 9. NMFS Habitat Data Layers ###
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
nmfs_efhca_dir <- "data/a_raw_data/westcoast_efha"
study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Submodel directory
natural_resources_submodel <- "data/c_submodel_data/oregon_natural_resources_submodel.gpkg"

#### Intermediate directories
study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
efhca_gpkg <- "data/b_intermediate_data/oregon_efhca.gpkg"

#####################################
#####################################

# Set parameters
## setback (buffer) distance
buffer <- 500

## designate region name
region <- "california"

## layer names
layer <- "efhca"

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

## NMFS Essential Fish Habitat Conservation Areas (source: https://www.habitat.noaa.gov/protection/efh/newInv/data/west_coast/westcoast_efha.zip)
### Text: https://www.ecfr.gov/current/title-50/chapter-VI/part-660/subpart-C/section-660.76
### Alternative download source: https://www.fisheries.noaa.gov/s3/2021-02/EFH-HAPC-EFHCA-shapefiles-AM19-2006-AM28-2020.zip
### Essential Fish Habitat Mapper: https://www.habitat.noaa.gov/apps/efhmapper/?data_id=dataSource_13-EFHA_7887%3A102&page=page_4
#### ***Note: In particular, see the "EFH Areas Protected From Fishing"
nmfs_efhca_data <- sf::st_read(dsn = nmfs_efhca_dir, layer = "westcoast_efha") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

#####################################
#####################################

# NMFS Essential Fish Habitat Conservation Area in Oregon call areas
oregon_nmfs_efhca <- nmfs_efhca_data %>%
  # filter out 700 fathom areas and deep-sea conservation ecosystem conservation areas
  dplyr::filter(!OBJECTID == c(1016, 1017)) %>%
  # obtain NMFS EFCA within Oregon call areas
  rmapshaper::ms_clip(target = .,
                      clip = oregon_call_areas) %>%
  # add a 500-meter buffer (setback distance)
  sf::st_buffer(dist = buffer)

#####################################
#####################################

# NMFS EFHCA hex grids
oregon_hex_efhca <- oregon_hex[oregon_nmfs_efhca, ] %>%
  # spatially join NMFS EFHCA values to Oregon hex cells
  sf::st_join(x = .,
              y = oregon_nmfs_efhca,
              join = st_intersects) %>%
  # add field "layer" and populate with "efhca"
  dplyr::mutate(layer = "efhca") %>%
  # select fields of importance
  dplyr::select(index, layer)

#####################################
#####################################

# Export data
## Analysis geopackage
sf::st_write(oregon_hex_efhca, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", layer, "_", buffer, "m"), append = F)

## EFHCA geopackage
sf::st_write(nmfs_efhca_data, dsn = efhca_gpkg, layer = "nmfs_efhca", append = F)

sf::st_write(oregon_nmfs_efhca, dsn = efhca_gpkg, layer = "oregon_nmfs_efhca_500m", append = F)
sf::st_write(oregon_hex_efhca, dsn = efhca_gpkg, layer = "oregon_hex_efhca_500m", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
