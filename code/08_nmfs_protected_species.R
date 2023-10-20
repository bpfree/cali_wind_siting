#################################
### 8. NMFS Protected Species ###
#################################

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
prd_species <- "data/a_raw_data/PRD_combined_layer.gdb"

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Analysis directories
natural_resources_submodel <- "data/c_submodel_data/oregon_natural_resources_submodel.gpkg"

#####################################
#####################################

# Inspect available layers and names
sf::st_layers(dsn = prd_species,
              do_count = TRUE)

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

## layer names
export_leatherback <- "leatherback"
export_humpback_ca <- "humpback_ca_dps"
export_humpack_mx <- "humpback_mexico_dps"
export_killer_whale <- "killer_whale"
export_blue_whale <- "blue_whale"

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

## Species data
leatherback <- sf::st_read(dsn = prd_species, layer = "PRD_Leatherback") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")
killer_whale <- sf::st_read(dsn = prd_species, layer = "PRD_KillerWhale") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")
humpback_ca <- sf::st_read(dsn = prd_species, layer = "PRD_Humpback_CEN") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")
humpback_mx <- sf::st_read(dsn = prd_species, layer = "PRD_Humpback_Mex") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")
blue_whale <- sf::st_read(dsn = prd_species, layer = "PRD_BlueWhale") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

#####################################
#####################################

# Oregon hex grid by species
oregon_hex_leatherback <- oregon_hex[leatherback, ] %>%
  # spatially join protected species values to Oregon hex cells 
  sf::st_join(x = .,
              y = leatherback,
              join = st_intersects)

oregon_hex_humpback_ca_dps <- oregon_hex[humpback_ca, ] %>%
  # spatially join protected species values to Oregon hex cells 
  sf::st_join(x = .,
              y = humpback_ca,
              join = st_intersects) %>%
  # due to overlapping exclusion areas there are a few duplicated indexes
  ## group by unique indexes
  dplyr::group_by(index) %>%
  ## summarise to remove duplicates
  dplyr::summarise()

oregon_hex_humpback_mexico_dps <- oregon_hex[humpback_mx, ] %>%
  # spatially join protected species values to Oregon hex cells 
  sf::st_join(x = .,
              y = humpback_mx,
              join = st_intersects) %>%
  # due to overlapping exclusion areas there are a few duplicated indexes
  ## group by unique indexes
  dplyr::group_by(index) %>%
  ## summarise to remove duplicates
  dplyr::summarise()

oregon_hex_killer_whale <- oregon_hex[killer_whale, ] %>%
  # spatially join protected species values to Oregon hex cells 
  sf::st_join(x = .,
              y = killer_whale,
              join = st_intersects) %>%
  # due to overlapping exclusion areas there are a few duplicated indexes
  ## group by unique indexes
  dplyr::group_by(index) %>%
  ## summarise to remove duplicates
  dplyr::summarise()

oregon_hex_blue_whale <- oregon_hex[blue_whale, ] %>%
  # spatially join protected species values to Oregon hex cells 
  sf::st_join(x = .,
              y = blue_whale,
              join = st_intersects)

#####################################

# Non-protected species areas
oregon_protected_species_areas <- leatherback %>%
  rbind(killer_whale,
        humpback_ca,
        humpback_mx,
        blue_whale)

oregon_hex_protected_species <- oregon_hex[oregon_protected_species_areas, ] %>%
  # spatially join protected species values to Oregon hex cells 
  sf::st_join(x = .,
              y = oregon_protected_species_areas,
              join = st_intersects) %>%
  # due to overlapping exclusion areas there are a few duplicated indexes
  ## group by unique indexes
  dplyr::group_by(index) %>%
  ## summarise to remove duplicates
  dplyr::summarise()

oregon_hex_non_protected <- oregon_hex %>%
  rmapshaper::ms_erase(oregon_hex_protected_species)

#####################################
#####################################

# Export data
## Natural resources submodel
### Species
sf::st_write(obj = oregon_hex_leatherback, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", export_leatherback), append = F)
sf::st_write(obj = oregon_hex_humpback_ca_dps, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", export_humpback_ca), append = F)
sf::st_write(obj = oregon_hex_humpback_mexico_dps, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", export_humpack_mx), append = F)
sf::st_write(obj = oregon_hex_killer_whale, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", export_killer_whale), append = F)
sf::st_write(obj = oregon_hex_blue_whale, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", export_blue_whale), append = F)

### Non-protected species areas
sf::st_write(obj = oregon_hex_non_protected, dsn = natural_resources_submodel, layer = "oregon_hex_non_protected_species", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
