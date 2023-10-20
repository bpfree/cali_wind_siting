################################
### 18. Constraints Submodel ###
################################

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
constraints_submodel <- "data/c_submodel_data/oregon_constraints.gpkg"

## Output directories
oregon_constraints_gpkg <- "data/c_submodel_data/oregon_constraints.gpkg"
suitability_models <- "data/d_suitability_data/oregon_suitability_model.gpkg"

#####################################

# Inspect available layers and names within constraints submodel geopackage
sf::st_layers(dsn = constraints_submodel,
              do_count = T)

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

## submodel
submodel <- "constraints"

## designate date
date <- format(Sys.time(), "%Y%m%d")

#####################################
#####################################

# Load data
## Oregon hex areas (original data)
oregon_hex <- sf::st_read(dsn = study_area_gpkg,
                          layer = paste(sf::st_layers(dsn = study_area_gpkg,
                                                      do_count = TRUE)[[1]][4]))

## Constraints
### DoD OPNAV
oregon_hex_dod_opnav <- sf::st_read(dsn = constraints_submodel, layer = "oregon_hex_dod_opnav") %>%
  dplyr::mutate(dod_value = 0) %>%
  as.data.frame()

### PACPARS
oregon_hex_pacpars <- sf::st_read(dsn = constraints_submodel, layer = "oregon_hex_pacpars") %>%
  dplyr::mutate(pacpars_value = 0) %>%
  as.data.frame()

#####################################
#####################################

# Create Oregon constraints submodel
oregon_hex_constraints <- oregon_hex %>%
  dplyr::left_join(x = .,
                   y = oregon_hex_dod_opnav,
                   by = "index") %>%
  dplyr::left_join(x = .,
                   y = oregon_hex_pacpars,
                   by = "index") %>%
  dplyr::select(index,
                dod_value,
                pacpars_value)

### ***Warning: there are duplicates of the index
duplicates_verify <- oregon_hex_constraints %>%
  # create frequency field based on index
  dplyr::add_count(index) %>%
  # see which ones are duplicates and verify that values for DoD and PACPARS are equal
  dplyr::filter(n>1) %>%
  # show distinct options (otherwise will get 128 results [duplicated cells x frequency of duplication])
  dplyr::distinct()

# Keep only one result per cell
oregon_constraints <- oregon_hex_constraints %>%
  # group by key fields to reduce duplicates
  dplyr::group_by(index, dod_value, pacpars_value) %>%
  # return only distinct rows (remove duplicates)
  dplyr::distinct() %>%
  # create a field called "constraints" that populates with 0 if either DoD or PACPARS values are 0
  dplyr::mutate(constraints = case_when(dod_value == 0 ~ 0, # when DoD value is 0, constraints gets a 0
                                        # if PACPARS values is 0, constraints get a 0
                                        pacpars_value == 0 ~ 0))

#####################################
#####################################

# Export data
## Suitability
sf::st_write(obj = oregon_constraints, dsn = suitability_models, layer = paste0(region, "_", submodel), append = F)

## Constraints
sf::st_write(obj = oregon_hex_dod_opnav, dsn = oregon_constraints_gpkg, layer = "oregon_hex_dod_opnav", append = F)
sf::st_write(obj = oregon_hex_pacpars, dsn = oregon_constraints_gpkg, layer = "oregon_hex_pacpars", append = F)
sf::st_write(obj = oregon_hex_constraints, dsn = suitability_models, layer = "oregon_hex_constraints", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
