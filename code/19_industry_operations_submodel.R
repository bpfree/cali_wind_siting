############################################
### 19. Industry and Operations Submodel ###
############################################

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
industry_operations_submodel <- "data/c_submodel_data/oregon_industry_operations_submodel.gpkg"

## Output directories
### Oregon suitability geopackage
oregon_suitability_gpkg <- "data/d_suitability_data/oregon_suitability_model.gpkg"

### Industry and operations directory
suitability_dir <- "data/d_suitability_data"
dir.create(paste0(suitability_dir, "/",
                  "industry_operations_suitability"))

oregon_industry_operations_dir <- "data/d_suitability_data/industry_operations_suitability"
oregon_industry_operations_suitability <- "data/d_suitability_data/industry_operations_suitability/oregon_industry_operations_suitability.gpkg"

#####################################

# Inspect available layers and names within industry and operations submodel geopackage
sf::st_layers(dsn = industry_operations_submodel,
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
submodel <- "industry_operations"

## designate date
date <- format(Sys.time(), "%Y%m%d")

## geometric mean weight
io_wt <- 1/2

#####################################
#####################################

# Load data
## Oregon hex areas (original data)
oregon_hex <- sf::st_read(dsn = study_area_gpkg,
                          layer = paste(sf::st_layers(dsn = study_area_gpkg,
                                                      do_count = TRUE)[[1]][4]))

## Industry and Operations
### Submarine cables
oregon_hex_submarine_cable500 <- sf::st_read(dsn = industry_operations_submodel, layer = "oregon_hex_submarine_cable500m") %>%
  dplyr::mutate(sc500_value = 0.6) %>%
  as.data.frame() %>%
  dplyr::select(-geom)

oregon_hex_submarine_cable1000 <- sf::st_read(dsn = industry_operations_submodel, layer = "oregon_hex_submarine_cable1000m") %>%
  dplyr::mutate(sc1000_value = 0.8) %>%
  as.data.frame() %>%
  dplyr::select(-geom)

### Scientific surveys
oregon_hex_eastwest_survey_corridors <- sf::st_read(dsn = industry_operations_submodel,
                                                    layer = "oregon_hex_eastwest_survey_corridors_4nm") %>%
  dplyr::mutate(eastwest_value = 0.01) %>%
  as.data.frame() %>%
  dplyr::select(-geom)

oregon_hex_additional_eastwest_survey_corridors <- sf::st_read(dsn = industry_operations_submodel,
                                                               layer = "oregon_hex_additional_eastwest_survey_corridors_4nm") %>%
  dplyr::mutate(eastwest_add_value = 0.5) %>%
  as.data.frame() %>%
  dplyr::select(-geom)

oregon_hex_survey_stations <- sf::st_read(dsn = industry_operations_submodel, layer = "oregon_hex_survey_stations_2nm") %>%
  dplyr::mutate(sstat_value = 0.5) %>%
  as.data.frame() %>%
  dplyr::select(-geom)

oregon_hex_survey_transects <- sf::st_read(dsn = industry_operations_submodel, layer = "oregon_hex_survey_transects_1nm") %>%
  dplyr::mutate(stransect_value = 0.5) %>%
  as.data.frame() %>%
  dplyr::select(-geom)

#####################################
#####################################

# Calculate geometric mean for industry and operations submodel
oregon_industry_operations <- oregon_hex %>%
  # join the submarine cable (0 - 500m) values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = oregon_hex_submarine_cable500,
                   by = "index") %>%
  # join the submarine cable (501 - 1000m) values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = oregon_hex_submarine_cable1000,
                   by = "index") %>%
  # join the east-west survey corridors values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = oregon_hex_eastwest_survey_corridors,
                   by = "index") %>%
  # join the additional east-west survey corridors values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = oregon_hex_additional_eastwest_survey_corridors,
                   by = "index") %>%
  # join the survey station values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = oregon_hex_survey_stations,
                   by = "index") %>%
  # join the survey transects values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = oregon_hex_survey_transects,
                   by = "index") %>%
  # select fields of interest
  dplyr::select(index,
                sc500_value,
                sc1000_value,
                eastwest_value,
                eastwest_add_value,
                sstat_value,
                stransect_value) %>%
  
  # add value of 1 for datasets when hex cell has value of NA
  ## for hex cells not impacted by a particular dataset, that cell gets a value of 1
  ### this indicates  suitability with wind energy development
  dplyr::mutate(across(2:7, ~replace(x = .,
                                     list = is.na(.),
                                     # replacement values
                                     values = 1))) %>%
  
  # create combined submarine cable field
  ## fill with values of 500-m buffer, 500-1000-m buffer
  ## when a hex has both values give it the lower value (500-m buffer = 0.6)
  dplyr::mutate(sub_cable_value = sc500_value + sc1000_value) %>%
  dplyr::mutate(sub_cable_value = case_when(sub_cable_value == 2 ~ 1,
                                            sub_cable_value == 1.8 ~ 0.8,
                                            sub_cable_value == 1.6 ~ 0.6,
                                            sub_cable_value == 1.4 ~ 0.6,
                                            sub_cable_value == 0.8 ~ 0.8,
                                            sub_cable_value == 0.6 ~ 0.6)) %>%
                                      
  dplyr::relocate(sub_cable_value,
                  .after = sc1000_value) %>%
  
  # calculate a summary value for scientific surveys
  ## when a hex grid cell has multiple values the minimum
  ## value across the datasets is assigned to the new
  ## summarized field
  dplyr::mutate(sci_survey_value = pmin(eastwest_value,
                                        eastwest_add_value,
                                        sstat_value,
                                        stransect_value,
                                        # remove any values that are NA when new field
                                        na.rm = T)) %>%
  dplyr::relocate(sci_survey_value,
                  .after = stransect_value) %>%

  ## geometric mean = nth root of the product of the variable values
  dplyr::mutate(io_geom_mean = (sub_cable_value ^ io_wt) * (sci_survey_value ^ io_wt)) %>%
  
  # relocate the industry and operations geometric mean field
  dplyr::relocate(io_geom_mean,
                  .after = sci_survey_value)

#list(unique(oregon_industry_operations$sci_survey)

### Check to see if there are any duplicates of the indices
### There are none
io_duplicates <- oregon_industry_operations %>%
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
sf::st_write(obj = oregon_industry_operations, dsn = oregon_suitability_gpkg, layer = paste0(region, "_", submodel, "_suitability"), append = F)

## Submodel
saveRDS(obj = oregon_hex_submarine_cable500, file = paste(oregon_industry_operations_dir, "oregon_hex_submarine_cable_500m.rds", sep = "/"))
saveRDS(obj = oregon_hex_submarine_cable1000, file = paste(oregon_industry_operations_dir, "oregon_hex_submarine_cable_1000m.rds", sep = "/"))

saveRDS(obj = oregon_hex_eastwest_survey_corridors, file = paste(oregon_industry_operations_dir, "oregon_hex_eastwest_survey_corridors_4nm.rds", sep = "/"))
saveRDS(obj = oregon_hex_additional_eastwest_survey_corridors, file = paste(oregon_industry_operations_dir, "oregon_hex_additional_eastwest_survey_corridors_4nm.rds", sep = "/"))
saveRDS(obj = oregon_hex_survey_stations, file = paste(oregon_industry_operations_dir, "oregon_hex_survey_stations_2nm.rds", sep = "/"))
saveRDS(obj = oregon_hex_survey_transects, file = paste(oregon_industry_operations_dir, "oregon_hex_submarine_cable_1nm.rds", sep = "/"))

sf::st_write(obj = oregon_industry_operations, dsn = oregon_industry_operations_suitability, layer = "oregon_industry_operations_suitability", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
