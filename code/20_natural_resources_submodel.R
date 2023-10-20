######################################
### 20. Natural Resources Submodel ###
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
## Input directories
study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
natural_resources_submodel <- "data/c_submodel_data/oregon_natural_resources_submodel.gpkg"

## Output directories
### Oregon suitability geopackage
oregon_suitability_gpkg <- "data/d_suitability_data/oregon_suitability_model.gpkg"

### Natural resources directory
suitability_dir <- "data/d_suitability_data"
dir.create(paste0(suitability_dir, "/",
                  "natural_resources_suitability"))

oregon_natural_resources_dir <- "data/d_suitability_data/natural_resources_suitability"
dir.create(paste0(oregon_natural_resources_dir, "/",
                  "protected_species"))
dir.create(paste0(oregon_natural_resources_dir, "/",
                  "habitat"))
dir.create(paste0(oregon_natural_resources_dir, "/",
                  "marine_bird"))

oregon_nr_species_dir <- "data/d_suitability_data/natural_resources_suitability/protected_species"
oregon_nr_habitat_dir <- "data/d_suitability_data/natural_resources_suitability/habitat"
oregon_nr_marine_bird_dir <- "data/d_suitability_data/natural_resources_suitability/marine_bird"

#####################################

# Inspect available layers and names within natural resources submodel geopackage
sf::st_layers(dsn = natural_resources_submodel,
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
submodel <- "natural_resources"

## designate date
date <- format(Sys.time(), "%Y%m%d")

## geometric mean weight
nr_wt <- 1/3

#####################################
#####################################

# Load data
## Oregon hex areas (original data)
oregon_hex <- sf::st_read(dsn = study_area_gpkg,
                          layer = paste(sf::st_layers(dsn = study_area_gpkg,
                                                      do_count = TRUE)[[1]][4]))

## Natural Resources
### Protected species
leatherback <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_leatherback")
killer_whale <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_killer_whale")
humpback_ca_dps <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_humpback_ca_dps")
humpback_mexico_dps <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_humpback_mexico_dps")
blue_whale <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_blue_whale")

non_protected_species_areas <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_non_protected_species")

## Habitat layers
efhca <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_efhca_500m")
rocky_reef_mapped <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_rocky_reef_mapped_500m")
rocky_reef_probable <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_rocky_reef_probable_500m")
deep_sea_coral_sponge <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_high_habitat_coral_sponge_500m")
continental_shelf <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_continental_shelf_10km")
methane_bubble <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_methane_bubble_streams_1km")

## Marine bird
marine_bird <- sf::st_read(dsn = natural_resources_submodel, layer = "oregon_hex_marine_bird")

#####################################
#####################################

# Calculate score
## Protected species
### Individual species
#### Leatherback sea turtle
leatherback_value <- leatherback %>%
  dplyr::mutate(leatherback_value = 0.1) %>%
  clean_function()

#### Killer whale
killer_whale_value <- killer_whale %>%
  dplyr::mutate(killerwhale_value = 0.1) %>%
  clean_function()

#### Humpback whale (Central America DPS)
humpback_ca_value <- humpback_ca_dps %>%
  dplyr::mutate(humpback_ca_value = 0.3) %>%
  clean_function()

#### Humpback whale (Mexico DPS)
humpback_mx_value <- humpback_mexico_dps %>%
  dplyr::mutate(humpback_mx_value = 0.5) %>%
  clean_function()

#### Blue whale
bluewhale_value <- blue_whale %>%
  dplyr::mutate(bluewhale_value = 0.2) %>%
  clean_function()

#### Non-protected species areas
non_protected_species_value <- non_protected_species_areas %>%
  dplyr::mutate(non_protected_value = 1) %>%
  clean_function

### Summary value
protected_species <- oregon_hex %>%
  # join the leatherback sea turtle values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = leatherback_value,
                   by = "index") %>%
  # join the killer whale values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = killer_whale_value,
                   by = "index") %>%
  # join the humpback whale (Central America DPS) values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = humpback_ca_value,
                   by = "index") %>%
  # join the humpback whale (Mexico DPS) values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = humpback_mx_value,
                   by = "index") %>%
  # join the blue whale values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = bluewhale_value,
                   by = "index") %>%
  # select only fields of interest
  dplyr::select(index,
                leatherback_value,
                killerwhale_value,
                humpback_ca_value,
                humpback_mx_value,
                bluewhale_value) %>%
  
  # add value of 1 for datasets when hex cell has value of NA
  ## for hex cells not impacted by a particular dataset, that cell gets a value of 1
  ### this indicates  suitability with wind energy development
  dplyr::mutate(across(2:6, ~replace(x = .,
                                     list = is.na(.),
                                     # replacement values
                                     values = 1))) %>%

  # calculate the product of all protected species values
  dplyr:::mutate(species_product_value = leatherback_value *killerwhale_value * humpback_ca_value * humpback_mx_value * bluewhale_value) %>%
  # select all the key fields
  dplyr::select(index,
                leatherback_value,
                killerwhale_value,
                humpback_ca_value,
                humpback_mx_value,
                bluewhale_value,
                species_product_value) %>%
  as.data.frame()

#####################################

## Habitat layers
### EFHCA
efhca_value <- efhca %>%
  dplyr::mutate(efhca_value = 0.01) %>%
  clean_function()

### Rocky reef (mapped)
rreef_map_value <- rocky_reef_mapped %>%
  dplyr::mutate(rreef_map_value = 0.01) %>%
  clean_function()

### Rocky reef (probable)
rreef_prob_value <- rocky_reef_probable %>%
  dplyr::mutate(rreef_prob_value = 0.2) %>%
  clean_function()

### Deep sea coral and sponge
coral_sponge_value <- deep_sea_coral_sponge %>%
  dplyr::mutate(deep_coralsponge_value = z_value) %>%
  clean_function()

### Continental shelf break
continental_shelf_value <- continental_shelf %>%
  dplyr::mutate(continental_shelf_value = 0.6) %>%
  clean_function()

### Methane bubble streams
methane_bubble_value <- methane_bubble %>%
  dplyr::mutate(methane_bubble_value = 0.8) %>%
  clean_function()

### Summary value
habitat_values <- oregon_hex %>%
  # join the EFHCAs values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = efhca_value,
                   by = "index") %>%
  # join the rocky reef mapped values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = rreef_map_value,
                   by = "index") %>%
  # join the rocky reef probable values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = rreef_prob_value,
                   by = "index") %>%
  # join the deep-sea coral sponge values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = coral_sponge_value,
                   by = "index") %>%
  # join the contiental shelf values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = continental_shelf_value,
                   by = "index") %>%
  # join the methane bubble stream values by index field to full Oregon call area hex grid
  dplyr::left_join(x = .,
                   y = methane_bubble_value,
                   by = "index") %>%
  # select only fields of interest
  dplyr::select(index,
                efhca_value,
                rreef_map_value,
                rreef_prob_value,
                deep_coralsponge_value,
                continental_shelf_value,
                methane_bubble_value) %>%
  
  # add value of 1 for datasets when hex cell has value of NA
  ## for hex cells not impacted by a particular dataset, that cell gets a value of 1
  ### this indicates  suitability with wind energy development
  dplyr::mutate(across(2:7, ~replace(x = .,
                                     list = is.na(.),
                                     # replacement values
                                     values = 1))) %>%
  
  # calculate minimum value across the habitat subdatasets
  dplyr::mutate(habitat_value = pmin(efhca_value,
                                     rreef_map_value,
                                     rreef_prob_value,
                                     deep_coralsponge_value,
                                     continental_shelf_value,
                                     methane_bubble_value,
                                     # remove NA values from the minimum calculation
                                     na.rm = T)) %>%
  dplyr::select(index,
                efhca_value,
                rreef_map_value,
                rreef_prob_value,
                deep_coralsponge_value,
                continental_shelf_value,
                methane_bubble_value,
                habitat_value) %>%
  # convert to data frame to later join with Oregon hex grid
  as.data.frame()

#####################################

### Marine bird values
marine_bird_value <- marine_bird %>%
  clean_function() %>%
  dplyr::rename(marine_bird_value = marine_bird_index)

#####################################
#####################################

# Calculate geometric mean for natural resources submodel
oregon_natural_resources <- oregon_hex %>%
  dplyr::left_join(x = .,
                   y = protected_species,
                   by = "index") %>%
  dplyr::left_join(x = .,
                   y = habitat_values,
                   by = "index") %>%
 dplyr::left_join(x = .,
                  y = marine_bird_value,
                  by = "index") %>%

  # calculate the geometric mean
  ## geometric mean = nth root of the product of the variable values
  dplyr::mutate(nr_geom_mean = (species_product_value ^ nr_wt) * (habitat_value ^ nr_wt) * (marine_bird_value ^ nr_wt)) %>%
  
  # select the fields of interest
  dplyr::select(index,
                leatherback_value,
                killerwhale_value,
                humpback_ca_value,
                humpback_mx_value,
                bluewhale_value,
                species_product_value,
                efhca_value,
                rreef_map_value,
                rreef_prob_value,
                deep_coralsponge_value,
                continental_shelf_value,
                methane_bubble_value,
                habitat_value,
                marine_bird_value,
                nr_geom_mean) %>%
  # rename the geometry field
  dplyr::rename(geom = geom.x)

### Check to see if there are any duplicates of the indices
### There are none
nr_duplicates <- oregon_natural_resources %>%
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
sf::st_write(obj = oregon_natural_resources, dsn = oregon_suitability_gpkg, layer = paste0(region, "_", submodel, "_suitability"), append = F)

## Submodel
### Protected species
base::saveRDS(object = protected_species, file = paste(oregon_nr_species_dir, "oregon_protected_species_suitability.rds", sep = "/"))
base::saveRDS(object = leatherback_value, file = paste(oregon_nr_species_dir, "oregon_leatherback_value.rds", sep = "/"))
base::saveRDS(object = killer_whale_value, file = paste(oregon_nr_species_dir, "oregon_killer_whale_value.rds", sep = "/"))
base::saveRDS(object = humpback_mx_value, file = paste(oregon_nr_species_dir, "oregon_humpback_mexico_value.rds", sep = "/"))
base::saveRDS(object = humpback_ca_value, file = paste(oregon_nr_species_dir, "oregon_humpback_central_america_value.rds", sep = "/"))
base::saveRDS(object = bluewhale_value, file = paste(oregon_nr_species_dir, "oregon_blue_whale_value.rds", sep = "/"))

### Habitat layers
base::saveRDS(object = habitat_values, file = paste(oregon_nr_habitat_dir, layer = "oregon_habitat_suitability.rds", sep = "/"))
base::saveRDS(object = efhca_value, file = paste(oregon_nr_habitat_dir, layer = "oregon_efhca_value.rds", sep = "/"))
base::saveRDS(object = rreef_map_value, file = paste(oregon_nr_habitat_dir, layer = "oregon_rocky_reef_mapped_500m_value.rds", sep = "/"))
base::saveRDS(object = rreef_prob_value, file = paste(oregon_nr_habitat_dir, layer = "oregon_rocky_reef_probable_500m_value.rds", sep = "/"))
base::saveRDS(object = coral_sponge_value, file = paste(oregon_nr_habitat_dir, layer = "oregon_deep_sea_coral_sponge_value.rds", sep = "/"))
base::saveRDS(object = continental_shelf_value, file = paste(oregon_nr_habitat_dir, layer = "oregon_continental_shelf_value.rds", sep = "/"))
base::saveRDS(object = methane_bubble_value, file = paste(oregon_nr_habitat_dir, layer = "oregon_methane_bubble_streams_value.rds", sep = "/"))

### Seabird
base::saveRDS(object = marine_bird_value, file = paste(oregon_nr_marine_bird_dir, layer = "oregon_marine_bird_value.rds", sep = "/"))

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
