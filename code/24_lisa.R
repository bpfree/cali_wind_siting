##############################################
### 24. Local Index of Spatial Association ###
##############################################

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
suitability_models <- "data/d_suitability_data/oregon_suitability_model.gpkg"

## Output directories
### LISA directory
lisa_gpkg <- "data/e_rank_data/oregon_lisa.gpkg"

#####################################
#####################################

# Inspect available layers and names within suitability model geopackage
sf::st_layers(dsn = suitability_models,
              do_count = F)

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

## lisa
layer <- "lisa"
classification <- "_highhigh"

## designate date
date <- format(Sys.time(), "%Y%m%d")

#####################################
#####################################

# Load data
## Oregon call area suitability areas
oregon_model_areas_sf <- sf::st_read(dsn = suitability_models, layer = "oregon_model_areas")

oregon_model_areas <- sf::st_read(dsn = suitability_models, layer = "oregon_model_areas") %>%
  as_Spatial(.)

#####################################
#####################################

# Calculate local indicator of spatial autocorrelation (LISA)
## ***NOTE: Multiple packages exist for calculating in R. rgeoda was
##          selected as it was developed by the team of Luc Anselin,
##          who developed the LISA methodology and GeoDA software.

## Create distance weights (8400 meters)
### Examine the minimum distance for the layer to have at least
### a single neighbor adjacent
dist_thres <- rgeoda::min_distthreshold(oregon_model_areas_sf)
dist_thres # every hex will have at least a single neighbor at a distance of 216.1691 meters

### In Oregon, it was previously determine that a distance of 8400 meters
### returned the ideal results
weights_8400 <- rgeoda::distance_weights(sf_obj = oregon_model_areas_sf,
                                         # set distance weight to be 8400 meters
                                         dist_thres = 8400,
                                         # set false to calculate distance in a Euclidean space
                                         is_arc = F,
                                         # set false to keep in metric system (units are in miles)
                                         is_mile = F)

## Create the LISA product
start <- Sys.time() # set start time
lisa <- rgeoda::local_moran(w = weights_8400, # weight is equal to distance weight of 8400 meters
                            # analyze using the final model geometric mean values
                            df = oregon_model_areas_sf["model_geom_mean"],
                            # run 9999 times to generate pseudo-p-values (default permutations is 999)
                            permutations = 9999,
                            # set cutoff to be 0.05 significance
                            significance_cutoff = 0.05,
                            # set seed so future runs can compare (default is 123456789)
                            seed = 561974)
print(Sys.time() - start) # print how long it takes to calculate the LISA results

#####################################

## Gather important fields
### p-values
oregon_lisa_pvalues <- lisa$p_vals %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("p_vals" = ".")

### cluster values
oregon_lisa_cvalues <- lisa$c_vals %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("c_vals" = ".")

### LISA values
oregon_lisa_lisa_values <- lisa$lisa_vals %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("lisa_vals" = ".")

### number of neighbors
oregon_lisa_neighbors <- lisa$nn_vals %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("num_neigh" = ".")

### labels
oregon_lisa_labels <- lisa$labels %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("labels" = ".") %>%
  # add field to join with cluster values
  ## Predefined values
  ### 0 Not significant
  ### 1 High-High
  ### 2 Low-Low
  ### 3 Low-High
  ### 4 High-Low
  ### 5 Undefined
  ### 6 Isolated
  dplyr::mutate(c_vals = c(0, 1, 2, 3, 4, 5, 6))

### colors
oregon_lisa_colors <- lisa$colors %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("colors" = ".") %>%
  # add field to join with cluster values
  ## Predefined values
  ### 0 #eeeeee (not significant)
  ### 1 #FF0000 (high-high)
  ### 2 #0000FF (low-low)
  ### 3 #a7adf9 (low-high)
  ### 4 #f4ada8 (high-low)
  ### 5 #464646 (undefined)
  ### 6 #999999 (isolated)
  dplyr::mutate(c_vals = c(0, 1, 2, 3, 4, 5, 6))

### joined cluster values with respective labels and colors
oregon_lisa_cvalue_labels_colors <- oregon_lisa_cvalues %>%
  # join cluster values with labels using cluster value
  dplyr::left_join(x = .,
                   # label dataset
                   y = oregon_lisa_labels,
                   # join field is cluster values
                   by = "c_vals") %>%
  # join cluster values with colors using cluster value
  dplyr::left_join(x = .,
                   # colors dataset
                   y = oregon_lisa_colors,
                   # join field is cluster values
                   by = "c_vals")

#####################################

# false discovery rate
oregon_fdr <- rgeoda::lisa_fdr(gda_lisa = lisa, current_p = 0.05)
oregon_bo <- rgeoda::lisa_bo(gda_lisa = lisa, current_p = 0.05)

#####################################

# Oregon hex grid joined with all new fields for LISA
oregon_hex_lisa <- oregon_model_areas_sf %>%
  cbind(oregon_lisa_pvalues,
        oregon_lisa_cvalue_labels_colors,
        oregon_lisa_lisa_values,
        oregon_lisa_neighbors)

#####################################
#####################################

combinations <- oregon_hex_lisa %>%
  dplyr::group_by(labels,
                  colors) %>%
  dplyr::summarise()
combinations

oregon_lisa_highhigh <- oregon_hex_lisa %>%
  dplyr::filter(labels == "High-High")

significance <- oregon_hex_lisa %>%
  dplyr::mutate(fdr = oregon_fdr) %>%
  dplyr::mutate(fdr_cluster = rgeoda::lisa_clusters(gda_lisa = lisa, cutoff = oregon_fdr))

oregon_new_labels <- lisa$labels %>%
  # convert to data frame to join to hex grid
  as.data.frame() %>%
  # rename field
  dplyr::rename("labels" = ".") %>%
  # add field to join with cluster values
  ## Predefined values
  ### 0 Not significant
  ### 1 High-High
  ### 2 Low-Low
  ### 3 Low-High
  ### 4 High-Low
  ### 5 Undefined
  ### 6 Isolated
  dplyr::mutate(fdr_cluster = c(0, 1, 2, 3, 4, 5, 6))

significance_new_label <- significance %>%
  dplyr::left_join(x = .,
                   y = oregon_new_labels,
                   by = "fdr_cluster")

#####################################
#####################################

label_colors <- c("#eeeeee", # Not significant
                  "#f4ada8", # High-Low
                  "#0000FF", # Low-Low
                  "#FF0000", # High-High
                  "#a7adf9") # Low-High

g <- ggplot2::ggplot() + 
  ggplot2::geom_sf(data = oregon_lisa_highhigh, aes(fill = "#FF0000"), color = NA)
g

#####################################
#####################################

# Export data
## LISA
sf::st_write(obj = oregon_lisa_highhigh, dsn = lisa_gpkg, layer = paste(region, layer, classification, sep = "_"), append = F)
sf::st_write(obj = oregon_hex_lisa, dsn = lisa_gpkg, layer = paste(region, "hex", layer, sep = "_"), append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
