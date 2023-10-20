#############################################
### 27. Sensitivity Analysis -- Jackknife ###
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
## Input directories
suitability_models <- "data/d_suitability_data/suitability_model.gpkg"

## Output directories
### Sensitivity directory
sensitivity_gpkg <- "data/f_sensitivity_data/oregon_sensitivity.gpkg"
sensitivity_dir <- "data/f_sensitivity_data"

#####################################

# Inspect available layers and names within suitability model geopackage
sf::st_layers(dsn = suitability_models,
              do_count = T)

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

## designate date
date <- format(Sys.time(), "%Y%m%d")

#####################################
#####################################

# Function to calculate product of non-NA values
calculate_geometric_mean <- function(row) {
  # subset all rows that do not have NA values
  non_na_values <- row[!is.na(row)]
  # when that subset is longer than 0 (i.e., 1 or more columns
  # in that row)
  if (length(non_na_values) > 0) {
    # calculate product of all non-NA values across the row
    prod_values <- prod(non_na_values)
    # calculate geometric mean (value products raised to 1
    # over the number of fields with non-NA values)
    geometric_mean <- prod_values^(1/length(non_na_values))
    # return with the geometric mean value
    return(geometric_mean)
  } else {
    # when the above conditions are not met, so all values
    # are NA, then also return NA values
    return(NA)
  }
}

# Function to calculate product of non-NA values (combined species data)
calculate_species_product <- function(...) {
  # the ... signifies that there can take multiple arguments
  # that will get passed to the object
  values <- c(...)
  # get all the non-NA values from the data
  non_na_values <- values[!is.na(values)]
  # when that subset is longer than 0 (i.e., 1 or more columns
  # in that row)
  if (length(non_na_values) > 0) {
    # calculate product of all non-NA values across the row 
    return(prod(non_na_values))
  } else {
    # when the above conditions are not met, so all values
    # are NA, then also return NA values
    return(NA)
  }
}

#####################################
#####################################

# load data
## Oregon suitability data
oregon_suitability <- sf::st_read(dsn = suitability_models, layer = "oregon_model_areas") %>%
  dplyr::select(index,
                # constraints
                dod_value, pacpars_value,
                # industry and operations
                sc500_value, sc1000_value, sub_cable_value, eastwest_value, eastwest_add_value, sstat_value, stransect_value, sci_survey_value,
                # natural resources
                ## species product
                leatherback_value, killerwhale_value, humpback_ca_value, humpback_mx_value, bluewhale_value, species_product_value,
                ## minimum habitat
                efhca_value, rreef_map_value, rreef_prob_value, deep_coralsponge_value, continental_shelf_value,
                methane_bubble_value, habitat_value,
                ## marine seabird
                marine_bird_value,
                # fisheries
                fisheries_value,
                # wind
                wind_value,
                # submodel geometric values
                constraints, io_geom_mean, nr_geom_mean, fish_geom_mean, wind_geom_mean,
                # model geometric value
                model_geom_mean)

#####################################
#####################################

# Create a reference dataframe where data will get added
sensitivity_jackknife <- oregon_suitability

#####################################
#####################################

# Jackknife analysis across all datasets to determine how suitability scores change when dataset is removed
## Datasets cross fields 4 - 27
for (i in 4:27){
  start2 <- Sys.time()
  
  # if wanting to test a particular dataset
  #i <- 14
  
  name <- names(oregon_suitability)[i]
  
  sensitivity_iteration <- oregon_suitability %>%
    # convert to dataframe so calculations can get added to reference dataframe
    as.data.frame() %>%
    
    # when field is selected (column i) fill with NA values so as to "remove" it from analysis
    dplyr::mutate(across(.cols = i,
                         # for the column, when the data are not already NA, replace with NA values
                         .fns = ~replace(i, !is.na(i), NA))) %>%
    
    # industry and operations
    ### fill sub_cable_values with the sum of 500-m buffer and 500-1000-m buffer values
    ### when a hex has both values give it the lower value (500-m buffer = 0.6)
    #### when iterating, if the column is the sub_cable_value (i = 6), the field gets NA
    #### otherwise, follow the minimum value gets precedence
    dplyr::mutate(sub_cable_value = rowSums((.[,4:5]),
                                            # remove any NA values to ensure sum is correct
                                            na.rm = T)) %>%
    dplyr::mutate(sub_cable_value = case_when(i == 6 ~ NA,
                                              # when sum is equal to 2 (both buffers have
                                              # values as one, change back to 1)
                                              sub_cable_value == 2 ~ 1,
                                              # when sum is equal to 1.8 (500m buffer is
                                              # equal to 1 and 1000m buffer is equal to
                                              # 0.8, hex grid gets value of 0.8)
                                              sub_cable_value == 1.8 ~ 0.8,
                                              # when sum is equal to 1.6 (500m buffer is
                                              # equal to 0.6 and 1000m buffer is equal to
                                              # 1, hex grid gets value of 0.6)
                                              sub_cable_value == 1.6 ~ 0.6,
                                              # when sum is equal to 1.4 (500m buffer is
                                              # equal to 0.6 and 1000m buffer is equal to
                                              # 0.8, hex grid gets value of 0.6)
                                              sub_cable_value == 1.4 ~ 0.6,
                                              sub_cable_value == 1.0 ~ 1.0,
                                              # when sum is equal to 0.8 (500m buffer is
                                              # equal to 0 and 1000m buffer is equal to
                                              # 0.8, hex grid gets value of 0.8)
                                              sub_cable_value == 0.8 ~ 0.8,
                                              # when sum is equal to 0.8 (500m buffer is
                                              # equal to 0.6 and 1000m buffer is equal to
                                              # 0, hex grid gets value of 0.6)
                                              sub_cable_value == 0.6 ~ 0.6)) %>%       
    
    ## calculate a summary value for scientific surveys
    ### when a hex grid cell has multiple values the minimum
    ### value across the datasets is assigned to the new
    ### summarized field
    dplyr::mutate(sci_survey_value = case_when(i == 11 ~ NA,
                                               # in all other cases when i is not the field for combined
                                               # scientific surveys, it gets the minimum value across all
                                               # scientific survey values
                                               i != 11 ~ pmin(eastwest_value,
                                                              eastwest_add_value,
                                                              sstat_value,
                                                              stransect_value,
                                                              # remove any values that are NA when new field
                                                              na.rm = T))) %>%
                                          
    ## calculate the geometric mean (geometric mean = nth root of the product of the variable values)
    dplyr::mutate(!!paste("io_geom_mean", name, sep ="_") := apply(X = .[c("sub_cable_value",
                                                                           "sci_survey_value")],
                                                                   # 1 = rows, 2 = columns
                                                                   MARGIN = 1,
                                                                   # use the calculate geometric mean function
                                                                   FUN = calculate_geometric_mean)) %>%
    
    #####################################
  
    # natural resources
    ## calculate the product of all protected species values
    dplyr::mutate(species_product_value = case_when(i == 17 ~ NA,
                                                    # in all other cases when i is not the field for combined
                                                    # species, it gets the product of all the species that do
                                                    # not have NA values
                                                    i != 17 ~ mapply(FUN = calculate_species_product,
                                                                     # fields to apply function to
                                                                     leatherback_value,
                                                                     killerwhale_value,
                                                                     humpback_ca_value,
                                                                     humpback_mx_value,
                                                                     bluewhale_value))) %>%
    
    ## calculate minimum value across the habitat subdatasets
    dplyr::mutate(habitat_value = case_when(i == 24 ~ NA,
                                            # in all other cases when i is not the field for combined
                                            # habitat value, it gets the minimum value across the
                                            # habitat fields
                                            i != 24 ~ pmin(efhca_value,
                                                           rreef_map_value,
                                                           rreef_prob_value,
                                                           deep_coralsponge_value,
                                                           continental_shelf_value,
                                                           methane_bubble_value,
                                                           # remove NA values from the minimum calculation
                                                           na.rm = T))) %>%
                                        
    ## calculate the geometric mean (geometric mean = nth root of the product of the variable values)
    dplyr::mutate(!!paste("nr_geom_mean", name, sep= "_") := apply(X = .[c("species_product_value",
                                                                           "habitat_value",
                                                                           "marine_bird_value")],
                                                                   # 1 = rows, 2 = columns
                                                                   MARGIN = 1,
                                                                   # use the calculate geometric mean function
                                                                   FUN = calculate_geometric_mean)) %>%
    
    #####################################
  
    # fisheries
    ## calculate the geometric mean (geometric mean = nth root of the product of the variable values)
    dplyr::mutate(!!paste("fish_geom_mean", name, sep = "_") := apply(X = .[c("fisheries_value")],
                                                                      # 1 = rows, 2 = columns
                                                                      MARGIN = 1,
                                                                      # use the calculate geometric mean function
                                                                      FUN = calculate_geometric_mean)) %>%
    
    #####################################
  
    # wind
    ## calculate the geometric mean (geometric mean = nth root of the product of the variable values)
    dplyr::mutate(!!paste("wind_geom_mean", name, sep = "_") := apply(X = .[c("wind_value")],
                                                                      # 1 = rows, 2 = columns
                                                                      MARGIN = 1,
                                                                      # use the calculate geometric mean function
                                                                      FUN = calculate_geometric_mean)) %>%
    
    #####################################
  
    # final model score
    ## calculate the geometric mean (geometric mean = nth root of the product of the variable values)
    dplyr::mutate(!!paste("model_geom_mean", name, sep = "_") := apply(X = .[c(paste("io_geom_mean", name, sep = "_"),
                                                                               paste("nr_geom_mean", name, sep = "_"),
                                                                               paste("fish_geom_mean", name, sep = "_"),
                                                                               paste("wind_geom_mean", name, sep = "_"))],
                                                                       # 1 = rows, 2 = columns
                                                                       MARGIN = 1,
                                                                       # use the calculate geometric mean function
                                                                       FUN = calculate_geometric_mean)) %>%
    
    #####################################
    
    # convert to dataframe so geometry is dropped and not duplicated when bound to the reference dataframe
    as.data.frame() %>%
    
    #####################################
    
    # select the fields created by the iteration
    dplyr::select(paste("io_geom_mean", name, sep = "_"),
                  paste("nr_geom_mean", name, sep = "_"),
                  paste("fish_geom_mean", name, sep = "_"),
                  paste("wind_geom_mean", name, sep = "_"),
                  paste("model_geom_mean", name, sep = "_"))
  
    # add the results from the iteration to the sensitivity jackknife dataframe
    sensitivity_jackknife <- cbind(sensitivity_jackknife, sensitivity_iteration)
    
    #####################################
  
    # print how long it takes to calculate
    print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating and adding", name, "data to dataframe", sep = " "))
}

#####################################
#####################################

# Export data
## Jackknife
sf::st_write(obj = sensitivity_jackknife, dsn = sensitivity_gpkg, layer = paste(region, "sensitivity", sep = "_"), append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
