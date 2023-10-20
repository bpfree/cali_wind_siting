################################################
### 28. Sensitivity Analysis -- Change Table ###
################################################

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

### Figures
figure_dir <- "figure"

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

plot_theme <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.title.x = element_text(size = 10),
                    axis.text.x = element_text(size = 8),
                    axis.title.y = element_text(size = 10),
                    axis.text.y = element_text(size = 8),
                    plot.subtitle = element_text(size = 8),
                    plot.caption = element_text(size = 6,
                                                lineheight = 0.5),
                    legend.title = element_text(size = 8),
                    legend.text = element_text(size = 5))

#####################################
#####################################

# load data
sensitivity_jackknife <- sf::st_read(dsn = sensitivity_gpkg, layer = paste(region, "sensitivity", sep = "_"))

#####################################
#####################################

# set template table for populating with nominal and percent changes
# between the different jackknifed analyses
sensitivity_nominal_percent_change <- sensitivity_jackknife %>%
  dplyr::select(index,
                model_geom_mean)

# Calculate nominal and percent changes between model iterations
for (i in 1:24){
  start2 <- Sys.time()
  
  # if wanting to test a particular dataset
  #i <- 8
  
  # identify indexed field for the original model's score
  original_model_score <- 33
  
  # identify all jackknifed final scores
  ## final scores are every 5th column
  ## 4 submodels + 1 final model
  p <- original_model_score + i*5
  
  # grab the name of the dataset
  name_change <- names(sensitivity_jackknife)[p] %>%
    # extract name of omitted dataset (field name contains
    # the omitted dataset)
    ## ?<= will look back and exclude the match (mean_)
    ## .*? will look for match until next pattern
    ## ?= will look ahead and exclude the match (_value)
    stringr::str_extract(string = ., pattern = "(?<=mean_).*?(?=_value)")
  
  #####################################
  
  # calculate the nominal and percent changes between dataset removal
  change_iteration <- sensitivity_jackknife %>%
    # create nominal change (new score - original score)
    dplyr::mutate(!!paste("model_geom_mean", name_change, "nominal_change", sep = "_") := .[[p]] - .[[original_model_score]],
                  # move field after the dataset model geometric mean field
                  .after = paste0("model_geom_mean_", name_change, "_value")) %>%
    # create percent change field ((new score - old score) / old score * 100)
    dplyr::mutate(!!paste("model_geom_mean", name_change, "percent_change", sep = "_") := round(((.[[p]] - .[[original_model_score]]) / .[[original_model_score]]) * 100, 2),
                  # move field after the dataset model nominal change field
                  .after = paste("model_geom_mean", name_change, "nominal_change", sep = "_")) %>%
    
    #####################################
    
    # convert to dataframe so geometry is dropped and not duplicated when bound to the reference dataframe
    as.data.frame() %>%
    
    #####################################
    
    # select the fields created by the iteration
    dplyr::select(paste("model_geom_mean", name_change, "nominal_change", sep = "_"),
                  paste("model_geom_mean", name_change, "percent_change", sep = "_"))
  
    # add the results from the iteration to the sensitivity nominal and percent change dataframe
    sensitivity_nominal_percent_change <- cbind(sensitivity_nominal_percent_change, change_iteration)
    
    #####################################
  
    # print how long it takes to calculate
    print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating and adding", name_change, "data to dataframe", sep = " "))
}

#####################################
#####################################

# create reference table to populate with
# minimum and maximum values for nominal
# and percent changes
table <- data.frame(dataset = character(),
                    n_min = numeric(),
                    n_max = numeric(),
                    p_min = numeric(),
                    p_max = numeric())

# Calculate minimum and maximum changes for each dataset on model geometric mean
for (i in 0:23){
  start2 <- Sys.time()
  
  # if wanting to test a particular dataset
  #i <- 4
  
  # identify all jackknifed changed score fields (nominal and percent changes)
  ## nominal change
  nominal_data <- 3 + i*2
  
  ## percent change
  percent_data <- 4 + i*2
  
  #####################################
  
  # grab the name of the datasent
  name_change <- names(sensitivity_nominal_percent_change)[nominal_data] %>%
    # extract name of leftout dataset
    ## ?<= will look back and exclude the match (mean_)
    ## .*? will look for match until next pattern
    ## ?= will look ahead and exclude the match (_value)
    stringr::str_extract(string = ., pattern = "(?<=mean_).*?(?=_nominal_change)")
  
  # add the dataset layer name to the dataset list
  ## this will allow for each row to have minimums and maximums
  ## for nominal and percent changes by dataset
  dataset_list <- c(name_change)
  
  #####################################
  
  # calculate minimum and maximum nominal changes
  ## minimum
  n_min <- min(sensitivity_nominal_percent_change[[nominal_data]])
  ## maximum
  n_max <- max(sensitivity_nominal_percent_change[[nominal_data]])
  
  # calculate minimum and maximum percent changes
  ## minimum
  p_min <- min(sensitivity_nominal_percent_change[[percent_data]])
  ## maximum
  p_max <- max(sensitivity_nominal_percent_change[[percent_data]])
  
  #####################################
  
  # get resulting table with nominal and percent change minimums and maximums
  result_table <- data.frame(
    dataset_list,
    n_min,
    n_max,
    p_min,
    p_max)
  
  # bind the resulting table to the original table
  table <- rbind(table, result_table)
  
  #####################################
  
  # print how long it takes to calculate
  print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating and adding", name_change, "data to dataframe", sep = " "))
}

#####################################
#####################################

table <- table %>%
  # calculate the ranges for nominal changes and percent changes
  dplyr::mutate(n_range = n_max - n_min,
                p_range = p_max - p_min) %>%
  # recode the layer names
  dplyr::mutate(dataset_list = recode(dataset_list,
                                      "sc500" = "Submarine cable (500m)",
                                      "sc1000" = "Submarine cable (1000m)",
                                      "sub_cable" = "Submarine cable (combined)",
                                      "eastwest" = "East-west surveys",
                                      "eastwest_add" = "East-west surveys (additional)",
                                      "sstat" = "Survey station",
                                      "stransect" = "Survey transect",
                                      "sci_survey" = "Scientific survey (combined)",
                                      "leatherback" = "Leatherback sea turtle",
                                      "killerwhale" = "Killer whale",
                                      "humpback_ca" = "Humpback whale (CA)",
                                      "humpback_mx" = "Humpback whale (MX)",
                                      "bluewhale" = "Blue whale",
                                      "species_product" = "Species product (combined)",
                                      "efhca" = "EFHCA",
                                      "rreef_map" = "Rocky reef (mapped)",
                                      "rreef_prob" = "Rocky reef (probable)",
                                      "deep_coralsponge" = "Deep sea coral-sponge",
                                      "continental_shelf" = "Contiental shelf (10km)",
                                      "methane_bubble" = "Methane bubble streams",
                                      "habitat" = "Habitat (combined)",
                                      "marine_bird" = "Marine seabird",
                                      "fisheries" = "Fisheries",
                                      "wind" = "Wind"))

#####################################
#####################################

# Tornado figures
## Nominal change figure
n <- ggplot(table, aes(x = reorder(dataset_list, n_range), y = n_range)) +
  geom_errorbar(aes(ymin = n_min,
                    ymax = n_max), width = 0.5, color = "darkblue") +
  coord_flip() +
  labs(x = "Dataset",
       y = "Change (Nominal)",
       title = paste("Nominal change for score fluctuations in", str_to_title(region), sep = " "), 
       subtitle = "Range of minimum and maximum changes across all hex grid cells in study area",
       caption = "Calculation: New score - Old score") +
  theme_bw() +
  plot_theme

print(n)

## Percent change figure
p <- ggplot(table, aes(x = reorder(dataset_list, p_range), y = p_range)) +
  geom_errorbar(aes(ymin = p_min,
                    ymax = p_max), width = 0.5, color = "darkblue") +
  coord_flip() +
  labs(x = "Dataset",
       y = "Change (Percent)",
       title = paste("Percent change for score fluctuations in", str_to_title(region), sep = " "), 
       subtitle = "Range of minimum and maximum changes across all hex grid cells in study area",
       caption = "Calculation: ((New score - Old score) / Old Score) * 100") +
  scale_y_continuous(limits = c(-75, 500), breaks = seq(-100, 500, 50)) +
  theme_bw() +
  plot_theme

print(p)

#####################################
#####################################

# Export data
## Nominal / percent change
sf::st_write(obj = sensitivity_nominal_percent_change, dsn = sensitivity_gpkg, layer = paste(region, "sensitivity_nominal_percent_change", sep = "_"), append = F)

## Table
base::saveRDS(object = table, file = paste(sensitivity_dir, paste(region, "sensitivity_table.rds", sep = "_"), sep = "/"))

## Figures
ggsave(n, filename=file.path(figure_dir, paste(region, "nominal_change_tornado_plot.tiff", sep = "_")), width=6.5,
       height=4.5, units="in", dpi=600, compression = "lzw")
ggsave(p, filename=file.path(figure_dir, paste(region, "percent_change_tornado_plot.tiff", sep = "_")), width=6.5,
       height=4.5, units="in", dpi=600, compression = "lzw")

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
