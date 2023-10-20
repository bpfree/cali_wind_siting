#######################################################################
### 31. Sensitivity Analysis -- Quantile Classifications (datasets) ###
#######################################################################

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

color_palette <- c(colorRampPalette(rev(brewer.pal(9, "BuPu")))(14),
                   # 0 values
                   "grey",
                   # positive values
                   colorRampPalette(brewer.pal(9, "YlGn"))(10))

color_palette <- c(rev(paletteer_c("ggthemes::Purple", 14)),
                   # 0 values
                   "grey",
                   # positive values
                   colorRampPalette(brewer.pal(9, "YlGn"))(10))

us_boundary <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf") %>%
  dplyr::filter(name %in% c("California",
                            "california"))

#####################################
#####################################

# load data
## Oregon suitability data
oregon_suitability <- sf::st_read(dsn = suitability_models, layer = "oregon_model_areas") %>%
  dplyr::select(index,
                # industry and operations
                sub_cable_value, sci_survey_value,
                # natural resources
                ## species product
                species_product_value,
                ## minimum habitat
                habitat_value,
                ## marine seabird
                marine_bird_value,
                # fisheries
                fisheries_value,
                # wind
                wind_value)

#####################################
#####################################

# Prepare clean quantitative classification dataframe for the datasets
oregon_dataset_quant_class <- oregon_suitability %>%
  dplyr::select(index)

# Calculate quantile classifications for each model iteration
for (i in 2:(length(oregon_suitability) - 1)){
  start2 <- Sys.time()
  
  # if need to test iterations
  #i <- 2
  
  # get the name of the dataset
  name_change <- names(oregon_suitability)[i] %>%
    # extract name of leftout dataset
    ## ?<= will look back and exclude the match (mean_)
    ## .*? will look for match until next pattern
    ## ?= will look ahead and exclude the match (_value)
    stringr::str_extract(string = ., pattern = ".*?(?=_value)")
  
  # calculate quantile
  quantile_calculation <- quantile(x = oregon_suitability[[i]],
                                   # probabilities between (0 - 1)
                                   ## these will be 0%, 12.5%, 25%, 37.5%, 50%, 67.5%, 75%, 87.5%, and 100%
                                   probs = c(0.0, 0.125, 0.25,
                                             0.375, 0.5,
                                             0.675, 0.75,
                                             0.875, 1.0))
  
  # reclassify quantiles
  quant_class_iteration <- oregon_suitability %>%
    # create new score classification and dictate which quantile a hex grid gets classified
    dplyr::mutate(!!paste(name_change, "quantile_classification", sep = "_") := dplyr::case_when(.[[i]] <= quantile_calculation[[1]] ~ 0,
                                                                                                 # 12.5% percentile
                                                                                                 .[[i]] > quantile_calculation[[1]] &
                                                                                                   .[[i]] <= quantile_calculation[[2]] ~ 1,
                                                                                                 # 25% percentile
                                                                                                 .[[i]] > quantile_calculation[[2]] & 
                                                                                                   .[[i]] <= quantile_calculation[[3]] ~ 2,
                                                                                                 # 37.5% percentile
                                                                                                 .[[i]] > quantile_calculation[[3]] & 
                                                                                                   .[[i]] <= quantile_calculation[[4]] ~ 3,
                                                                                                 # 50% percentile
                                                                                                 .[[i]] > quantile_calculation[[4]] & 
                                                                                                   .[[i]] <= quantile_calculation[[5]] ~ 4,
                                                                                                 # 67.5% percentile
                                                                                                 .[[i]] > quantile_calculation[[5]] & 
                                                                                                   .[[i]] <= quantile_calculation[[6]] ~ 5,
                                                                                                 # 75% percentile
                                                                                                 .[[i]] > quantile_calculation[[6]] & 
                                                                                                   .[[i]] <= quantile_calculation[[7]] ~ 6,
                                                                                                 # 87.5% percentile
                                                                                                 .[[i]] > quantile_calculation[[7]] & 
                                                                                                   .[[i]] <= quantile_calculation[[8]] ~ 7,
                                                                                                 # 100% percentile
                                                                                                 .[[i]] > quantile_calculation[[8]] & 
                                                                                                   .[[i]] <= quantile_calculation[[9]] ~ 8),
                  # move field after the dataset model geometric mean field
                  .before = geom) %>%
    
    # convert to dataframe so geometry is dropped and not duplicated when binded to the reference dataframe
    as.data.frame() %>%
    
    # select the fields created by the iteration
    dplyr::select(!!paste(name_change, "quantile_classification", sep = "_"))
  
  # add the results from the iteration to the template quantile classification dataframe
  oregon_dataset_quant_class <- cbind(oregon_dataset_quant_class, quant_class_iteration)
  
  # print how long it takes to calculate
  print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating and adding", name_change, "data to dataframe", sep = " "))
}

#####################################
#####################################

# Calculate quantile classifications change for each model iteration
oregon_dataset_quant_class_change <- oregon_dataset_quant_class %>%
  dplyr::select(index)

for (i in 2:(length(oregon_dataset_quant_class)-2)){
  start2 <- Sys.time()
  
  #i <- 2
  
  dataset1_name <- names(oregon_dataset_quant_class)[i] %>%
    # extract name of leftout dataset
    ## ?<= will look back and exclude the match (mean_)
    ## .*? will look for match until next pattern
    ## ?= will look ahead and exclude the match (_value)
    stringr::str_extract(string = ., pattern = ".*?(?=_quantile)")
  
  for (p in (i+1):(length(oregon_dataset_quant_class)-1)){
    
    #p <- i + 1
    
    dataset2_name <- names(oregon_dataset_quant_class)[p] %>%
      # extract name of leftout dataset
      ## ?<= will look back and exclude the match (mean_)
      ## .*? will look for match until next pattern
      ## ?= will look ahead and exclude the match (_value)
      stringr::str_extract(string = ., pattern = ".*?(?=_quantile)")
    
    quant_change_iteration <- oregon_dataset_quant_class %>%
      dplyr::mutate(!!paste(dataset1_name, dataset2_name, "quant_score_change", sep = "_") := (oregon_dataset_quant_class[[p]] - oregon_dataset_quant_class[[i]])) %>%
      
      # convert to dataframe so geometry is dropped and not duplicated when binded to the reference dataframe
      as.data.frame() %>%
      
      # select the fields created by the iteration
      dplyr::select(!!paste(dataset1_name, dataset2_name, "quant_score_change", sep = "_"))
    
    # add the results from the iteration to the template quantile classification change dataframe
    oregon_dataset_quant_class_change <- cbind(oregon_dataset_quant_class_change, quant_change_iteration)
    
    # stop when it gets to the final two datasets (no need to compare the final dataset against itself)
    #if (i == (length(oregon_dataset_quant_class)-2) & p == (length(oregon_dataset_quant_class)-1)) break
    
    # print how long it takes to calculate
    print(paste("Iteration", p, "takes", Sys.time() - start2, "minutes to complete creating and adding changes between", dataset1_name, "and", dataset2_name, "data to dataframe", sep = " "))
  }
  
  # print how long it takes to calculate
  print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating and adding changes between for all of", dataset1_name, "data to dataframe", sep = " "))
}

#####################################
#####################################

combination_names <- oregon_dataset_quant_class_change %>%
  as.data.frame() %>%
  dplyr::select(-index,
                -geom) %>%
  base::colnames(.) %>%
  # substitute nothing ("") for all the parts of the names that contain "_quant_score_change"
  base::gsub(pattern = "_quant_score_change",
             # replacement (i.e., "")
             replacement = "",
             # vector where the matches will be compared against
             x = .) %>%
  base::gsub(pattern = "sub_cable",
             replacement = "Submarine cable",
             x = .) %>%
  base::gsub(pattern = "sci_survey",
             replacement = "Scientific survey",
             x = .) %>%
  base::gsub(pattern = "species_product",
             replacement = "Species Product",
             x = .) %>%
  base::gsub(pattern = "habitat",
             replacement = "Habitat",
             x = .) %>%
  base::gsub(pattern = "marine_bird",
             replacement = "Marine seabird",
             x = .) %>%
  base::gsub(pattern = "fisheries",
             replacement = "Fisheries",
             x = .) %>%
  base::gsub(pattern = "wind",
             replacement = "Wind",
             x = .) %>%
  base::gsub(pattern = "_",
             replacement = " & ",
             x = .)

#####################################

data_combo_tbl_compare <- as.data.frame(matrix(nrow = length(combination_names),
                             ncol = 5)) %>%
  dplyr::rename(layer_combination = V1,
                no_change = V2,
                change = V3,
                total = V4,
                change_pct = V5)

oregon_combo_change <- oregon_dataset_quant_class_change %>%
  as.data.frame() %>%
  dplyr::select(-index,
                -geom) %>%
  as.matrix()

#####################################

# Initialize vectors to store occurrences
no_change_occurrences <- numeric((length(oregon_dataset_quant_class_change)-2))
change_occurrences <- numeric((length(oregon_dataset_quant_class_change)-2))

# Quantile change table
for (i in 1:(length(oregon_dataset_quant_class_change)-2)){
  #i <- 1
  
  # Calculate occurrences for each column
  ## retrieve column with data combination
  quant_combo <- oregon_combo_change[, i]
  
  ## get sum of all occurrences when quantitative classifications are equal (i.e., 0)
  no_change_occurrences[i] <- sum(quant_combo == 0)
  
  ## get sum of all values not equal to zero and add to the list
  change_occurrences[i] <- sum(quant_combo != 0)
  
  #####################################
  
  # Populate comparison table
  data_combo_tbl_compare[i,1] <- combination_names[[i]]
  
  data_combo_tbl_compare[i,2] <- no_change_occurrences[[i]]
  
  data_combo_tbl_compare[i,3] <- change_occurrences[[i]]
  
  data_combo_tbl_compare[i,4] <- sum(data_combo_tbl_compare[i, 2:3])
  
  data_combo_tbl_compare[i,5] <- round(data_combo_tbl_compare[i, 3] / data_combo_tbl_compare[i, 4], 3)
}

#####################################

n <- ggplot(data_combo_tbl_compare, aes(x = reorder(layer_combination, change_pct), y = change_pct)) +
  geom_bar(stat = "identity",
           fill = "lightblue") +
  coord_flip() +
  labs(x = "Dataset combination",
       y = "Change proportion",
       title = paste("Comparison of quantitative classification change in", str_to_title(region), "between two datasets", sep = " "), 
       subtitle = "Aggregated hexes with a changed quantitative classification",
       caption = "Calculation: Total changed hexes / Total hexes") +
  theme_bw() +
  plot_theme

print(n)

#####################################
#####################################

for (i in 2:(length(oregon_dataset_quant_class_change) - 1)){
  start2 <- Sys.time()
  
  # i <- 8
  
  values <- data.frame((unique(oregon_dataset_quant_class_change[[i]])))
  decrease <- length(values[values < 0])
  increase <- length(values[values > 0])
  
  #####################################
  
  dataset_name <- colnames(oregon_dataset_quant_class_change)[[i]] %>%
    # extract name of leftout dataset
    ## ?<= will look back and exclude the match (mean_)
    ## .*? will look for match until next pattern
    ## ?= will look ahead and exclude the match (_value)
    stringr::str_extract(string = ., pattern = ".*?(?=_quant)")
  
  #####################################
  
  dataset_comp_map <- ggplot() +
    geom_sf(data = oregon_dataset_quant_class_change,
            # map the minimum classification for the hex grid
            mapping = aes(fill = factor(oregon_dataset_quant_class_change[[i]])),
            # line width is 0 (i.e., no line)
            lwd = 0) +
    
    # add US state boundary data
    geom_sf(data = us_boundary,
            # fill with grey
            fill = "grey90",
            # line with a darker grey
            color = "grey30",
            # line type is dot-dash
            linetype = "dotdash",
            # transparency at 50%
            alpha = 0.5) +
    
    # rescale the fill colors between purple and red (stretch from 9 colors to 10)
    scale_fill_manual(values = c(rev(paletteer_c("ggthemes::Purple", decrease)),
                                 # 0 values
                                 "grey",
                                 # positive values
                                 colorRampPalette(brewer.pal(9, "YlGn"))(increase)),
                      # have legend name be: "Minimum quantile classification"
                      name = "Quantile\nclassification\nchange") +
    
    # alter labels
    labs(x = "",
         # no y-axis label
         y = "",
         # title
         title = paste("Quantile classification change between\n", dataset_name, "in", str_to_title(region), sep = " "),
         # subtitle
         subtitle = paste("Quantile classification change by hex grid between\n", dataset_name, sep = " ")) +
    
    # expand view so US coastal boundary is visible
    coord_sf(xlim = c(xmin = st_bbox(oregon_dataset_quant_class_change)$xmin,
                      # show more eastern longitudes
                      xmax = st_bbox(oregon_dataset_quant_class_change)$xmax+50000),
             # latitudes
             ylim = c(ymin = st_bbox(oregon_dataset_quant_class_change)$ymin,
                      ymax = st_bbox(oregon_dataset_quant_class_change)$ymax)) +
    
    # add label for Oregon
    geom_sf_text(data = us_boundary[2,],
                 mapping = aes(label = name),
                 # text size
                 size = 4,
                 # shift from the central location to the west
                 nudge_x = -295000,
                 # shift from the central location to the south
                 nudge_y = -140000) +
    theme_bw() +
    plot_theme
  
  #####################################
  
  # print(dataset_comp_map)
  
  ggplot2::ggsave(plot = dataset_comp_map, filename = file.path(figure_dir, paste(region, dataset_name, "comparison_map.tiff", sep = "_")),
                  width = 6, height = 8, units = "in", dpi = 600, compression = "lzw")
  
  #####################################
  
  # print how long it takes to calculate
  print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating maps between", dataset_name, sep = " "))
}

#####################################
#####################################

# Export data
## Tables
base::saveRDS(object = oregon_dataset_quant_class, file = paste(sensitivity_dir, paste(region, "combined_dataset_comparison_quantile_classifications.rds", sep = "_"), sep = "/"))
base::saveRDS(object = oregon_dataset_quant_class_change, file = paste(sensitivity_dir, paste(region, "combined_dataset_comparison_quantile_classifications_change.rds", sep = "_"), sep = "/"))
base::saveRDS(object = data_combo_tbl_compare, file = paste(sensitivity_dir, paste(region, "dataset_combination_comparison_table.rds", sep = "_"), sep = "/"))

## Map
ggplot2::ggsave(plot = n, filename = file.path(figure_dir, paste(region, "dataset_combination_comparison_table.tiff", sep = "_")),
                width = 10, height = 4.5, units = "in", dpi = 600, compression = "lzw")

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
