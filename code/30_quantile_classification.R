###############################################################################
### 30. Sensitivity Analysis -- Quantile Classifications (model iterations) ###
###############################################################################

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

## data layer names
data_names_list <- list("Submarine cable (500m)",
                        "Submarine cable (1000m)",
                        "Submarine cable (combined)",
                        "East-west surveys",
                        "East-west surveys (additional)",
                        "Survey station",
                        "Survey transect",
                        "Scientific survey (combined)",
                        "Leatherback sea turtle",
                        "Killer whale",
                        "Humpback whale (CA)",
                        "Humpback whale (MX)",
                        "Blue whale",
                        "Species product (combined)",
                        "EFHCA",
                        "Rocky reef (mapped)",
                        "Rocky reef (probable)",
                        "Deep sea coral-sponge",
                        "Contiental shelf (10km)",
                        "Methane bubble streams",
                        "Habitat (combined)",
                        "Marine seabird",
                        "Fisheries",
                        "Wind")

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

# color_palette <- c(colorRampPalette(rev(brewer.pal(9, "BuPu")))(14),
#                    # 0 values
#                    "grey",
#                    # positive values
#                    colorRampPalette(brewer.pal(9, "YlGn"))(10))

color_palette <- c(rev(paletteer_c("ggthemes::Purple", 14)),
                   # 0 values
                   "grey",
                   # positive values
                   colorRampPalette(brewer.pal(9, "YlGn"))(10))

iteration_change_palette <- c(# 0 values
                   "grey",
                   # positive values
                   colorRampPalette(brewer.pal(9, "YlOrRd"))(15))

us_boundary <- rnaturalearth::ne_states(country = "United States of America", returnclass = "sf") %>%
  dplyr::filter(name %in% c("California",
                            "california"))

#####################################
#####################################

# load data
sensitivity_jackknife <- sf::st_read(dsn = sensitivity_gpkg, layer = paste(region, "sensitivity", sep = "_"))

#####################################
#####################################

# Calculate sample quantiles
## Template sensitivity quantile classification
sensitivity_quant_class <- sensitivity_jackknife %>%
  # select fields of interest
  dplyr::select(index,
                model_geom_mean)

### Calculate the quantile probabilities for 8 classifications
model_quant_class <- quantile(x = sensitivity_quant_class[["model_geom_mean"]],
                              # probabilities between (0 - 1)
                              ## these will be 0%, 12.5%, 25%, 37.5%, 50%, 67.5%, 75%, 87.5%, and 100%
                              probs = c(0.0, 0.125, 0.25,
                                        0.375, 0.5,
                                        0.675, 0.75,
                                        0.875, 1.0))

### Reclassification of probabilities too classifications (1 - 8)
sensitivity_quant_class <- sensitivity_quant_class %>%
  # recode quantile scores
  dplyr::mutate(model_score_classification = dplyr::case_when(.[[2]] <= model_quant_class[[1]] ~ 0,
                                                              # when value is above the 0% quantile and below the 12.5% quantile,
                                                              # give it a value of 1
                                                              .[[2]] > model_quant_class[[1]] & 
                                                                .[[2]] <= model_quant_class[[2]] ~ 1,
                                                              # when value is above the 12.5% quantile and below the 25% quantile,
                                                              # give it a value of 2
                                                              .[[2]] > model_quant_class[[2]] & 
                                                                .[[2]] <= model_quant_class[[3]] ~ 2,
                                                              # when value is above the 25% quantile and below the 37.5% quantile,
                                                              # give it a value of 3
                                                              .[[2]] > model_quant_class[[3]] & 
                                                                .[[2]] <= model_quant_class[[4]] ~ 3,
                                                              # when value is above the 37.5% quantile and below the 50% quantile,
                                                              # give it a value of 4
                                                              .[[2]] > model_quant_class[[4]] & 
                                                                .[[2]] <= model_quant_class[[5]] ~ 4,
                                                              # when value is above the 50% quantile and below the 67.5% quantile,
                                                              # give it a value of 5
                                                              .[[2]] > model_quant_class[[5]] & 
                                                                .[[2]] <= model_quant_class[[6]] ~ 5,
                                                              # when value is above the 67.5% quantile and below the 75% quantile,
                                                              # give it a value of 6
                                                              .[[2]] > model_quant_class[[6]] & 
                                                                .[[2]] <= model_quant_class[[7]] ~ 6,
                                                              # when value is above the 75% quantile and below the 87.5% quantile,
                                                              # give it a value of 7
                                                              .[[2]] > model_quant_class[[7]] & 
                                                                .[[2]] <= model_quant_class[[8]] ~ 7,
                                                              # when value is above the 87.5% quantile and below the 100% quantile,
                                                              # give it a value of 8
                                                              .[[2]] > model_quant_class[[8]] & 
                                                                .[[2]] <= model_quant_class[[9]] ~ 8),
                # move classification field after the model geometric mean field
                .after = model_geom_mean) %>%
  # remove the model geometric mean field
  dplyr::select(-model_geom_mean)

#####################################
#####################################

# Calculate quantile classifications for each model iteration
for (i in 1:24){
  start2 <- Sys.time()
  
  # if wanting to test a particular dataset
  #i <- 1
  
  # iterate over fields (p = parameter) of final geometric mean scores
  p <- 33 + i*5
  
  # get the name of the dataset
  name_change <- names(sensitivity_jackknife)[p] %>%
    # extract name of leftout dataset
    ## ?<= will look back and exclude the match (mean_)
    ## .*? will look for match until next pattern
    ## ?= will look ahead and exclude the match (_value)
    stringr::str_extract(string = ., pattern = "(?<=mean_).*?(?=_value)")
  
  # calculate quantile
  quantile_calculation <- quantile(x = sensitivity_jackknife[[p]],
                                   # probabilities between (0 - 1)
                                   ## these will be 0%, 12.5%, 25%, 37.5%, 50%, 67.5%, 75%, 87.5%, and 100%
                                   probs = c(0.0, 0.125, 0.25,
                                             0.375, 0.5,
                                             0.675, 0.75,
                                             0.875, 1.0))
  
  # reclassify quantiles
  quant_class_iteration <- sensitivity_jackknife %>%
    # create new score classification and dictate which quantile a hex grid gets classified
    dplyr::mutate(!!paste("model_geom_mean", name_change, "score_classification", sep = "_") := dplyr::case_when(.[[p]] <= quantile_calculation[[1]] ~ 0,
                                                                                                                 # when value is above the 0% quantile and below the 12.5% quantile,
                                                                                                                 # give it a value of 1
                                                                                                                 .[[p]] > quantile_calculation[[1]] &
                                                                                                                   .[[p]] <= quantile_calculation[[2]] ~ 1,
                                                                                                                 # when value is above the 12.5% quantile and below the 25% quantile,
                                                                                                                 # give it a value of 2
                                                                                                                 .[[p]] > quantile_calculation[[2]] & 
                                                                                                                   .[[p]] <= quantile_calculation[[3]] ~ 2,
                                                                                                                 # when value is above the 25% quantile and below the 37.5% quantile,
                                                                                                                 # give it a value of 3
                                                                                                                 .[[p]] > quantile_calculation[[3]] & 
                                                                                                                   .[[p]] <= quantile_calculation[[4]] ~ 3,
                                                                                                                 # when value is above the 37.5% quantile and below the 50% quantile,
                                                                                                                 # give it a value of 4
                                                                                                                 .[[p]] > quantile_calculation[[4]] & 
                                                                                                                   .[[p]] <= quantile_calculation[[5]] ~ 4,
                                                                                                                 # when value is above the 50% quantile and below the 67.5% quantile,
                                                                                                                 # give it a value of 5
                                                                                                                 .[[p]] > quantile_calculation[[5]] & 
                                                                                                                   .[[p]] <= quantile_calculation[[6]] ~ 5,
                                                                                                                 # when value is above the 67.5% quantile and below the 75% quantile,
                                                                                                                 # give it a value of 6
                                                                                                                 .[[p]] > quantile_calculation[[6]] & 
                                                                                                                   .[[p]] <= quantile_calculation[[7]] ~ 6,
                                                                                                                 # when value is above the 75% quantile and below the 87.5% quantile,
                                                                                                                 # give it a value of 7
                                                                                                                 .[[p]] > quantile_calculation[[7]] & 
                                                                                                                   .[[p]] <= quantile_calculation[[8]] ~ 7,
                                                                                                                 # when value is above the 87.5% quantile and below the 100% quantile,
                                                                                                                 # give it a value of 8
                                                                                                                 .[[p]] > quantile_calculation[[8]] & 
                                                                                                                   .[[p]] <= quantile_calculation[[9]] ~ 8),
                  # move field after the dataset model geometric mean field
                  .after = paste0("model_geom_mean_", name_change, "_value")) %>%
    
    # convert to dataframe so geometry is dropped and not duplicated when binded to the reference dataframe
    as.data.frame() %>%
    
    # select the fields created by the iteration
    dplyr::select(!!paste("model_geom_mean", name_change, "score_classification", sep = "_"))
  
  # add the results from the iteration to the template quantile classification dataframe
  sensitivity_quant_class <- cbind(sensitivity_quant_class, quant_class_iteration)
  
  #####################################
  
  # print how long it takes to calculate
  print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating and adding", name_change, "data to dataframe", sep = " "))
}

#####################################
#####################################

# Calculate quantile classifications change for each model iteration
sensitivity_quant_change <- sensitivity_quant_class %>%
  dplyr::select(index)

for (i in 3:26){
  start2 <- Sys.time()
  
  #i <- 26
  
  name_change <- names(sensitivity_quant_class)[i] %>%
    # extract name of leftout dataset
    ## ?<= will look back and exclude the match (mean_)
    ## .*? will look for match until next pattern
    ## ?= will look ahead and exclude the match (_value)
    stringr::str_extract(string = ., pattern = "(?<=mean_).*?(?=_score)")
  
  #####################################
  
  quant_change_iteration <- sensitivity_quant_class %>%
    dplyr::mutate(!!paste("model_geom_mean", name_change, "quant_score_change", sep = "_") := (sensitivity_quant_class[[i]] - sensitivity_quant_class[[2]])) %>%
    
    # convert to dataframe so geometry is dropped and not duplicated when bound to the reference dataframe
    as.data.frame() %>%
    
    # select the fields created by the iteration
    dplyr::select(!!paste("model_geom_mean", name_change, "quant_score_change", sep = "_"))
    
  # add the results from the iteration to the template quantile classification change dataframe
  sensitivity_quant_change <- cbind(sensitivity_quant_change, quant_change_iteration)
  
  #####################################
  
  # print how long it takes to calculate
  print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating and adding", name_change, "data to dataframe", sep = " "))
}

#####################################

# Create total change across all datasets
sensitivity_quant_total_iteration_change <- as.data.frame(sensitivity_quant_change) %>%
  # calculate the total change and iteration change
  dplyr::mutate(total_change = base::rowSums(x = .[2:25]),
                # total number of times a hex grid changed classification in any iteration
                # no change = 0, so anytime change does not equal 0
                iteration_change = base::rowSums(x = .[2:25] != 0)) %>%
  # select fields of interest
  dplyr::select(index,
                total_change,
                iteration_change,
                geom) %>%
  # set back as sf
  sf::st_as_sf(x = .)

## Create histogram of total change
tc_hist <- ggplot() +
  geom_histogram(data = sensitivity_quant_total_iteration_change,
                 aes(x = total_change),
                 # fill color of histogram
                 fill = "#869BC4",
                 # line color of histogram
                 color = "#869BC4",
                 # transparency of 80%
                 alpha = 0.8,
                 # bin size is 25
                 bins = 25) +
  labs(x = "Quantile change",
       y = "Occurances",
       title = paste("Total quantile change in", str_to_title(region), sep = " "), 
       subtitle = "Aggregated quantile classification change by hex grid") +
  theme_bw() +
  plot_theme

tc_hist

## Create map of total change
tc_plot <- ggplot() +
  # add the total change data
  geom_sf(data = sensitivity_quant_total_iteration_change,
          # what will be mapped is the value of total change for each hex grid
          mapping = aes(fill = factor(total_change)),
          # line width is 0 (i.e., no line)
          lwd = 0) +
  # add US state boundary data
  geom_sf(data = us_boundary,
          # fill with grey
          fill = "grey90",
          color = "grey30",
          linetype = "dotdash",
          alpha = 0.5) +
  # fill the total change by the pre-defined color palette
  scale_fill_manual(values = color_palette,
                    name = "Total change") +
  labs(x = "",
       y = "",
       title = paste("Total quantile change in", str_to_title(region), sep = " "),
       subtitle = "Aggregated quantile classification change by hex grid") +
  # change the coordinate box so it goes beyond just the offshore wind farms
  coord_sf(xlim = c(xmin = st_bbox(sensitivity_quant_total_iteration_change)$xmin,
                    xmax = st_bbox(sensitivity_quant_total_iteration_change)$xmax+50000),
           # latitudes
           ylim = c(ymin = st_bbox(sensitivity_quant_total_iteration_change)$ymin,
                    ymax = st_bbox(sensitivity_quant_total_iteration_change)$ymax)) +
  # have legend get two columns
  guides(fill = guide_legend(ncol = 2)) +
  geom_sf_text(data = us_boundary[2,], mapping = aes(label = name)) +
  theme_bw() +
  plot_theme

print(tc_plot)

#####################################

## Create histogram of total change
ic_hist <- ggplot() +
  geom_histogram(data = sensitivity_quant_total_iteration_change,
                 aes(x = iteration_change),
                 # fill color of histogram
                 fill = "#869BC4",
                 # line color of histogram
                 color = "#869BC4",
                 # transparency of 80%
                 alpha = 0.8,
                 # bin size is 25
                 bins = 25) +
  labs(x = "Quantile change",
       y = "Occurances",
       title = paste("Iteration quantile change in", str_to_title(region), sep = " "), 
       subtitle = "Aggregated quantile classification change by hex grid") +
  theme_bw() +
  plot_theme

ic_hist

## Create map of iteration change
ic_plot <- ggplot() +
  # add the iteration change data
  geom_sf(data = sensitivity_quant_total_iteration_change,
          # what will be mapped is the value of iteration change for each hex grid
          mapping = aes(fill = factor(iteration_change)),
          # line width is 0 (i.e., no line)
          lwd = 0) +
  # add US state boundary data
  geom_sf(data = us_boundary,
          # fill with grey
          fill = "grey90",
          color = "grey30",
          linetype = "dotdash",
          alpha = 0.5) +
  # fill the iteration change by the pre-defined color palette
  scale_fill_manual(values = iteration_change_palette,
                    name = "Iteration change") +
  labs(x = "",
       y = "",
       title = paste("Iteration quantile change in", str_to_title(region), sep = " "),
       subtitle = "Aggregated quantile classification change by hex grid") +
  # change the coordinate box so it goes beyond just the offshore wind farms
  coord_sf(xlim = c(xmin = st_bbox(sensitivity_quant_total_iteration_change)$xmin,
                    xmax = st_bbox(sensitivity_quant_total_iteration_change)$xmax+50000),
           # latitudes
           ylim = c(ymin = st_bbox(sensitivity_quant_total_iteration_change)$ymin,
                    ymax = st_bbox(sensitivity_quant_total_iteration_change)$ymax)) +
  # have legend get two columns
  guides(fill = guide_legend(ncol = 2)) +
  geom_sf_text(data = us_boundary[2,], mapping = aes(label = name)) +
  theme_bw() +
  plot_theme

print(ic_plot)

#####################################
#####################################

# minimum quantitative classification for each index
sensitivity_quant_minmax <- sensitivity_quant_class %>%
  dplyr::mutate(quant_min = purrr::pmap(.l = across(starts_with("model")),
                                        .f = min),
                quant_max = purrr::pmap(.l = across(starts_with("model")),
                                        .f = max)) %>%
  dplyr::select(index,
                quant_min,
                quant_max)

sensitivity_quant_minmax <- as.data.frame(lapply(sensitivity_quant_minmax, unlist)) %>%
  dplyr::group_by(index) %>%
  dplyr::reframe(quant_min = min(quant_min),
                 quant_max = max(quant_max)) %>%
  dplyr::left_join(x = .,
                   y = sensitivity_quant_total_iteration_change,
                   by = "index") %>%
  st_as_sf( x = .)

#####################################

min_plot <- ggplot() +
  geom_sf(data = sensitivity_quant_minmax,
          # map the minimum classification for the hex grid
          mapping = aes(fill = factor(quant_min)),
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
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "PuRd"))(10),
                    # have legend name be: "Minimum quantile classification"
                    name = "Minimum\nquantile\nclassification") +
  # alter labels
  labs(x = "",
       # no y-axis label
       y = "",
       # title
       title = paste("Minimum quantile classification in", str_to_title(region), sep = " "),
       # subtitle
       subtitle = "Minimum quantile classification change by hex grid") +
  # have legend be across 3 columns
  #guides(fill = guide_legend(ncol = 3)) +
  # expand view so US coastal boundary is visible
  coord_sf(xlim = c(xmin = st_bbox(sensitivity_quant_minmax)$xmin,
                    # show more eastern longitudes
                    xmax = st_bbox(sensitivity_quant_minmax)$xmax+50000),
           # latitudes
           ylim = c(ymin = st_bbox(sensitivity_quant_minmax)$ymin,
                    ymax = st_bbox(sensitivity_quant_minmax)$ymax)) +
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

print(min_plot)

#####################################
#####################################

list <- list()

# Contingency matrix
## create dataframes and list of 
for (i in 3:(length(sensitivity_quant_class) - 1)){
  start2 <- Sys.time()
  
  # i <- 6
  
  original_data <- as.data.frame(sensitivity_quant_class[, c(2,i)]) %>%
    dplyr::select(-geom)
  
  #####################################
  
  # get the dataset name
  layer_name <- names(original_data)[2] %>%
    # extract name of leftout dataset
    ## ?<= will look back and exclude the match (mean_)
    ## .*? will look for match until next pattern
    ## ?= will look ahead and exclude the match (_value)
    stringr::str_extract(string = ., pattern = "(?<=mean_).*?(?=_score)")
  
  matrix_name <- paste("contingency_matrix", layer_name, sep = "_")
  
  #####################################
  
  # create a contigency matrix for the desired dataset compared to the original model
  contigency_matrix <- assign(matrix_name, as.data.frame(unclass(table(original_data))))
  
  # add each contigency matrix to the reference list
  list[[length(list) + 1]] <- contigency_matrix
  
  #####################################
  
  # save each dataset's contigency matrix as an RDS file
  base::saveRDS(object = contigency_matrix,
                file = paste(sensitivity_dir, paste(region, layer_name, "consistency_matrix.rds", sep = "_"), sep = "/"))

  #####################################
  
  # print how long it takes to calculate and export as RDS file
  print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating and adding", layer_name, "data to dataframe and export as RDS file", sep = " "))
}

## create plot of consistency matrices for each dataset
for (i in 1:length(list)){
  start2 <- Sys.time()
  
  #i <- 1
  
  layer_name <- data_names_list[[i]]
  
  #####################################
  
  contigency_matrix_species <- list[[i]] %>%
    dplyr::mutate(class = row_number() - 1) %>%
    dplyr::relocate(class, .before = 0) %>%
    reshape2::melt(.,
                   id.var = "class")
  
  #####################################
  
  matrix_plot <- ggplot(contigency_matrix_species,
                        mapping = aes(x=variable,
                                      y=class,
                                      fill=factor(value))) +
    geom_tile() +
    scale_fill_manual(values = colorRampPalette(brewer.pal(9, "BuGn"))(32)) +
    geom_text(aes(label = value)) +
    scale_y_reverse(breaks = c(0:8)) +
    # add labels
    labs(x = paste("Reclassified quantiles after", layer_name, "removed", sep = " "),
         # y-axis
         y = "Original quantiles",
         # title
         title = paste("Classification reclassifications in", str_to_title(region), sep = " "),
         # subtitle
         subtitle = "Aggregated quantile reclassifications") +
    theme_bw() +
    plot_theme +
    # remove legend
    theme(legend.position = "none")
  
  # print(matrix_plot)
  
  #####################################
  
  ggplot2::ggsave(plot = matrix_plot, filename = file.path(figure_dir, paste(region, layer_name, "consistency_matrix.tiff", sep = "_")),
                  width = 6.5, height = 4.5, units = "in", dpi = 600, compression = "lzw")
  
  #####################################
  
  # print how long it takes to generate figures and export
  print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating and adding", layer_name, "data to dataframe and export as plot", sep = " "))
}


#####################################

# Create quantile change table
quantile_table <- lapply(1:length(list), function(i) {

  # create function that will find places where column and row equal each other
  ## these will designate when no change has occurred
  add_diag <- function(x,t){
    sum(x[row(x)+col(x)==t])
  }
  
  # designate contigency matrix for analyzing
  data <- list[[i]]
  
  #####################################
  
  # create a dataframe for each summarised contigency matrix
  quant_sum <- data.frame(
    # designate dataset name
    layer_name = data_names_list[[i]],
    # calculate when the quantile classifications did not change
    no_change = add_diag(x = data, t = seq(2, 18, 2)),
    # calculate when the quantile classification increased
    # from the original iteration and the removed dataset iteration
    ## will occur when the column number is greater than the row number
    increase = sum(sapply(2:9, function(j) sum(data[c(1:j - 1), j]))),
    # calculate when the quantile classification decreased
    # from the original iteration and the removed dataset iteration
    ## will occur when the row number is greater than the column number
    decrease = sum(sapply(2:9, function(j) sum(data[j, c(1:j - 1)]))),
    # give a placeholder value of 0 for the dataset's total
    total = 0)
  
  #####################################
  
  # verify that all no change, increased, and decreased hex grids
  # equal the total number of hex grids in study area
  quant_sum$total <- sum(quant_sum$no_change,
                         quant_sum$increase,
                         quant_sum$decrease)
  quant_sum
})

# add all the results to a final table
final_quantile_table <- bind_rows(quantile_table)
View(final_quantile_table)

#####################################
#####################################

# Export data
## Tables
base::saveRDS(object = sensitivity_quant_class, file = paste(sensitivity_dir, paste(region, "sensitivity_quantile_classifications.rds", sep = "_"), sep = "/"))
base::saveRDS(object = sensitivity_quant_change, file = paste(sensitivity_dir, paste(region, "sensitivity_quantile_classification_change.rds", sep = "_"), sep = "/"))
base::saveRDS(object = sensitivity_quant_total_iteration_change, file = paste(sensitivity_dir, paste(region, "sensitivity_quantile_classification_total_iteration_change.rds", sep = "_"), sep = "/"))
base::saveRDS(object = sensitivity_quant_minmax, file = paste(sensitivity_dir, paste(region, "sensitivity_quantile_classification_minmax.rds", sep = "_"), sep = "/"))
base::saveRDS(object = final_quantile_table, file = paste(sensitivity_dir, paste(region, "sensitivity_quantile_change_table.rds", sep = "_"), sep = "/"))

## Figures
### total change
ggplot2::ggsave(plot = tc_hist, filename = file.path(figure_dir, paste(region, "quantile_total_change_histogram.tiff", sep = "_")),
                width = 6.5, height = 4.5, units = "in", dpi = 600, compression = "lzw")
ggplot2::ggsave(plot = tc_plot, filename = file.path(figure_dir, paste(region, "quantile_total_change_plot.tiff", sep = "_")),
                width = 5, height = 8, units = "in", dpi = 600, compression = "lzw")
ggplot2::ggsave(plot = min_plot, filename = file.path(figure_dir, paste(region, "minimum_quantile_value_map.tiff", sep = "_")),
                width = 5, height = 8, units = "in", dpi = 600, compression = "lzw")

### total change
ggplot2::ggsave(plot = ic_hist, filename = file.path(figure_dir, paste(region, "quantile_iteration_change_histogram.tiff", sep = "_")),
                width = 6.5, height = 4.5, units = "in", dpi = 600, compression = "lzw")
ggplot2::ggsave(plot = ic_plot, filename = file.path(figure_dir, paste(region, "quantile_iteration_change_plot.tiff", sep = "_")),
                width = 5, height = 8, units = "in", dpi = 600, compression = "lzw")

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
print(Sys.time() - start) # print how long it takes to calculate
