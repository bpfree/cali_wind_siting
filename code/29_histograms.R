################################################
### 29. Sensitivity Analysis -- Histograms #####
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

sf::st_layers(dsn = sensitivity_gpkg,
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
sensitivity_nominal_percent_change <- sf::st_read(dsn = sensitivity_gpkg,
                                                  layer = paste(region, "sensitivity_nominal_percent_change", sep = "_"))

table <- base::readRDS(file = paste("data/f_sensitivity_data", "oregon_sensitivity_table.rds", sep = "/"))

#####################################
#####################################

# histogram of original analysis
h <- ggplot() +
  geom_histogram(data = sensitivity_nominal_percent_change, aes(x = sensitivity_nominal_percent_change[["model_geom_mean"]]),
                 fill = "#869BC4", color = "#869BC4", alpha = 0.8) +
  labs(x = "Relative score",
       y = "Occurances",
       title = paste("Relative suitability in", str_to_title(region), sep = " "), 
       subtitle = "Base analysis",
       caption = "All datasets included in model run") +
  theme_bw() +
  plot_theme

print(h)

#####################################
#####################################

# # Create a list for the plots to be placed
# plot_list <- list()

# Produce histograms of the minimum and maximum changes 
for (i in 0:23){
  start2 <- Sys.time()
  
  # if wanting to test a particular dataset
  #i <- 0
  
  dataset_list <- list(table$dataset_list)
  
  # Nominal histograms
  ## designate the field of interest
  nominal_field <- 3 + i*2
  
  ## obtain dataset name
  dataset_name <- dataset_list[[1]][i+1]
  
  ## create histogram
  h_nominal <- ggplot() +
    geom_histogram(data = sensitivity_nominal_percent_change, aes(x = sensitivity_nominal_percent_change[[nominal_field]]),
                   fill = "#869BC4", color = "#869BC4", alpha = 0.8) +
    labs(x = "Score change (nominal)",
         y = "Occurances",
         title = paste("Sensitivity Analysis in", str_to_title(region), sep = " "), 
         subtitle = "Nominal change in final suitability score with dataset removal from original score",
         caption = paste("Dataset Removed:", dataset_name, sep = " ")) +
    theme_bw() +
    plot_theme
  
  print(h_nominal)
  
  ggsave(h_nominal, filename = file.path(figure_dir, paste(i+1, region, str_to_lower(dataset_name), "nominal_change_plot.tiff", sep = "_")),
         width=6.5, height = 4.5, units = "in", dpi = 600, compression = "lzw")
  
  #####################################
  
  # Percent histograms
  percent_field <- 4 + i*2
  
  h_percent <- ggplot() +
    geom_histogram(data = sensitivity_nominal_percent_change, aes(x = sensitivity_nominal_percent_change[[percent_field]]),
                   fill = "#869BC4", color = "#869BC4", alpha = 0.8) +
    labs(x = "Score change (%)",
         y = "Occurances",
         title = paste("Sensitivity Analysis in", str_to_title(region), sep = " "), 
         subtitle = "Percent change in final suitability score with dataset removal from original score",
         caption = paste("Dataset Removed:", dataset_name, sep = " ")) +
    theme_bw() +
    plot_theme
  
  print(h_percent)
  
  ggsave(h_percent, filename = file.path(figure_dir, paste(i+1, region, str_to_lower(dataset_name), "percent_change_plot.tiff", sep = "_")),
         width = 6.5, height = 4.5, units = "in", dpi = 600, compression = "lzw")
  
  #####################################
  
  # Add plots to the list
  #plot_list[[i + 1]] <- h_nominal + h_percent
  
  # print how long it takes to calculate
  print(paste("Iteration", i, "takes", Sys.time() - start2, "minutes to complete creating and adding", dataset_name, "histogram to plot", sep = " "))
}

#####################################

# Create a list of ggplot objects for each page
# pages <- lapply(seq(1, length(plot_list), 2), function(start_idx) {
#   end_idx <- min(length(plot_list), start_idx + 1)
#   do.call(patchwork::wrap_plots, plot_list[start_idx:end_idx])
# })
# 
# # Create the final document with cowplot
# doc <- plot_grid(plotlist = pages, ncol = 6)

#####################################
#####################################

# Export figures
ggsave(plot = h, filename = file.path(figure_dir, paste(region, "relative_suitability.tiff", sep = "_")),
       width = 6.5, height = 4.5, units = "in", dpi = 600, compression = "lzw")

# Save the document
# ggplot2::ggsave(plot = doc, filename = file.path(figure_dir, paste(region, "nominal_percent_changes_combined.pdf", sep = "_")),
#                 width = 20, height = 24)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
