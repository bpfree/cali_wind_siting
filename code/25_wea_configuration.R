##########################################
### 25. Wind Energy Area Configuration ###
##########################################


# Clear environment
rm(list = ls())

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
lisa_gpkg <- "data/e_rank_data/lisa.gpkg"
aliquot_dir <- "data/a_raw_data/PC_ALIQUOTS"

## Output directories
### configuration directory
wea_config_gpkg <- "data/e_rank_data/wea_configuration.gpkg"

#####################################
#####################################

# Load data
## LISA high-high clusters
oregon_lisa_highhigh <- sf::st_read(dsn = lisa_gpkg, layer = "oregon_hex_lisa_highhigh")

## Pacific aliquot
pacific_aliquot <- sf::st_read(dsn = aliquot_dir, layer = "Aliquots_WestCoast_MMC")

#####################################
#####################################


