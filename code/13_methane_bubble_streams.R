##################################
### 13. Methane Bubble Streams ###
##################################

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
## Define data directory (as this is an R Project, pathnames are simplified)
### Input directories
data_dir <- "data/a_raw_data"
methane_gdb <- "data/a_raw_data/NCCOS_OR_Modeling-Habitat_20220825/v107/sh_22_09_express_sep2022.gdb"

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Submodel directory
natural_resources_submodel <- "data/c_submodel_data/oregon_natural_resources_submodel.gpkg"

#### Intermediate directories
methane_bubble_gpkg <- "data/b_intermediate_data/oregon_methane_bubble_streams.gpkg"

#####################################

# Inspect available layers and names
sf::st_layers(dsn = methane_gdb,
              do_count = T)

#####################################
#####################################

# Set parameters
## setback (buffer) distance
buffer <- 1000

## designate region name
region <- "california"

## layer names
layer <- "methane_bubble_streams"

## designate date
date <- format(Sys.time(), "%Y%m%d")

#####################################
#####################################

# Load data
## Oregon call areas
oregon_call_areas <- sf::st_read(dsn = wind_area_gpkg,
                                 layer = paste(sf::st_layers(dsn = wind_area_gpkg,
                                                             do_count = TRUE)))

## Oregon hex areas (original data)
oregon_hex <- sf::st_read(dsn = study_area_gpkg,
                          layer = paste(sf::st_layers(dsn = study_area_gpkg,
                                                      do_count = TRUE)[[1]][4]))

#####################################

# Methane bubble streams
## ***Note: These data were compiled by Curt Whitmire (curt.whitmire@noaa.gov)
### Two major datasets contributed to the methane bubble streams:
#### 1.) Merle et al. (2021): https://www.pmel.noaa.gov/eoi/Cascadia/Supplemental-Tables-US-only-revised-dec30-2020.xlsx
####     These data came from the supplemental table in the paper: https://www.frontiersin.org/articles/10.3389/feart.2021.531714/full
####     They also contain data from two other datasets: Reidel et al. (2018) and Johnson et al. (2015)
####     i.) Reidel et al. (2018): https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-05736-x/MediaObjects/41467_2018_5736_MOESM4_ESM.xlsx
#####     ***Note: data come from the supplementary data 2
#####     Paper: https://www.nature.com/articles/s41467-018-05736-x
####     ii.) Johnson et al. (2015): https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2015GC005955&file=ggge20859-sup-0001-2015GC005955-SupInfo.docx
#####     ***Note: data come from the supporting information document (see S2 and S3)
#####     ***Note: S3 does not contain any sites that fall within Oregon call areas
#####     Paper: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GC005955

#### 2.) Rainier (H13118 (2018)) data ahead of a 2019 EXPRESS cruise aboard the NOAA Ship Lasker (RL-19-05)
####     These are unpublished location data and came from Nancy Prouty at USGS (nprouty@usgs.gov)

methane_bubble_streams <- sf::st_read(dsn = methane_gdb, layer = "USGS_seeps_merge") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910") %>%
  # drop M geometry (methane data has XYZM geometries)
  sf::st_zm()

#####################################
#####################################

# Buffered (1km) methane bubble streams within Oregon call areas
oregon_methane_bubbles <- methane_bubble_streams %>%
  # apply a 1-kilometer (1000m) buffer on all methane bubble streams
  sf::st_buffer(dist = buffer) %>%
  # limit to only Oregon call areas
  rmapshaper::ms_clip(target = .,
                      clip = oregon_call_areas) %>%
  # create field called "layer" and fill with "methane bubble streams" for summary
  dplyr::mutate(layer = "methane bubble streams") %>%
  # group by layer and summarise to have single feature
  dplyr::group_by(layer) %>%
  dplyr::summarise()

#####################################
#####################################

# Methane bubble stream hex grid
oregon_hex_methane_bubbles <- oregon_hex[oregon_methane_bubbles, ] %>%
  # spatially join continental shelf values to Oregon hex cells
  sf::st_join(x = .,
              y = oregon_methane_bubbles,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(index, layer)

#####################################
#####################################

# Subdatasets within the methane bubble streams dataset (***WARNING: these exclude the Rainier subdataset)

## Methane bubble streams (Merle et al. 2021) (source: https://www.pmel.noaa.gov/eoi/Cascadia/Supplemental-Tables-US-only-revised-dec30-2020.xlsx)
### ***Note: data come from the supplemental table
### Paper: https://www.frontiersin.org/articles/10.3389/feart.2021.531714/full
#### These data also contain data from earlier research surveys
#### They include:
#####   1.) Reidel et al. 2018: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-018-05736-x/MediaObjects/41467_2018_5736_MOESM4_ESM.xlsx
#####       ***Note: data come from the supplementary data 2
#####       Paper: https://www.nature.com/articles/s41467-018-05736-x
#####   2.) Johnson et al. 2015: https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2F2015GC005955&file=ggge20859-sup-0001-2015GC005955-SupInfo.docx
#####       ***Note: data come from the supporting information document (see S2 and S3)
#####       ***Note: S3 does not contain any sites that fall within Oregon call areas
#####       Paper: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015GC005955
methane_merle_all <- readxl::read_xlsx(path = paste(data_dir, "methane_bubble_streams_merle.xlsx", sep = "/"),
                                       # designate the sheet with the data
                                       sheet = 2) %>%
  # convert to simple feature
  sf::st_as_sf(coords = c("Longitude", "Latitude"),
               # set the coordinate reference system to WGS84
               crs = 4326) %>% # EPSG 4326 (https://epsg.io/4326)
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

#####################################

## Datasets included in the Merle et al. (2021)
## ***Note: These data are separately downloaded during script 1 during the data acquisition phase.
##          Thus, you will find the raw data in the a_raw_data subdirectory within the data directory

## Methane bubble streams (Reidel et al. 2018)
methane_reidel <- readxl::read_xlsx(path = paste(data_dir, "methane_bubble_streams_reidel.xlsx", sep = "/"),
                                   # designate the sheet with the data
                                   sheet = 1)

# set column names to be values of first row
colnames(methane_reidel) <- methane_reidel[1, ]

methane_reidel <- methane_reidel %>%
  # delete the first two rows as they do not contain data
  dplyr::filter(!row_number() %in% c(1:2)) %>%
  # delete last 12 rows as they do not contain data
  # This works by saying remove any that have a "," in the Longitude field (those are citations)
  dplyr::filter(!grepl(",", Longitude)) %>%
  # remove any sites that do not contain longitude data
  dplyr::filter(!is.na(.$Longitude)) %>%
  # convert to simple feature
  sf::st_as_sf(coords = c("Longitude", "Latitude"),
               # set the coordinate reference system to WGS84
               crs = 4326) %>% # EPSG 4326 (https://epsg.io/4326)
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

#####################################

## Johnson et al. (2015)
methane_johnson_doc <- docxtractr::read_docx(file.path(data_dir, "methane_bubble_streams_johnson.docx"))

### extract data from supporting information table 2 (S2)
methane_johnson_s2 <- docxtractr::docx_extract_tbl(methane_johnson_doc,
                                     tbl_number = 1)

# set column names to be values of first row
colnames(methane_johnson_s2) <- methane_johnson_s2[1, ]

methane_johnson_s2 <- methane_johnson_s2 %>%
  # delete the first two rows as they do not contain data
  dplyr::filter(!row_number() %in% c(1:2)) %>%
  # convert to simple feature
  sf::st_as_sf(coords = c("Lon", "Lat"),
               # set the coordinate reference system to WGS84
               crs = 4326) %>% # EPSG 4326 (https://epsg.io/4326)
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

#####################################
#####################################

# Export data
## Natural Resources submodel
sf::st_write(obj = oregon_hex_methane_bubbles, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", layer, "_1km"), append = F)

## Methane geopackage
sf::st_write(obj = methane_bubble_streams, dsn = methane_bubble_gpkg, layer = "methane_bubble_streams", append = F)
sf::st_write(obj = oregon_methane_bubbles, dsn = methane_bubble_gpkg, layer = "oregon_methane_bubble_streams", append = F)
sf::st_write(obj = oregon_hex_methane_bubbles, dsn = methane_bubble_gpkg, layer = "oregon_hex_methane_bubble_streams_1km", append = F)

sf::st_write(obj = methane_merle_all, dsn = methane_bubble_gpkg, layer = "methane_bubble_stream_merle_all", append = F)
sf::st_write(obj = methane_reidel, dsn = methane_bubble_gpkg, layer = "methane_bubble_stream_reidel", append = F)
sf::st_write(methane_johnson_s2, dsn = methane_bubble_gpkg, layer = "methane_bubble_stream_johnson_s2", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
