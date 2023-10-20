#################################
### 10. Rocky Reef Groundfish ###
#################################

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
rocky_reef_mapped_dir <- "data/a_raw_data/OSU_ATSML_SGH_V4_0_WA_OR_NCA"
rocky_reef_probable_gpkg <- "data/a_raw_data/rocky_reef_probable.gpkg"
study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Submodel directory
natural_resources_submodel <- "data/c_submodel_data/oregon_natural_resources_submodel.gpkg"

#### Intermediate directories
rocky_reef_gpkg <- "data/b_intermediate_data/oregon_rocky_reef.gpkg"

#####################################
#####################################

# Set parameters
## setback (buffer) distance
buffer <- 500

## designate region name
region <- "california"

## layer names
layer_mapped <- "rocky_reef_mapped"
layer_probable <- "rocky_reef_probable"

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

## ***Note: These data were provided by Curt Whitmire (curt.whitmire@noaa.gov) [NWFSC] for the data are not publicly available for download
### Original report (Goldfinger et al. 2014): https://espis.boem.gov/final%20reports/5453.pdf
### Metadata: https://www.webapps.nwfsc.noaa.gov/server7/rest/services/FRAM/USWestCoast_SeafloorInduration_v4_0/MapServer
### Can view the data on the West Coast Data Portal: https://portal.westcoastoceans.org/visualize/#x=-124.32&y=43.47&z=7&logo=true&controls=true&dls%5B%5D=true&dls%5B%5D=0.5&dls%5B%5D=778&basemap=ocean&themes%5Bids%5D%5B%5D=4&tab=active&legends=false&layers=true
### Also view data on FRAM (under habitat area): https://www.webapps.nwfsc.noaa.gov/data/map
rocky_reef_mapped <- sf::st_read(dsn = rocky_reef_mapped_dir, layer = "V4_0_SGH_WA_OR_NCA") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

## ***Note: Curt Whitmire (curt.whitmire@noaa.gov) generated the probable rocky reef habitat for groundfish using
##          two datasets:
##              (1) multibeam acoustic backscatter imagery from the NOAA Ship Fairweather (research cruise WOO474
##                  https://www.ngdc.noaa.gov/nos/W00001-W02000/W00474.html), and
##              (2) unpublished SeaBED-class AUV image annotations from two dives conducted during an 2019 EXPRESS
##                  cruise aboard the NOAA Ship Lasker in 2019 (RL-19-05).

### These data were manually created
rocky_reef_probable <- sf::st_read(dsn = rocky_reef_probable_gpkg, layer = "rockyreef_hapc_probable") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

#####################################
#####################################

# Hard and mixed rocky reef habitats
rocky_reef_mapped_hard_mixed <- rocky_reef_mapped %>%
  sf::st_make_valid() %>%
  # subset to hard and mixed induration (IND)
  dplyr::filter(IND %in% c("hard", "mixed")) %>%
  # add field "layer" and populate with "rocky reef (mapped)"
  dplyr::mutate(layer = "rocky reef (mapped)") %>%
  # group by the "layer" field
  dplyr::group_by(layer) %>%
  # summarise all features by "layer" field to get single object
  dplyr::summarise()

#####################################

## ***Note: The analysis calls for a 500-meter buffer around the mapped rocky reef habitats.
##          Normally, the buffer distance would get added to the areas and then clipped to the call areas.
##          However, the processing for that took too long.
##          Instead a buffer of a similar distance was added to the call area to get all possible habitat areas,
##          then those areas had a 500-meter buffer added before focused on areas that would venture into the
##          call areas.

## ***Note: this is important as habitats fall just beyond the Oregon call areas but when distances are
##          applied, they then cross into the call areas.

call_area_500m <- oregon_call_areas %>%
  # add the 500-meter buffer to the call areas
  sf::st_buffer(dist = buffer)

rocky_reef_oregon_hard_mixed <- rocky_reef_mapped_hard_mixed %>%
  # extract habitat areas that fall within the extend call areas
  rmapshaper::ms_clip(target = .,
                      clip = call_area_500m) %>%
  # add a 500-meter buffer to those areas
  sf::st_buffer(dist = buffer) %>%
  # limit back to the original Oregon call areas
  rmapshaper::ms_clip(target = .,
                      clip = oregon_call_areas)

oregon_hex_rocky_reef_mapped <- oregon_hex[rocky_reef_oregon_hard_mixed, ] %>%
  # spatially join NMFS EFHCA values to Oregon hex cells
  sf::st_join(x = .,
              y = rocky_reef_oregon_hard_mixed,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(index, layer)

#####################################

# Probable rocky reef habitat
rocky_reef_probable500 <- rocky_reef_probable %>%
  # apply 500-meter buffer onto probable area
  sf::st_buffer(dist = buffer)

oregon_hex_rocky_reef_probable <- oregon_hex[rocky_reef_probable500, ] %>%
  # spatially join NMFS EFHCA values to Oregon hex cells
  sf::st_join(x = .,
              y = rocky_reef_probable500,
              join = st_intersects) %>%
  # add field "layer" and populate with "rocky reef (probable)"
  dplyr::mutate(layer = "rocky reef (probable)") %>%
  # select fields of importance
  dplyr::select(index, layer)

#####################################
#####################################

# Export data
## Analysis geopackage
sf::st_write(oregon_hex_rocky_reef_mapped, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", layer_mapped, "_", buffer, "m"), append = F)
sf::st_write(oregon_hex_rocky_reef_probable, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", layer_probable, "_", buffer, "m"), append = F)

## Rocky Reef geopackage
sf::st_write(rocky_reef_mapped, dsn = rocky_reef_gpkg, layer = "sgh_v4_rocky_reef_mapped", append = F)
sf::st_write(rocky_reef_mapped_hard_mixed, dsn = rocky_reef_gpkg, layer = "sgh_v4_rocky_reef_mapped_hard_mixed", append = F)  
sf::st_write(rocky_reef_oregon_hard_mixed, dsn = rocky_reef_gpkg, layer = "oregon_sgh_v4_rocky_reef_mapped_hard_mixed", append = F)
sf::st_write(oregon_hex_rocky_reef_mapped, dsn = rocky_reef_gpkg, layer = "oregon_hex_rocky_reef_mapped", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
