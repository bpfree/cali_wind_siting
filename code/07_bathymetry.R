##############################
### 7. Bathymetry Contours ###
##############################

# Clear environment
rm(list = ls())

# Calculate start time of code (determine how long it takes to complete all code)
start <- Sys.time()

#####################################
#####################################

# Load packages
## Need to install a development version of terra to open the netCDF
### ***Note: May need restart R upon installing (stop running after first installation)
install.packages('terra', repos='https://rspatial.r-universe.dev')

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

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Intermediate directories
intermediate_dir <- "data/b_intermediate_data"

#### Bathymetry directory
dir.create(paste0(intermediate_dir, "/",
                  "bathymetry"))

bathymetry_dir <- "data/b_intermediate_data/bathymetry"
bathymetry_gpkg <- "data/b_intermediate_data/bathymetry/oregon_bathymetry.gpkg"

#####################################
#####################################

# Set parameters
## designate region name
region <- "california"

## layer names
layer <- "contours"

## contour distance
contour <- 50

## designate date
date <- format(Sys.time(), "%Y%m%d")

#####################################
#####################################

# Load data

## ***Note: Can find bathymetry data for the United States and associated areas here: https://www.ncei.noaa.gov/maps/bathymetry/
## 1.) In the right panel, mark "DEM Footprints" (can also uncheck any other layers so easier to read)
## 2.) Zoom into area(s) of interest
## 3.) Click on the map within area(s) of interest
## 4.) Within the pop-up window in the map panel, expand the NCEI Digital Elevation Models directory
## 5.) Click the source that is desired (clicking magnifying glass to right will center map view on that data source and extent)
## 6.) Click "Link to Metadata" to open new tab for data source and download data

## Central Oregon (source: https://www.ngdc.noaa.gov/thredds/fileServer/regional/central_oregon_13_navd88_2015.nc)
### Metadata: https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ngdc.mgg.dem:11500;view=iso
### ***Note: mean high water data are not available for download, but the NAVD88 are available
central_oregon <- terra::rast(paste(data_dir, "central_oregon_13_navd88_2015.nc", sep = "/"))

## Port Orford (source: https://www.ngdc.noaa.gov/thredds/fileServer/regional/port_orford_13_mhw_2008.nc)
### Metadata: https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ngdc.mgg.dem:410/html
### ***Note: no data exist for NAVD88 so the mean high water data are integrated
port_orford <- terra::rast(paste(data_dir, "port_orford_13_mhw_2008.nc", sep = "/"))

## Crescent City (source: https://www.ngdc.noaa.gov/thredds/fileServer/regional/crescent_city_13_navd88_2010.nc)
### Metadata: https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ngdc.mgg.dem:693/html
### ***Note: NAVD88 (https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ngdc.mgg.dem:724/html)
###          and mean high water data both exist; the mean high water data are downloaded in case the MHW data
###          for Central Oregon become available
crescent_city <- terra::rast(paste(data_dir, "crescent_city_13_navd88_2010.nc", sep = "/"))

## Coastal Relief Model volume 7 (source: https://www.ngdc.noaa.gov/thredds/fileServer/crm/crm_vol7.nc))
### Metadata: https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ngdc.mgg.dem:348
### ***Note: the northern portion of the Brookings call area fall outside of smaller, higher resolution areas.
###          These data cover that area, but also cover the smaller, higher resolution areas.
crm_v7 <- terra::rast(paste(data_dir, "crm_vol7.nc", sep = "/"))

#####################################

## Oregon call areas
oregon_call_areas <- sf::st_read(dsn = wind_area_gpkg,
                                 layer = paste(sf::st_layers(dsn = wind_area_gpkg,
                                                             do_count = TRUE))) %>%
  # add a 20-meter buffer to call area to help contours extend beyond real call areas
  sf::st_buffer(dist = 20) %>%
  # reproject data into the coordinate system that match the DEM data (some of the DEMs expand beyond UTM 10N)
  sf::st_transform("EPSG:4326")

#####################################
#####################################

# Inspect the data
## Coordinate reference systems
terra::crs(central_oregon) # EPSG:4326
terra::crs(port_orford) # EPSG:4326
terra::crs(crescent_city) # EPSG:4326
terra::crs(crm_v7) # EPSG:9122, yet the metadata page states that the horizontal datum is NAD83 (EPSG: 4269)
cat(crs(crm_v7))

## Resolution
terra::res(central_oregon) # 9.259259e-05 x 9.259259e-05
terra::res(port_orford) # 9.259259e-05 x 9.259259e-05
terra::res(crescent_city) # 9.259259e-05 x 9.259259e-05
terra::res(crm_v7) # 0.0008333333 x 0.0008333333 --> so this needs to be downscaled to match other datasets

#####################################
#####################################

# Prepare the coastal relief model volume 7 data
## Change the coordinate reference system to match other datasets
terra::crs(crm_v7) <- "EPSG:4326"

## Disaggregate the resolution to match the other datasets
### Calculate the factor between the resolutions (9)
factor <- terra::res(crm_v7)[1] / terra::res(crescent_city)[1] 

### Disaggreate the resolution
crm_v7_disagg <- terra::disagg(x = crm_v7,
                               # resolution factor
                               fact = factor)

### Set other aspects of the data
#### Units
units(crm_v7_disagg) <- "meters"

#### Variable names
varnames(crm_v7_disagg) <- "z"

### Reinspect data
cat(crs(crm_v7_disagg))

#####################################

# Crop coastal relief model data to Oregon call areas
## ***Note: the coastal relief model is large, thus to make
##          the data more manageable, it can be cropped and
##           masked to the area of interest (Oregon call areas)
crm_v7_call_area <- crm_v7_disagg %>%
  # crop the CRM volume 7 data to the Oregon call areas
  terra::crop(oregon_call_areas,
              # use the Oregon call areas to mask the data
              mask = T)

#####################################
#####################################

# Combine the bathymetry datasets together
## ***Note: due to the high resolutions and large areas,
##          to minimize crashing R, the datasets are combined piecemeal.
## ***Warning: This will affect the mean calculations for new datasets
##             get added to the previous mean before a new mean is calculated.

### Oregon bathymetry using 1/3-arc second datasets
oregon_bath <- terra::mosaic(central_oregon,
                             port_orford,
                             crescent_city,
                             fun = "mean")

### Oregon call area bathymetry with the addition of the CRM volume 7 dataset
oregon_bath <- terra::mosaic(oregon_bath,
                             crm_v7_call_area,
                             fun = "mean") %>%
  # crop the CRM volume 7 data to the Oregon call areas
  terra::crop(oregon_call_areas,
              # use the Oregon call areas to mask the data
              mask = T)

### ***Note: if computer can handle it, this combination will create the most accurate dataset
# oregon_bath <- terra::mosaic(central_oregon,
#                              port_orford,
#                              crescent_city,
#                              crm_v7_call_area,
#                              fun = "mean")

#####################################
#####################################

# Create bathymetry contours
oregon_contours <- terra::as.contour(x = oregon_bath,
                                     maxcells = 100500000,
                                     # have contours match desired distance (likely to be 50 meters) 
                                     ## minimum is multiplied by -1 so number of levels can be calculated
                                     nlevels = ((-1 * terra::minmax(oregon_bath)[1] + terra::minmax(oregon_bath)[2]) / contour)) %>%
  # set as a simple feature
  sf::st_as_sf() %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910") %>% # EPSG 26910 (https://epsg.io/26910) %>%
  # calculate length of objects to later help with filtering
  dplyr::mutate(length = as.numeric(sf::st_length(.)))

#####################################
#####################################

# Export data
## Bathymetry data
sf::st_write(obj = oregon_contours, dsn = bathymetry_gpkg, layer = paste0(region, "_", layer, "_", contour, "m"), append = F)

## Intermediate data
### ***WARNING: These files are large (e.g., oregon_bath = 2.3 GB, crm_v7_disagg = 33.5 GB)
### ***Suggestion: Export if want to inspect
# terra::writeRaster(oregon_bath, filename = file.path(bathymetry_dir, "oregon_study_area_bath.grd"), overwrite = T)
# terra::writeRaster(crm_v7_disagg, filename = file.path(bathymetry_dir, "noaa_crm_v7_disagg_bath.grd"), overwrite = T)
# terra::writeRaster(crm_v7_call_area, filename = file.path(bathymetry_dir, "oregon_call_area_crm_v7.grd"), overwrite = T)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
