####################################
### 4. Pacific Coast Port Access ###
####################################

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

### ***Note: running install below will provide message that it skipped install if the package has not changed since last install
remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"), INSTALL_opts = "--no-multiarch")

# Inspect packages and system
sessionInfo()

#####################################
#####################################

# Set directories
## Define data directory (as this is an R Project, pathnames are simplified)
### Input directories
data_dir <- "data/a_raw_data"
land_dir <- "data/a_raw_data/USGSEsriWCMC_GlobalIslandsv2_Continents"
oregon_dir <- "data/a_raw_data/or_state_boundary"

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Constraint directories
constraint_gpkg <- "data/c_submodel_data/constraints.gpkg"

#### Intermediate directories
intermediate_dir <- "data/b_intermediate_data"
oregon_gpkg <- "data/b_intermediate_data/oregon.gpkg"
pacpars_gpkg <- "data/b_intermediate_data/oregon_pacpars.gpkg"

#### PACPARS directory
dir.create(paste0(intermediate_dir, "/",
                  "pacpars"))

pacpars_dir <- "data/b_intermediate_data/pacpars"

#####################################
#####################################

# Set parameters
## setback (buffer) distance
buffer <- 3704

## designate region name
region <- "california"

## layer name
layer <- "pacpars"

## designate date
date <- format(Sys.time(), "%Y%m%d")

#####################################
#####################################

## Shoreline data (USGS Global Islands Vector dataset)
## Global Island Explorer has detailed information about the dataset and can link to the paper detailing the methods (https://rmgsc.cr.usgs.gov/gie/)
## For visual inspection, navigate here: https://rmgsc.cr.usgs.gov/gie/gie.shtml
## Previously able to download the data here: https://rmgsc.cr.usgs.gov/outgoing/ecosystems/Global/USGSEsriWCMC_GlobalIslands_v3.mpk
## Research paper: https://www.tandfonline.com/doi/full/10.1080/1755876X.2018.1529714
## Metadata: https://www.sciencebase.gov/catalog/file/get/63bdf25dd34e92aad3cda273?f=__disk__66%2F87%2F8e%2F66878e86e1aa205eb35e5a0bb641d469c67d7e2c&transform=1&allowOpen=true
### 1.) Navigate to this page: https://www.sciencebase.gov/catalog/item/63bdf25dd34e92aad3cda273
### 2.) Click the USGSEsriWCMC_GlobalIslands_v3.mpk.zip extension (***Note: there may be a newer version so name could be different -- Version 3 as of `October 28 March 2023)
### 4.) Where file is located change, unzip file
### 3.) Change .mpk (Esri mappackage) to .zip
### 4.) Unzip file
### 5.) Navigate to and copy geodatabase (as of 11 October 2021, v10 and v108 had no differences in their data)
### 6.) Paste to data dictionary

### Alternative:
#### 1.) Navigate to this page: https://www.sciencebase.gov/catalog/item/63bdf25dd34e92aad3cda273
#### 2.) Click USGSEsriWCMC_GlobalIslands_v3.mpk (size = 1.46)
#### 3.) On new page that opens, mark the box for "I'm not a robot"
#### 4.) Click Submit
#### 5.) On new page, click "Download File" (alternatively, click the temporarily generated hyperlink)
#### 6.) Where file is located change .mpk (Esri mappackage) to .zip
#### 7.) Unzip file
#### 8.) Navigate to and copy geodatabase (as of 11 October 2021, v10 and v108 had no differences in their data)
#### 9.) Paste to data dictionary

### Continental land data (uncomment below if used)
# continents <- sf::st_read(dsn = land_dir, layer = "USGSEsriWCMC_GlobalIslandsv2_Continents") %>%
#   # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
#   sf::st_transform("EPSG:26910")

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

## Oregon state boundary (source: http://navigator.state.or.us/sdl/data/shapefile/k24/or_state_boundary.zip)
### Oregon state boundary with 2-nautical mile buffer
#### ***Note: The 2-nautical mile buffer is added so points can be generated to later generate a polygon to meet PACPARS documentation
oregon_boundary <- st_read(dsn = oregon_dir, layer = "or_state_boundary") %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

#####################################

list.files(data_dir)

# PACPARS PDF (source: https://navcen.uscg.gov/sites/default/files/pdf/PARS/PAC_PARS_22/Draft%20PAC-PARS.pdf)
pdf <- paste(data_dir, "pacpars_draft_report.pdf", sep = "/")

#####################################
#####################################

# Oregon boundary
## 2-nautical mile boundary
oregon_boundary_2nm <- oregon_boundary %>%
  # filter out coastal water areas (Feature = 0)
  dplyr::filter(FEATURE != 0) %>%
  # set a buffer of 2 nautical miles to match PACPARS documentation (see page 33) -- 1 nautical mile = 1852 meters
  sf::st_buffer(dist = buffer)

## Oregon 3-nautical mile area (***Note: the data have areas between coast and 3-nautical mile boundary)
oregon_coast_area <- oregon_boundary %>%
  # get coastal nautical boundary (3 nautical miles)
  dplyr::filter(FEATURE == 0 & ORBNDY24_ == 2)

## Oregon nautical areas (linestring)
oregon_nautical_area <- oregon_coast_area %>%
  # remove area that is between coast and 2 nautical miles
  rmapshaper::ms_erase(oregon_boundary_2nm) %>%
  # convert from polygon to multistring
  sf::st_cast("MULTILINESTRING") %>%
  # convert to linestring
  sf::st_cast("LINESTRING") %>%
  # create field called "length"
  dplyr::mutate(length = sf::st_length(geometry)) %>%
  # arrange features by length from largest to smallest
  dplyr::arrange(desc(length))

### Oregon 2 - 3 nautical mile boundary
oregon_3nm <- oregon_nautical_area %>%
  # get only the longest feature (will be the north-south linestring)
  head(1)

### Oregon 2 - 3 nautical mile boundaries with points
oregon_3nm_points <- oregon_3nm %>%
  # create points along linestring
  sf::st_cast("POINT")

### Oregon eastern points (2 nautical mile linestring points)
oregon_2nm_east <- oregon_3nm_points %>%
  # create fields "fid" (to filter for eastern points) and "polygon"
  dplyr::mutate(fid = row_number(),
                polygon = "polygon") %>%
  # filter for only the eastern most points (top northeastern point = 78603)
  dplyr::filter(fid >= 78603) %>%
  # select field "polygon" (will be used to combine and generate data for PACPARS coastal area)
  dplyr::select(polygon)

#####################################
#####################################

# PACPARS
## Extract the two tables to create the data
### D13 Offshore Fairway and D13 Coastal Fairway
table <- tabulizer::extract_tables(pdf,
                                   # page tables are found (33)
                                   pages = 33)

## View D13 Offshore Fairway (page 33) -- raw
View(table[[1]])

## View D13 Coastal Fairway Zone (page 33) -- raw
View(table[[2]])

#####################################
#####################################

# Create offshore tables
## D13 Offshore Fairway table (page 33)
d13_offshore_table <- data.frame(table[[1]]) %>%
  # delete the 1st row as it does not relevant information
  dplyr::filter(!row_number() %in% c(1)) %>%
  # split out column X5 to be two columns (not calling "lat" and "lon" since other columns already have those names)
  tidyr::separate(X5, into=c("latitude", "longitude"), sep=" ", remove=T, convert = T) %>%
  # rename columns
  dplyr::rename(point = X1,
                # latitude
                lat = X2,
                # longitude
                lon = X3) %>%
  # delete the 1st row as it does not relevant information
  dplyr::filter(!row_number() %in% c(1)) %>%
  # replace any "" values with NA
  dplyr::mutate(across(where(is.character),
                       ~na_if(x = .,
                              y = ""))) %>%
  # split the data across columns by 3 (point, lat, lon)
  split.default(., (seq_along(.)-1) %/% 3) %>%
  # recombine so in long format (all features share same fields)
  dplyr::bind_rows()

### D13 offshore area (first set)
d13_offshore1 <- d13_offshore_table %>%
  # select columns of interest
  dplyr::select(point, lat, lon)

### D13 offshore area (second set)
d13_offshore2 <- d13_offshore_table %>%
  # select columns of interest
  dplyr::select(X4:longitude) %>%
  # rename columns
  dplyr::rename(point = X4,
                lat = latitude,
                lon = longitude)

### D13 offshore area (third set)
d13_offshore3 <- d13_offshore_table %>%
  # select columns of interest
  dplyr::select(X6:X8) %>%
  # rename columns
  dplyr::rename(point = X6,
                lat = X7,
                lon = X8)

### D13 offshore (combined)
d13_offshore <- dplyr::bind_rows(d13_offshore1,
                                 d13_offshore2,
                                 d13_offshore3) %>%
  # remove NA values
  drop_na(lat)

#####################################

## Create polygons
### D13 Offshore Fairway
#### Area 1
##### Point locations
d13_offshore_point1 <- d13_offshore %>%
  # filter the points to create fairway area affecting offshore call area
  dplyr::filter(point %in% c(10:11,
                             22:25)) %>%
  # convert to simple feature
  sf::st_as_sf(coords = c("lon", "lat"),
               # set the coordinate reference system to WGS84
               crs = 4326) %>% # EPSG 4326 (https://epsg.io/4326)
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

##### Polygon creation
d13_offshore_polygon <- d13_offshore_point1 %>%
  # add new field called "polygon"
  dplyr::mutate(polygon = "polygon") %>%
  # group by the points field
  dplyr::group_by(polygon) %>%
  # combine geometries without resolving borders to create multipoint feature
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  # convert back to sf
  sf::st_as_sf() %>%
  # convert to polygon simple feature
  sf::st_cast("POLYGON") %>%
  # convert back to sf
  sf::st_as_sf()

#####################################
#####################################

# Create coastal tables
## D13 Coastal Fairway table (page 33)
d13_coastal_table <- data.frame(table[[2]]) %>%
  # clean table
  janitor::row_to_names(row_number = 2,
                        # remove row after setting as column names
                        remove_row = TRUE,
                        # remove above rows
                        remove_rows_above = TRUE) %>%
  # clean column names
  janitor::clean_names()

### set cells in row 5 and 6 of column 6 as NA
d13_coastal_table[5,6] <- NA # row 5, column 6
d13_coastal_table[6,6] <- NA # row 6, column 6

### Split columns to have cleaned longitudes and latitudes
d13_coastal_table <- d13_coastal_table %>%
  # split out column to be two columns
  tidyr::separate(latitude_longitude, into=c("latitude2", "longitude2"), sep=" ", remove=T, convert = T) %>%
  # split out column to be two columns
  tidyr::separate(latitude_longitude_2, into=c("latitude3", "longitude3"), sep=" ", remove=T, convert = T) %>%
  # split the data across columns by 3 (point, lat, lon)
  split.default(., (seq_along(.)-1) %/% 3) %>%
  # recombine so in long format (all features share same fields)
  dplyr::bind_rows()

## D13 coast area (first set)
d13_coastal1 <- d13_coastal_table %>%
  # select the columns of interest
  dplyr::select(x:longitude) %>%
  # rename column names
  dplyr::rename(polygon = x,
                lat = latitude,
                lon = longitude) %>%
  # change values in "lat" and "lon" fields from character to numeric
  dplyr::mutate(across(lat:lon, as.numeric))

## D13 coast area (second set)
d13_coastal2 <- d13_coastal_table %>%
  # select the columns of interest
  dplyr::select(x_2:longitude2) %>%
  # rename columns
  dplyr::rename(polygon = x_2,
                lat = latitude2,
                lon = longitude2)

## D13 coast area (third set)
d13_coastal3 <- d13_coastal_table %>%
  # select the columns of interest
  dplyr::select(latitude3:longitude3) %>%
  # rename columns
  dplyr::rename(lat = latitude3,
                lon = longitude3) %>%
  # create field "polygon"
  dplyr::mutate(polygon = "polygon")

## D13 coast area (combined)
d13_coastal <- dplyr::bind_rows(d13_coastal1,
                                d13_coastal2,
                                d13_coastal3) %>%
  # remove NA values
  drop_na(lat) %>%
  # fill "polygon" with value polygon
  ## ***Note: values from from first two sets use numbers for points
  dplyr::mutate(polygon = "polygon")

#####################################

## Create polygons
### D13 Coastal Fairway
#### Convert table to simple feature
##### Point locations
d13_coastal_point <- d13_coastal %>%
  # convert to simple feature
  sf::st_as_sf(coords = c("lon", "lat"),
               # set the coordinate reference system to WGS84
               crs = 4326) %>% # EPSG 4326 (https://epsg.io/4326)
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910")

##### Polygon creation (D13 coast PACPARS polygon)
d13_coast_3nm <- d13_coastal_point %>%
  # combine with points marking the 2-nautical mile eastern boundary
  dplyr::bind_rows(oregon_2nm_east) %>%
  # group by the points field
  dplyr::group_by(polygon) %>%
  # combine geometries without resolving borders to create multipoint feature
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  # convert back to sf
  sf::st_as_sf() %>%
  # convert to polygon simple feature
  sf::st_cast("POLYGON") %>%
  # convert back to sf
  sf::st_as_sf() %>%
  # remove area up to 3 nautical miles
  rmapshaper::ms_erase(oregon_coast_area)

#####################################
#####################################

# PACPARS in call areas
## D13 Offshore fairway
d13_offshore_call_area <- d13_offshore_polygon %>%
  rmapshaper::ms_clip(target = .,
                      clip = oregon_call_areas)

oregon_hex_d13_offshore <- oregon_hex[d13_offshore_call_area, ] %>%
  # spatially join D13 offshore fairway values to Oregon hex cells
  sf::st_join(x = .,
              y = d13_offshore_call_area,
              join = st_intersects)

## D13 Coastal fairway
d13_coastal_call_area <- d13_coast_3nm %>%
  rmapshaper::ms_clip(target = .,
                      clip = oregon_call_areas)

oregon_hex_d13_coastal <- oregon_hex[d13_coastal_call_area, ] %>%
  # spatially join D13 coastal fairway values to Oregon hex cells
  sf::st_join(x = .,
              y = d13_coastal_call_area,
              join = st_intersects)

#####################################

# Oregon fairway areas in call area
oregon_hex_d13_fairway <- oregon_hex_d13_offshore %>%
  # combine offshore fairway area with coastal area
  dplyr::bind_rows(oregon_hex_d13_coastal) %>%
  # add field "layer" and populate with "fairway"
  dplyr::mutate(layer = "fairway") %>%
  # select fields of importance
  dplyr::select(index, layer)

#####################################
#####################################

# Export data
## Constraints geopackage
sf::st_write(oregon_hex_d13_fairway, dsn = constraint_gpkg, layer = paste0(region, "_hex_", layer), append = F)

## PACPARS
### D13 Offshore Fairway
saveRDS(object = d13_offshore_table, file = paste(pacpars_dir, "oregon_pacpars_ds13_offshore_table.rds", sep = "/"))
saveRDS(d13_offshore1, paste(pacpars_dir, "oregon_pacpars_ds13_offshore1.rds", sep = "/"))
saveRDS(d13_offshore2, paste(pacpars_dir, "oregon_pacpars_ds13_offshore2.rds", sep = "/"))
saveRDS(d13_offshore3, paste(pacpars_dir, "oregon_pacpars_ds13_offshore3.rds", sep = "/"))

saveRDS(d13_offshore, paste(pacpars_dir, "d13_offshore_points.rds", sep = "/"))
write.csv(d13_offshore, paste(intermediate_dir, "d13_offshore_points.csv", sep = "/"), row.names=FALSE)

sf::st_write(d13_offshore_point1, dsn = pacpars_gpkg, layer = "oregon_pacpars_ds13_offshore_points", append = F)
sf::st_write(d13_offshore_polygon, dsn = pacpars_gpkg, layer = "oregon_pacpars_ds13_offshore_polygon", append = F)

### D13 Coastal Fairway
saveRDS(d13_coastal_table, paste(pacpars_dir, "oregon_pacpars_ds13_coastal_table.rds", sep = "/"))
saveRDS(d13_coastal1, paste(pacpars_dir, "oregon_pacpars_d13_coastal1.rds", sep = "/"))
saveRDS(d13_coastal2, paste(pacpars_dir, "oregon_pacpars_d13_coastal2.rds", sep = "/"))
saveRDS(d13_coastal3, paste(pacpars_dir, "oregon_pacpars_d13_coastal3.rds", sep = "/"))
saveRDS(d13_coastal, paste(pacpars_dir, "oregon_pacpars_d13_coastal.rds", sep = "/"))

sf::st_write(d13_coastal_point, dsn = pacpars_gpkg, layer = "oregon_pacpars_d13_coastal_point", append = F)
sf::st_write(d13_coast_3nm, dsn = pacpars_gpkg, layer = "oregon_pacpars_d13_coast_3nm_polygon", append = F)

### D13 Fairway (combined) -- call area
sf::st_write(d13_offshore_call_area, dsn = pacpars_gpkg, layer = "oregon_pacpars_ds13_offshore_call_area", append = F)
sf::st_write(oregon_hex_d13_offshore, dsn = pacpars_gpkg, layer = "oregon_pacpars_ds13_offshore_hex", append = F)
sf::st_write(d13_coastal_call_area, dsn = pacpars_gpkg, layer = "oregon_pacpars_ds13_coastal_call_area", append = F)
sf::st_write(oregon_hex_d13_coastal, dsn = pacpars_gpkg, layer = "oregon_pacpars_ds13_coastal_hex", append = F)
sf::st_write(oregon_hex_d13_fairway, dsn = pacpars_gpkg, layer = "oregon_pacpars_ds13_fairway_hex", append = F)

## Oregon boundary
sf::st_write(oregon_boundary_2nm, dsn = oregon_gpkg, layer = "oregon_state_boundary_2nm_buffer", append = F)
sf::st_write(oregon_coast_area, dsn = oregon_gpkg, layer = "oregon_coastal_boundary", append = F)
sf::st_write(oregon_nautical_area, dsn = oregon_gpkg, layer = "oregon_nautical_boundary", append = F)
sf::st_write(oregon_3nm, dsn = oregon_gpkg, layer = "oregon_3nm", append = F)
sf::st_write(oregon_3nm_points, dsn = oregon_gpkg, layer = "oregon_3nm_points", append = F)

sf::st_write(oregon_2nm_east, dsn = oregon_gpkg, layer = "oregon_2nm_east_points", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
