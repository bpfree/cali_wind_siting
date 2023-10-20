############################################
### 15. Marine Bird Species Density Maps ###
############################################

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
marine_bird_dir <- "data/a_raw_data/marine_bird/0242882/1.1/data/0-data/model_output_predictions/"
marine_bird_species_dir <- "data/b_intermediate_data/marine_bird_species"

study_area_gpkg <- "data/b_intermediate_data/oregon_study_area.gpkg"
wind_area_gpkg <- "data/b_intermediate_data/oregon_wind_area.gpkg"

### Output directories
#### Submodel directory
natural_resources_submodel <- "data/c_submodel_data/oregon_natural_resources_submodel.gpkg"

#### Intermediate directories
##### Marine bird species geopackage (vector data)
marine_bird_gpkg <- "data/b_intermediate_data/marine_bird_species.gpkg"

##### Marine bird species subdirectory (raster data)
intermediate_dir <- "data/b_intermediate_data"
dir.create(paste0(intermediate_dir, "/",
                  "marine_bird_species"))

marine_bird_raster <- "data/b_intermediate_data/marine_bird_species"

#####################################
#####################################

# Data and file preparation
## see all files
list.files(marine_bird_dir)

# Set parameters
## designate region name
region <- "california"

## layer names
layer <- "marine_bird"

## designate date
date <- format(Sys.time(), "%Y%m%d")

#####################################
#####################################

# Functions
## Prepping annual datasets
annual_species_density_function <- function(directory, spp_code){
  # list all files that relate to marine bird species density
  species_season_files <- list.files(directory,
                                     # pattern for marine bird species density is "_density.tif"
                                     pattern = "_density.tif")
  
  # create species list of unique bird species
  species_list <- unique(sapply(strsplit(x = species_season_files,
                                         # split file names into elements by "_"
                                         split = "_"),
                                # function across all files is to return first element from the string
                                function(x) x[1]))
  
  # vector of species codes
  species_code <- spp_code
  
  # get seasons for species that are contained within the vector of species codes
  target_species_seasons <- list.files(directory, pattern = paste0(spp_code, "_"))
  
  # get population densities for all seasons available for target species
  annual_species_list <- target_species_seasons[sapply(strsplit(x = target_species_seasons,
                                                                # split file names into elements by "_"
                                                                split = "_"),
                                                       # function across all files is to return fourth element from the string
                                                       # when element is equal to "density.tif"
                                                       function(x) x[4] == "density.tif")]
  
  #####################################
  
  # create a summed raster across the seasons
  species_annual_raster <- sum(terra::rast(file.path(directory, annual_species_list)),
                               na.rm = T)
  
  #####################################
  
  # calculate annual total density
  species_annual_total <- terra::global(x = species_annual_raster, fun = "sum", na.rm = T)
  ## give total density to a new object for normalizing
  species_annual_sum <- species_annual_total$sum
  
  # Normalize densities
  ## Divide the annual density data by the total summed data
  species_annual_norm <- species_annual_raster[] / species_annual_sum
  ## Set the normalized data back to the original annual density dataset
  species_normalize <- terra::setValues(species_annual_raster, species_annual_norm)
}

## Cormorant species preparation
cormorant_species_function <- function(directory, spp_code){
  species_season_files <- list.files(directory,
                                     # pattern for marine bird species density is "_density.tif"
                                     pattern = "_density.tif")
  
  # create species list of unique bird species
  species_list <- unique(sapply(strsplit(x = species_season_files,
                                         # split file names into elements by "_"
                                         split = "_"),
                                # function across all files is to return first element from the string
                                function(x) x[1]))
  
  # vector of species codes
  species_code <- spp_code
  
  # get seasons for species that are contained within the vector of species codes
  target_species_seasons <- list.files(directory, pattern = paste0(spp_code, "_"))
  
  # get population densities for all seasons available for target species
  annual_species_list <- target_species_seasons[sapply(strsplit(x = target_species_seasons,
                                                                # split file names into elements by "_"
                                                                split = "_"),
                                                       # function across all files is to return fourth element from the string
                                                       # when element is equal to "density.tif"
                                                       function(x) x[4] == "density.tif")]

  #####################################
  
  # create a summed raster across the seasons
  species_annual_raster <- sum(terra::rast(file.path(directory, annual_species_list)),
                               na.rm = T)
}

## Create z-shape membership function
### Adapted from https://www.mathworks.com/help/fuzzy/zmf.html
zmf_function <- function(raster){
  # calculate minimum value
  min <- terra::minmax(raster)[1,]
  
  # calculate maximum value
  max <- terra::minmax(raster)[2,]
  
  # calculate z-score minimum value
  ## this ensures that no value gets a value of 0
  z_max <- max + (max * 1 / 1000)
  
  # calculate z-scores (more desired values get score of 1 while less desired will decrease till 0)
  z_value <- ifelse(raster[] == min, 1, # if value is equal to minimum, score as 1
                    # if value is larger than minimum but lower than mid-value, calculate based on reduction equation
                    ifelse(raster[] > min & raster[] < (min + z_max) / 2, 1 - 2 * ((raster[] - min) / (z_max - min)) ** 2,
                           # if value is larger than mid-value but lower than maximum, calculate based on equation
                           ifelse(raster[] >= (min + z_max) / 2 & raster[] < z_max, 2*((raster[] - z_max) / (z_max - min)) ** 2,
                                  # if value is equal to maximum, score min - (min * 1 / 1000); otherwise give NA
                                  ifelse(raster[] == z_max, 0, NA))))
  
  # set values back to the original raster
  zvalues <- terra::setValues(raster, z_value)
  
  # return the raster
  return(zvalues)
}

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

## Marine bird vulnerability indices
marine_bird_vulnerability <- base::readRDS(file = paste(marine_bird_species_dir, "species_vulnerabilities_normalized.rds", sep = "/"))

### Species files
#### All species and seasons
files_list <- list.files(marine_bird_dir, pattern = "_density.tif")

#### Only species names
species_list <- unique(sapply(strsplit(files_list, "_"), function(x) x[1]))

#### Only modeled species (so species not later in taxonomic group)
modeled_species <- species_list[!(species_list %in% c("POJA", "PAJA-LTJA", "RTLO", "COLO", "BRAC", "PECO", "DCCO"))]

weights <- marine_bird_vulnerability %>%
  # get species codes for modeled species
  dplyr::filter(species_code %in% modeled_species) %>%
  # select only species of interest
  dplyr::select(species_code,
                overall_vul) %>%
  # arrange alphabetically by species code
  dplyr::arrange(species_code)

# convert weights to list (only weights and not species_codes)
species_weights <- weights[[2]]

#####################################

## Marine bird species (source: https://www.ncei.noaa.gov/archive/archive-management-system/OAS/bin/prd/jquery/download/242882.1.1.tar.gz)
### ***NOTE: NOAA conducted an analysis on marine bird species that combined data from Leirness et al. (2021)
###          https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc:0242882. For more specific downloads
###          on just the model output predictors: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/

### ***NOTE: NOAA analyzed and quantified the data for 30 species and 12 taxonomic groups
###          Densities were summarized by season  (fall, spring, summer, winter) -- though not all seasons exist

#### Species
####   1.) South polar skua (Stercorarius maccormicki) -- SPSK:
####      Fall (https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/SPSK_fall_predicted_density.tif)
spsk <- annual_species_density_function(directory = marine_bird_dir, spp_code = "SPSK")

####   2.) Common murre (Uria aalge) -- COMU: 
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/COMU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/COMU_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/COMU_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/COMU_winter_predicted_density.tif
comu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "COMU")

####   3.) Pigeon guillemot (Cepphus columba) -- PIGU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PIGU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PIGU_summer_predicted_density.tif
pigu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "PIGU")

####   4.) Marbled murrelet (Brachyramphus marmoratus) -- MAMU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/MAMU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/MAMU_summer_predicted_density.tif
mamu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "MAMU")

####   5.) Ancient murrelet (Synthliboramphus antiquus) -- ANMU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/ANMU_spring_predicted_density.tif
anmu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "ANMU")

####   6.) Cassin's auklet (Ptychoramphus aleuticus) -- CAAU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CAAU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CAAU_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CAAU_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CAAU_winter_predicted_density.tif
caau <- annual_species_density_function(directory = marine_bird_dir, spp_code = "CAAU")

####   7.) Rhinoceros auklet (Cerorhinca monocerata) -- RHAU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/RHAU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/RHAU_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/RHAU_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/RHAU_winter_predicted_density.tif
rhau <- annual_species_density_function(directory = marine_bird_dir, spp_code = "RHAU")

####   8.) Tufted puffin (Fratercula cirrhata) -- TUPU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/TUPU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/TUPU_summer_predicted_density.tif
tupu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "TUPU")

####   9.) Black-legged kittiwake (Rissa tridactyla) -- BLKI:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BLKI_spring_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BLKI_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BLKI_winter_predicted_density.tif
blki <- annual_species_density_function(directory = marine_bird_dir, spp_code = "BLKI")

####   10.) Sabine's gull (Xema sabini) -- SAGU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/SAGU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/SAGU_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/SAGU_fall_predicted_density.tif
sagu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "SAGU")

####   11.) Bonaparte's gull (Chroicocephalus philadelphia) -- BOGU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BOGU_spring_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BOGU_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BOGU_winter_predicted_density.tif
bogu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "BOGU")

####   12.) Heermann's gull (Larus heermanni) -- HEEG:
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/HEEG_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/HEEG_fall_predicted_density.tif  
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/HEEG_winter_predicted_density.tif
heeg <- annual_species_density_function(directory = marine_bird_dir, spp_code = "HEEG")

####   13.) California gull (Larus californicus) -- CAGU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CAGU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CAGU_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CAGU_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CAGU_winter_predicted_density.tif
cagu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "CAGU")

####   14.) Caspian tern (Hydroprogne caspia) -- CATE:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CATE_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CATE_summer_predicted_density.tif
cate <- annual_species_density_function(directory = marine_bird_dir, spp_code = "CATE")

####   15.) Laysan albratross (Phoebastria immutabilis) -- LAAL:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/LAAL_spring_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/LAAL_winter_predicted_density.tif
laal <- annual_species_density_function(directory = marine_bird_dir, spp_code = "LAAL")

####   16.) Black-footed alaatross (Phoebastria nigripes) -- BFAL:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BFAL_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BFAL_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BFAL_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BFAL_winter_predicted_density.tif
bfal <- annual_species_density_function(directory = marine_bird_dir, spp_code = "BFAL")

####   17.) Fork-tailed storm-petrel (Hydrobates furcatus) -- FTSP:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/FTSP_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/FTSP_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/FTSP_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/FTSP_winter_predicted_density.tif
ftsp <- annual_species_density_function(directory = marine_bird_dir, spp_code = "FTSP")

####   18.) Leach's storm-petrel (Hydrobates leucorhous) -- LESP:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/LESP_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/LESP_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/LESP_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/LESP_winter_predicted_density.tif
lesp <- annual_species_density_function(directory = marine_bird_dir, spp_code = "LESP")

####   19.) Ashy storm-petrel (Hydrobates homochroa) -- ASSP:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/ASSP_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/ASSP_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/ASSP_fall_predicted_density.tif
assp <- annual_species_density_function(directory = marine_bird_dir, spp_code = "ASSP")

####   20.) Black storm-petrel (Hydrobates melania) -- BLSP:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BLSP_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BLSP_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BLSP_fall_predicted_density.tif
blsp <- annual_species_density_function(directory = marine_bird_dir, spp_code = "BLSP")

####   21.) Northern fulmar (Fulmarus glacialis) -- NOFU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/NOFU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/NOFU_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/NOFU_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/NOFU_winter_predicted_density.tif
nofu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "NOFU")

####   22.) Murphy's petrel (Pterodroma ultima) -- MUPE:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/MUPE_spring_predicted_density.tif
mupe <- annual_species_density_function(directory = marine_bird_dir, spp_code = "MUPE")

####   23.) Cook's petrel (Pterodroma cookii) -- COPE:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/COPE_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/COPE_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/COPE_fall_predicted_density.tif
cope <- annual_species_density_function(directory = marine_bird_dir, spp_code = "COPE")

####   24.) Buller's shearwater (Ardenna bulleri) -- BULS:
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BULS_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BULS_fall_predicted_density.tif
buls <- annual_species_density_function(directory = marine_bird_dir, spp_code = "BULS")

####   25.) Pink-footed shearwater (Ardenna creatopus) -- PFSH:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PFSH_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PFSH_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PFSH_fall_predicted_density.tif
pfsh <- annual_species_density_function(directory = marine_bird_dir, spp_code = "PFSH")

####   26.) Black-vented shearwater (Puffinus opisthomelas) -- BVSH:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BVSH_spring_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BVSH_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BVSH_winter_predicted_density.tif
bvsh <- annual_species_density_function(directory = marine_bird_dir, spp_code = "BVSH")

####   30.) Brown pelican (Pelecanus occidentalis) -- BRPE:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BRPE_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BRPE_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BRPE_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BRPE_winter_predicted_density.tif
brpe <- annual_species_density_function(directory = marine_bird_dir, spp_code = "BRPE")

#####################################

#### Taxonomic Groups:
####   1.) Scoter species (surf scoter, white-winged scoter, black scoter) -- SUSC, WWSC, BLSC:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/SCOT_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/SCOT_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/SCOT_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/SCOT_winter_predicted_density.tif
scot <- annual_species_density_function(directory = marine_bird_dir, spp_code = "SCOT")

####   2.) Western (Aechmophorus occidentalis) / Clark's Grebe (Aechmophorus clarkii) -- WEGR, CLGR: 
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/WEGR-CLGR_spring_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/WEGR-CLGR_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/WEGR-CLGR_winter_predicted_density.tif
wegr_clgr <- annual_species_density_function(directory = marine_bird_dir, spp_code = "WEGR-CLGR")

####   3.) Phalarope species (red-necked phalarope, red phalarope) -- RNPH, REPH: 
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PHAL_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PHAL_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PHAL_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PHAL_winter_predicted_density.tif
phal <- annual_species_density_function(directory = marine_bird_dir, spp_code = "PHAL")

####   4.) Jaeger species:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/JAEG_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/JAEG_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/JAEG_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/JAEG_winter_predicted_density.tif
jaeg <- annual_species_density_function(directory = marine_bird_dir, spp_code = "JAEG")

###### ***NOTE: There are other jaeger data based on species
#####     Pomarine jaeger (Stercorarius pomarinus) -- POJA:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/POJA_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/POJA_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/POJA_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/POJA_winter_predicted_density.tif

#####     Parasitic jaeger (Stercorarius parasiticus) & long-tailed jaeger (Stercorarius longicaudus) -- PAJA, LTJA:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PAJA-LTJA_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PAJA-LTJA_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PAJA-LTJA_fall_predicted_density.tif

####   5.) Scripps's (Synthliboramphus scrippsi) / Guadeloupe (Synthliboramphus hypoleucus) / Craveri's murrelet (Synthliboramphus craveri) -- SCMU, GUMU, CRMU: 
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/SCMU-GUMU-CRMU_spring_predicted_density.tif
scmu_gumu_crmu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "SCMU-GUMU-CRMU")

####   6.) Herring (Larus smithsonianus) / Iceland gull (Larus glaucoides) -- HERG, ICGU: 
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/HERG-ICGU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/HERG-ICGU_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/HERG-ICGU_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/HERG-ICGU_winter_predicted_density.tif
herg_icgu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "HERG-ICGU")

####   7.) Western (Larus occidentalis) / glaucous-winged gull (Larus glaucescens) -- WEGU, GWGU:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/WEGU-WGWH-GWGU_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/WEGU-WGWH-GWGU_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/WEGU-WGWH-GWGU_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/WEGU-WGWH-GWGU_winter_predicted_density.tif
wegu_wgwh_gwgu <- annual_species_density_function(directory = marine_bird_dir, spp_code = "WEGU-WGWH-GWGU")

####   8.) Common (Sterna hirundo) / Arctic tern (Sterna paradisaea) -- COTE, ARTE: 
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/COTE-ARTE_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/COTE-ARTE_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/COTE-ARTE_fall_predicted_density.tif
cote_arte <- annual_species_density_function(directory = marine_bird_dir, spp_code = "COTE-ARTE")

####   9.) Royal (Thalasseus maximus) / elegant tern (Thalasseus elegans) -- ROYT, ELTE: 
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/ROYT-ELTE_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/ROYT-ELTE_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/ROYT-ELTE_fall_predicted_density.tif
royt_elte <- annual_species_density_function(directory = marine_bird_dir, spp_code = "ROYT-ELTE")

####   10.) Loon species (Gavia spp.) -- RTLO, PALO, COLO, YBLO:
####      ***NOTE: paper references 4 species: red-throated (Gavia stellata), Pacific (Gavia pacifica), common (Gavia immer), and yellow-billed (Gavia adamsii)
####      ***WARNING: Cornell does not record the yellow-billed loon as part of the Gaviidae family: https://www.allaboutbirds.org/guide/browse/taxonomy/Gaviidae
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/LOON_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/LOON_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/LOON_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/LOON_winter_predicted_density.tif
loon <- annual_species_density_function(directory = marine_bird_dir, spp_code = "LOON")

####   11.) Short-tailed (Ardenna tenuirostris) / sooty (Ardenna grisea) / flesh-footed shearwater (Ardenna carneipes) -- STTS, SOSH, FFSH:
####      ***WARNING: the short-tailed shearwater has different codes between the data download (STTS) and the paper (SRTS)
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/STTS-SOSH-FFSH_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/STTS-SOSH-FFSH_summer_predicted_density.tif
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/STTS-SOSH-FFSH_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/STTS-SOSH-FFSH_winter_predicted_density.tif
stts_sosh_ffsh <- annual_species_density_function(directory = marine_bird_dir, spp_code = "STTS-SOSH-FFSH")

#####################################

## Cormorants
####   27.) Brandt's cormorant (Phalacrocorax penicillatus) -- BRAC:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BRAC_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/BRAC_summer_predicted_density.tif
brac <- cormorant_species_function(directory = marine_bird_dir, spp_code = "BRAC")

####   28.) Pelagic cormorant (Phalacrocorax pelagicus) -- PECO:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PECO_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/PECO_summer_predicted_density.tif
peco <- cormorant_species_function(directory = marine_bird_dir, spp_code = "PECO")

####   29.) Double-crested cormorant (Phalacrocorax auritus) -- DCCO:
####      Spring: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/DCCO_spring_predicted_density.tif
####      Summer: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/DCCO_summer_predicted_density.tif
dcco <- cormorant_species_function(directory = marine_bird_dir, spp_code = "DCCO")

####   12.) Cormorant (Phalacrocorax spp.):
####      Fall: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CORM_fall_predicted_density.tif
####      Winter: https://www.nodc.noaa.gov/archive/arc0193/0242882/1.1/data/0-data/model_output_predictions/CORM_winter_predicted_density.tif
corm <- cormorant_species_function(directory = marine_bird_dir, spp_code = "CORM")


corm_annual_raster <- c(brac,
                        peco,
                        dcco,
                        corm) %>%
  terra::app(sum, na.rm = T)

# calculate annual total density
corm_annual_total <- terra::global(x = corm_annual_raster, fun = "sum", na.rm = T)
corm_annual_sum <- corm_annual_total$sum

# Normalize densities
## Divide the annual density data by the total summed data
corm_annual_norm <- corm_annual_raster[] / corm_annual_sum
## Set the normalized data back to the original annual density dataset
corm_normalize <- terra::setValues(corm_annual_raster, corm_annual_norm)

#####################################
#####################################

# Combine all species together to formulate Pacific species group
# ***Note: Pacific species raster needs to function as a list given
#          the weighted mean will apply species weights to values
#          within respective species raster
pacific_species <- terra::rast(list(anmu,
                                    assp,
                                    bfal,
                                    blki,
                                    blsp,
                                    bogu,
                                    brpe,
                                    buls,
                                    bvsh,
                                    caau,
                                    cagu,
                                    cate,
                                    comu,
                                    cope,
                                    corm_normalize,
                                    cote_arte,
                                    ftsp,
                                    heeg,
                                    herg_icgu,
                                    jaeg,
                                    laal,
                                    lesp,
                                    loon,
                                    mamu,
                                    mupe,
                                    nofu,
                                    pfsh,
                                    phal,
                                    pigu,
                                    rhau,
                                    royt_elte,
                                    sagu,
                                    scmu_gumu_crmu,
                                    scot,
                                    spsk,
                                    stts_sosh_ffsh,
                                    tupu,
                                    wegr_clgr,
                                    wegu_wgwh_gwgu))

#####################################
#####################################

# Calculate the weighted mean of Pacific marine bird specie
marine_bird_wt_mean <- terra::weighted.mean(pacific_species, w = species_weights)
plot(marine_bird_wt_mean)
hist(marine_bird_wt_mean)
terra::minmax(marine_bird_wt_mean)
marine_bird_wt_mean

# 800m setback provides the minimum setback required to get all hex grids to show
# may affect the z-membership function scores slightly
oregon_call_area_800m <- oregon_call_areas %>%
  # apply a 800m buffer to the Oregon call areas
  sf::st_buffer(dist = 800) %>%
  # change projection to match the weighted mean marine bird species raster
  sf::st_transform(crs = terra::crs(marine_bird_wt_mean))

marine_bird_wt_mean_call_area_800m <- marine_bird_wt_mean %>%
  # extract marine bird weighted mean data to expanded Oregon call areas
  terra::crop(oregon_call_area_800m,
              # mask data to the expanded Oregon call areas
              mask = T)

#####################################

# rescale marine bird values using a z-membership function
marine_bird_wt_mean_800m_zmf <- marine_bird_wt_mean_call_area_800m %>%
  zmf_function()

#####################################

# convert to polygon
oregon_marine_bird_wt_mean_800m_zmf_polygon <- terra::as.polygons(x = marine_bird_wt_mean_800m_zmf,
                                                         aggregate = F,
                                                         values = T) %>%
  # change to simple feature (sf)
  sf::st_as_sf() %>%
  # reproject data into a coordinate system (NAD 1983 UTM Zone 10N) that will convert units from degrees to meters
  sf::st_transform("EPSG:26910") %>%
  # create field called "layer" and populate with "marine bird"
  dplyr::mutate(layer = "marine bird") %>%
  # rename field "sum" as "wt_mean"
  dplyr::rename(wt_mean = sum)

#####################################
#####################################

# Marine sea bird hex grid
start2 <- Sys.time()
oregon_hex_marine_bird <- oregon_hex[oregon_marine_bird_wt_mean_800m_zmf_polygon, ] %>%
  # spatially join continental shelf values to Oregon hex cells
  sf::st_join(x = .,
              y = oregon_marine_bird_wt_mean_800m_zmf_polygon,
              join = st_intersects) %>%
  # select fields of importance
  dplyr::select(index, layer,
                wt_mean) %>%
  # group by the index values as there are duplicates
  dplyr::group_by(index) %>%
  # summarise the fisheries score values
  ## take the minimum value of the fisheries score for any that overlap
  ## ***Note: this will provide the most conservation given that low
  ##          values are less desirable
  dplyr::summarise(marine_bird_index = min(wt_mean))
print(Sys.time() - start2) # print how long it takes to create hex grid

# check plot
terra::plot(marine_bird_wt_mean_800m_zmf, col = colorRampPalette(c("#2C7BB6", "#ABD9E9", "#FFFFBF", "#FDAE61", "#D7191C"))(256), colNA = "gray50", main = "All species")

#####################################
#####################################

# Export data
## Natural resources submodel
sf::st_write(obj = oregon_hex_marine_bird, dsn = natural_resources_submodel, layer = paste0(region, "_hex_", layer), append = F)

##### Marine bird species geopackage (vector data)


#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate
