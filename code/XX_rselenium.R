###################################
### 1. Webscrape Data -- NEXRAD ###
###################################

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
               RSelenium,
               sf,
               sp,
               stringr,
               terra, # is replacing the raster package
               tidyr,
               tidyverse)

#####################################
#####################################

# Set directories
## output directory
outputdir <- "data/a_raw_data"

## NEXRAD geopackage
nexrad_gpkg <- "data/a_raw_data/nexrad_sites.gpkg"

#####################################
#####################################

# Coastal states and regions
## All US coastal states and some other areas of interest (e.g., Guam, Korea, etc.)
states <- c("Alabama",
            "Alaska",
            #"Azores", # nothing returned
            "California",
            "Delaware",
            "Florida",
            "Georgia",
            "Guam",
            "Hawaii",
            "Japan",
            "Korea", # this is South Korea
            "Louisiana",
            "Maine",
            "Maryland",
            "Massachusetts",
            "Mississippi",
            "North Carolina",
            #"New Hampshire", # no NEXRAD sites
            "New Jersey",
            "New York",
            "Oregon",
            "Puerto Rico",
            #"Rhode Island", # no NEXRAD sites
            "South Carolina",
            "Texas",
            "Virginia",
            "Washington")

#####################################

# Zoom out Firefox (taken from here: https://stackoverflow.com/questions/60139836/zoom-out-of-website-when-using-rselenium-without-changing-page-size-resolution)
## This allows for all elements to be shown on the page
zoom_firefox <- function(client, percent){
  store_page <- client$getCurrentUrl()[[1]]
  client$navigate("about:preferences")
  webElem <- client$findElement("css", "#defaultZoom")
  webElem$clickElement()
  webElem$sendKeysToElement(list(as.character(percent)))
  webElem$sendKeysToElement(list(key = "return"))
  client$navigate(store_page)
}

#####################################
#####################################

# Webscrape set-up
## Process uses RSelenium package (learn more about basics here: https://cran.r-project.org/web/packages/RSelenium/vignettes/basics.html)
### Another helpful tutorial: https://joshuamccrain.com/tutorials/web_scraping_R_selenium.html
### Firefox profile (based on this link: https://yizeng.me/2014/05/23/download-pdf-files-automatically-in-firefox-using-selenium-webdriver/)
fprof <- RSelenium::makeFirefoxProfile(list(
  # detail level for download (0 = Desktop, 1 = systems default downloads location, 2 = custom folder.)
  browser.download.folderList = 2L,
  # location for data download
  browser.download.dir = outputdir,
  # stores a comma-separated list of MIME types to save to disk without asking what to use to open the file
  browser.helperApps.neverAsk.saveToDisk = "application/pdf",
  # disable PDF javascript so that PDFs are not displayed
  pdfjs.disabled = TRUE,
  # turn off scan and loading of any additionally added plug-ins
  plugin.scan.plid.all = FALSE,
  # high number defined for version of Adobe Acrobat
  plugin.scan.Acrobat = "99.0"))

#####################################

# Launch RSelenium server and driver
rD <- RSelenium::rsDriver(browser="firefox",
                          # set which version of browser
                          version = "latest",
                          # Chrome version (turn off as Firefox will be used)
                          chromever = NULL,
                          # set which version of Gecko to use
                          geckover = "latest",
                          # status messages
                          verbose = TRUE,
                          # populate with the Firefox profile
                          extraCapabilities = fprof)

## Remote driver
remDr <- rD[["client"]]
remDr$open(silent = TRUE)

## Set client (to be used in the zoom function)
client <- rD$client

#####################################
#####################################

# Obtain NEXRAD sites
## create reference table
nexrad_table <- data.frame(state = character(),
                           nexrad_sitename = character(),
                           site_id = character(),
                           agency = character(),
                           equip = character(),
                           lon_dd = numeric(),
                           lat_dd = numeric())

#####################################

# Loop that navigates to NEXRAD page (source: https://www.roc.noaa.gov/WSR88D/Program/SiteID.aspx)
# and scrapes the site data for the states and regions of interest

#i <- 3
for(i in 1:length(states)){
  start2 <- Sys.time()
  
  # Base URL
  base_url <- "https://www.roc.noaa.gov/WSR88D/Program/SiteID.aspx"
  
  # Navigate to page
  remDr$navigate(base_url)
  
  #####################################
  
  # Prepare window
  remDr$maxWindowSize()
  ## Zoom to 80%
  zoom_firefox(client, 80)
  
  #####################################
  
  # Click "Advanced Search" to search by state
  advanced_search <-remDr$findElement(using = "link text",
                                      value = "Advanced Search")
  advanced_search$clickElement()
  Sys.sleep(5)
  
  #####################################
  
  # Define web elements as an iframe so other elements can get selected
  webElem <- remDr$findElements("css", "iframe")
  ## Make the iframe the location to find the other elements
  remDr$switchToFrame(webElem[[1]])
  
  #####################################
  
  # Choose the states
  state <-remDr$findElement(using = "name",
                            value = "lstState")
  state$sendKeysToElement(list(states[i],
                               key = "enter"))
  Sys.sleep(5)
  
  # Select pertinent data fields
  ## state field
  state_field <- remDr$findElement(using = "id",
                                   value = "FrmFld12")
  Sys.sleep(3)
  state_field$clickElement()
  Sys.sleep(2)
  
  ## latitude field
  latitude <- remDr$findElement(using = "id",
                                value = "FrmFld18")
  Sys.sleep(3)
  latitude$clickElement()
  Sys.sleep(2)
  
  ## longitude field
  longitude <- remDr$findElement(using = "id",
                                 value = "FrmFld19")
  Sys.sleep(3)
  longitude$clickElement()
  Sys.sleep(2)
  
  #####################################
  
  # Get search results
  search <- remDr$findElement(using = "name",
                              value = "Submit")
  search$clickElement()
  
  #####################################
  
  # Extract table
  table <- remDr$findElement(using = "xpath",
                             value = "/html/body/center/table")
  
  # check_url <- search$getCurrentUrl()
  
  source <- remDr$getPageSource()[[1]]
  
  #####################################
  
  # Make table from the results of the search
  ## Read HTML page to create the table
  table_clean <- rvest::read_html(source) %>%
    # obtain the table
    rvest::html_element(css = "table") %>%
    # read the table to get it as a data frame
    rvest::html_table() %>%
    
    # set as a data frame
    as.data.frame() %>%
    
    # clean data table
    ## make first row the names
    janitor::row_to_names(row_number = 1) %>%
    ## make names all lowercase
    janitor::clean_names() %>%
    # separate out latitude data components
    tidyr::separate(latitude,
                    into = c("lat_d", "lat_m", "lat_s"),
                    sep = " ",
                    remove = T,
                    convert = T) %>%
    # separate out longitude data components
    tidyr::separate(longitude,
                    into = c("lon_d", "lon_m", "lon_s"),
                    sep = " ",
                    remove = T,
                    convert = T) %>%
    # remove any sites without longitude and latitude data
    na.omit() %>%
    # remove "+" from lat_d
    dplyr::mutate(lat_d = str_replace(lat_d, "\\+","")) %>%
    # make longitude and latitude values numeric to calculate decimal degrees
    dplyr::mutate_at(c("lat_d", "lon_d"),
                     as.numeric) %>%
    
    # convert to longitude and latitude into decimal degrees (degrees + minutes / 60 + seconds / (60 * 60))
    # ***Note: longitude values are multiplied by -1 as they are in the west of the Prime Meridian
    dplyr::mutate(lon_dd = -1 * (-1 * lon_d + lon_m /60 + lon_s/60^2),
                  lat_dd = lat_d + lat_m /60 + lat_s/60^2) %>%
    # select fields of interest
    dplyr::select(state,
                  nexrad_sitename,
                  site_id,
                  agency,
                  equip,
                  lon_dd,
                  lat_dd) %>%
    
    # convert to simple feature
    sf::st_as_sf(coords = c("lon_dd", "lat_dd"),
                 # set the coordinate reference system to WGS84
                 crs = 4326) %>% # EPSG 4326 (https://epsg.io/4326)
    # reproject the coordinate reference system to match study area data (EPSG:5070)
    sf::st_transform("EPSG:5070") # EPSG 5070 (https://epsg.io/5070)
  
  nexrad_table <- rbind(nexrad_table, table_clean)
  
  print(paste("Iteration", i, "of", length(states), "takes", Sys.time() - start2, units(Sys.time() - start2), "to complete creating and adding", states[i], "data to dataframe and export as plot", sep = " "))
}

#####################################
#####################################

# Close RSelenium servers
remDr$close()
rD$server$stop()

#####################################
#####################################

# Export data
sf::st_write(obj = nexrad_table, dsn = nexrad_gpkg, layer = "nexrad_sites", append = F)

#####################################
#####################################

# calculate end time and print time difference
print(Sys.time() - start) # print how long it takes to calculate