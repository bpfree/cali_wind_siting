# Clear workspace
rm(list = ls())

# Setup
################################################################################

# Packages
library(RSelenium)
library(tidyverse)
library(httr)

# Directories
outputdir <- "data/zz_miscellaneous"

# States
states <- c("North Carolina")

# Loop through
################################################################################

# Firefox profile
# Based on this link:
# https://yizeng.me/2014/05/23/download-pdf-files-automatically-in-firefox-using-selenium-webdriver/
fprof <- RSelenium::makeFirefoxProfile(list(
  browser.download.folderList = 2L,
  browser.download.dir = outputdir,
  browser.helperApps.neverAsk.saveToDisk = "application/pdf",
  pdfjs.disabled = TRUE,
  plugin.scan.plid.all = FALSE,
  plugin.scan.Acrobat = "99.0"))

# Launch RSelenium server and driver
rD <- rsDriver(browser="firefox",
               version = "latest",
               chromever = NULL,
               geckover = "latest",
               verbose = T)

remDr <- rD[["client"]]
remDr$open(silent=T)

# Base URL
base_url <- "https://www.roc.noaa.gov/WSR88D/Program/SiteID.aspx"

# Navigate to page
remDr$navigate(base_url)


# Click "Advanced Search" to search by state
advanced_search <-remDr$findElement(using = "link text",
                                    value = "Advanced Search")
advanced_search$clickElement()

remDr$maxWindowSize()
Sys.sleep(20)


check_url <- advanced_search$getCurrentUrl()
check_url

source <- remDr$getPageSource()[[1]]

# Choose the states
webElem <- remDr$findElements("css", "iframe")
remDr$switchToFrame(webElem[[1]])


state <-remDr$findElement(using = "name",
                          value = "lstState")
state$sendKeysToElement(list("North Carolina",
                             key = "enter"))
Sys.sleep(10)


# Select pertinent data fields
## state field
state_field <- remDr$findElement(using = "id",
                                 value = "FrmFld12")
state_field$clickElement()
  
latitude <- remDr$findElement(using = "id",
                               value = "FrmFld18")
latitude$clickElement()

longitude <- remDr$findElement(using = "id",
                               value = "FrmFld19")
longitude$clickElement()


# Get search results
search <- remDr$findElement(using = "name",
                            value = "Submit")
search$clickElement()


# Extract table
table <- remDr$findElement(using = "xpath",
                           value = "/html/body/center/table")

check_url <- search$getCurrentUrl()

source <- remDr$getPageSource()[[1]]


table_clean <- rvest::read_html(source) %>%
  # obtain the table
  rvest::html_element(css = "table") %>%
  # read the table to get it as a data frame
  rvest::html_table() %>%
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

remDr$close()
rD$server$stop()


