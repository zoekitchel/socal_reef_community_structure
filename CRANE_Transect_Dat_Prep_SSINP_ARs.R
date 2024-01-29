# CREATION DATE 27 Jul 2023

# AUTHOR: jtclaisse@cpp.edu

# PURPOSE: SSINP Artificial reef analysis
# Creates data tables of subsetted CRANE data
# Create a location table of AR sites (for use with AR habitat classification and mapping)

# TO DO - link metric table for AR modules from habitat classification and mapping to CRANE data


## Questions:
# Which Horseshoe Kelp sites are ARs or NRs?

# Source CRANE data prep --------------------------------------------------

source("~/Dropbox/VRG Files/R Code/Integrated Dive General/CRANE_data_prep.R")

CRANE_data_prep(SELECT = "SSINP_AR_NR_Comps", add_0s = "sp_0s") #dat_fish_t 337280 obs, 16 columns; #I need to create a new select

glimpse(dat_fish_t)
event_counts_check <- dat_event |> 
  count(Region, Site, DepthZone, TransectType)

#include current date in file name
write_csv(dat_fish_t, paste("data/processed_crane/dat_fish_t_",format(Sys.time(), "%d_%b_%Y"),".csv", sep = ""))

#include current date in file name
write_csv(dat_event, paste("data/processed_crane/dat_event_",format(Sys.time(), "%d_%b_%Y"),".csv", sep = ""))


# Extract AR Site Lat/Longs from dat_event (for Matt to use in GIS) -----------
library(tidyverse)

glimpse(dat_event)
dat_event_locations <- dat_event |>
  # filter(DepthZone == "ARM",
  #        TransectType == "Bottom",
  #        !str_starts(Site, "PVR")) |>
  dplyr::select(Site, SampleDate,SampleYear, Latitude, Longitude) |>
  distinct() |>
  mutate(SampleDate = dmy(SampleDate)) |>
  arrange(Site, SampleDate)|>
  mutate(SampleDate = format(SampleDate, "%d-%b-%Y"))

glimpse(dat_event_locations)

#include current date in file name
#write_csv(dat_event_locations, paste("Data Tables/SMB_AR_Lat_Longs_",format(Sys.time(), "%d_%b_%Y"),".csv", sep = ""))

# test out Zoe's enviro/habitat variables functions based on lat/longs --------

source("https://raw.githubusercontent.com/zoekitchel/CA_environmental_data/main/scripts/bathy_join_function.R")
#you will need the following packages loaded
library(sf)
# library(data.table)
# library(dplyr)
library(raster)
library(marmap) #to pull depth data

# must remove NAs or causes error with: 
# ETOPO = T
dat_event_locations <- dat_event_locations |>
  na.omit() |> 
  filter(Longitude > -180) #typo in Dat_event PVR 2B 20-Sep-2022 -188.3500 (should be -118)

#ETOPO = T causes error when run with island sites (I am assuming this is because it includes islands there is no depth data there?)

dat_event_locations

# ETOPO_test <- add_depth_columns(dat_event_locations, ETOPO = T, CDFW = F, USGS_socal=F, dist_200m = F)
# CDFW_test <- add_depth_columns(dat_event_locations, ETOPO = F, CDFW = T, USGS_socal=F, dist_200m = F)
# USGS_socal_test <- add_depth_columns(dat_event_locations, ETOPO = F, CDFW = F, USGS_socal=T, dist_200m = F)
# dist_200m_test <- add_depth_columns(dat_event_locations, ETOPO = F, CDFW = F, USGS_socal=F, dist_200m = T)

ALL4_add_depth_columns_test <- add_depth_columns(dat_event_locations, ETOPO = F, CDFW = T, USGS_socal=T, dist_200m = T)

# Test out in-situ Relief and Substrate functions ----------------------------
source("https://raw.githubusercontent.com/zoekitchel/CA_environmental_data/main/UPC_in_situ_habitat.R")
upc_relief <- get_relief(UPC_complete)
upc_substrate <- get_substrate(UPC_complete)
in_situ_metrics <- upc_relief |>  
  left_join(upc_substrate)

# correlation matrix of in-situ metrics
library(corrr)
in_situ_cor <- in_situ_metrics |> 
  dplyr::select(-c(1:7)) |> 
  correlate() |> 
  stretch() |> 
  arrange(desc(abs(r)))

