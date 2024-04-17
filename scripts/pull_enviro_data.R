# CREATION DATE 14 Jan 2024
# MODIFIED DATE 14 April 2024

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Pull environmental data, both in situ and external

#############################
##Setup
#############################
library(ggplot2)
library(data.table)
library(dplyr)
library(sf)
library(raster)
library(sdmpredictors)
library(stringr)
library(rvest)

#pull in functions
source("https://raw.githubusercontent.com/zoekitchel/CA_environmental_data/main/UPC_in_situ_habitat.R")


#############################
#Get correct lat lon values
#############################

#Chelsea recommends using Lat/Lon from "2023 Dive Site Priority" instead of All Sites File
dive_site_priority_list <- fread("dive_site_priority_list.csv")

#extract depth zone and site from location column
dive_site_priority_list[,DepthZone := word(Location,-1)][,Site := word(Location, start = 1, end = -2)]

#change column names
colnames(dive_site_priority_list) <- c("Location","Latitude_fix","Longitude_fix","DepthZone","Site")

#limit to lat, lon, depthzone, site
dive_site_priority_list.r <- dive_site_priority_list[,.(Site, DepthZone, Latitude_fix, Longitude_fix)]

#load dat_event.r
dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))

#reduce to unique lat, lon, site, depth zone
lat_lon_site_orig <- unique(dat_event.r[,.(Site, DepthZone, Latitude, Longitude)])

#merge event data with fixed lat and lon from dive site priority list
lat_lon_site_fix <- lat_lon_site_orig[dive_site_priority_list.r, on = c("Site","DepthZone")]

#Delete any rows without values, and use this as key for Site, Lat and Long
lat_lon_site_fix <- lat_lon_site_fix[complete.cases(lat_lon_site_fix),]

#only keep unique values
lat_lon_site_fix <- unique(lat_lon_site_fix) #245 sites

#delete old lat lon columns from all sites
lat_lon_site_fix <- lat_lon_site_fix[,c(1,2,5,6)]

#change col names
colnames(lat_lon_site_fix) <- c("Site","DepthZone","Latitude","Longitude")

#convert to spatial points
lat_lon_site.sf <- st_as_sf(lat_lon_site_fix,
                            coords = c("Longitude","Latitude"),
                            crs = 4326)

#save csv
fwrite(lat_lon_site_fix, file.path("data","processed_crane","lat_lon_site_fix.csv"))

#############
#Prep in situ habitat characteristics
#############

#pull in UPC_complete data for substrate and relief calculations
UPC_complete <- readRDS(file.path("data","full_crane","UPC_complete.rds"))

#######Calculate average macrocystis stipe density and plant density per site
#first, only macrocystis
dat_kelp_site_averages <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages.rds"))

#stipes
stipe_density_bysite <- dat_kelp_site_averages[SpeciesGroupF == "giant kelp stipes",.(Site,DepthZone, mean_density_m2)]

colnames(stipe_density_bysite) <- c("Site","DepthZone","giantkelp_stipe_density_m2")

#kelp plants
macro_density_bysite <- dat_kelp_site_averages[BenthicReefSpecies == "Macrocystis pyrifera"][,.(Site, DepthZone, mean_density_m2)]

colnames(macro_density_bysite) <- c("Site","DepthZone","giantkelp_density_m2")


###########relief
relief <- get_relief(UPC_complete)

#average to get single value per Site and Depth Zone
relief <- relief[,.(Site,DepthZone, Relief_index, Relief_SD, Relief_simpson)] #limit to columns we need

relief_bysite <-  relief[ , lapply(.SD, mean) , by=c("Site", "DepthZone")] #single value for each site

###########substrate
substrate <- get_substrate(UPC_complete)

#average to get single value per Site and Depth Zone
substrate <- substrate[,.(Site,DepthZone, Substrate_index, Substrate_SD, Substrate_simpson)] #limit to columns we need

substrate_bysite <-  substrate[ , lapply(.SD, mean) , by=c("Site", "DepthZone")] #single value for each site


##################################################
#Prep large scale habitat characteristics
##################################################
#event data
source("https://raw.githubusercontent.com/zoekitchel/CA_environmental_data/main/scripts/bathy_join_function.R")

#######TEMP

BO_sst <- load_layers(layercodes = c("BO_sstmax", "BO_sstmean", "BO_sstmin", "BO_sstrange") ,
                      equalarea=FALSE, rasterstack=TRUE) 



#reduce to CA only
CA_sst <- crop(BO_sst, extent(min(dat_event.r$Longitude, na.rm = T)-0.1, max(dat_event.r$Longitude, na.rm = T)+0.1,
                              min(dat_event.r$Latitude, na.rm = T)-0.1, max(dat_event.r$Latitude, na.rm = T)+0.1))

rm(BO_sst)

#create multifocal function to apply focal function to each layer of the raster brick

multiFocal <- function(x, w=matrix(1, nr=3, nc=3), ...) {
  
  if(is.character(x)) {
    x <- brick(x)
  }
  # The function to be applied to each individual layer
  fun <- function(ind, x, w, ...){
    focal(x[[ind]], w=w, ...)
  }
  
  n <- seq(nlayers(x))
  list <- lapply(X=n, FUN=fun, x=x, w=w, ...)
  
  out <- stack(list)
  return(out)
}


#fill in two missing cells with neighborhood averages
#important to note that this does also add erroneous cells with values on the edge, but this does not impact our extractions because these are points on land and no sampling points are on land
CA_sst.t <- multiFocal(CA_sst, w = matrix(1, 3, 3), fun = "mean", NAonly = TRUE, na.rm = T)

#make points into sf in same projection
VRG_lat_lon_only <- copy(dat_event.r)

#we only need unique lat lon coordinates, so we exclude time and reduce to unique longitude and latitude values
VRG_lat_lon_only <- VRG_lat_lon_only[complete.cases(Longitude, Latitude),.(Site, DepthZone,Longitude, Latitude)]


#convert list of lat lons to a simple feature object
VRG_lat_lon_only.sf <- st_as_sf(VRG_lat_lon_only, coords = c("Longitude","Latitude"), crs = 4326)

#transform to crs of raster stack
VRG_lat_lon_only.t <- st_transform(VRG_lat_lon_only.sf, crs = st_crs(CA_sst))

#add new columns with extracted temperature data
#maximum
VRG_lat_lon_only[,BO_sstmax := raster::extract(CA_sst.t[[1]], VRG_lat_lon_only.t)]

#mean
VRG_lat_lon_only[,BO_sstmean := raster::extract(CA_sst.t[[2]], VRG_lat_lon_only.t)]

#minimum
VRG_lat_lon_only[,BO_sstmin := raster::extract(CA_sst.t[[3]], VRG_lat_lon_only.t)]

#range
VRG_lat_lon_only[,BO_sstseas := raster::extract(CA_sst.t[[4]], VRG_lat_lon_only.t)]

#save
saveRDS(VRG_lat_lon_only, file.path("data","enviro_predictors","VRG_lat_lon_only.rds"))

#######TEMP san diego 1km satellite
#Link: http://spg-satprojects.ucsd.edu
#Link: https://spg-satdata.ucsd.edu
#Metadata for HDF files

#Scrape links to monthly temp and chlorophyll values from UCSD website

#first, navigate to year page
folder_urls <- readLines(file.path("data","enviro_predictors","ucsd_data_download_1km_chl_sst_2016_2023.txt"))

#empty data.table
full_link.dt <- data.table()

#make list of all links to download
for(i in 1:length(folder_urls)){
  url <- folder_urls[i]
  
  variable <- ifelse(str_detect(url,"chl"),"chl","sst")
  year <- as.numeric(substr(url,30,33))
 
   #individual month links using rvest
  monthly_link_list <- url %>%
    read_html() %>%
    html_nodes("table") %>% html_nodes("tr") %>% html_nodes("a") %>%
    html_attr("href") #identify all links on page
  
  monthly_link_list.r <- monthly_link_list |>
    keep(~ str_detect(.x, "comp.hdf"))
  
  #add full url to all links
  monthly_link_list.full <- paste0(url,monthly_link_list.r)
  
  #make datatable
  subset.dt <- data.table(variable = variable, year = year, month = seq(1,12),
                          file_name = monthly_link_list.r, link_long = monthly_link_list.full)
  
  full_link.dt <- rbind(full_link.dt, subset.dt)
  
}

#empty data table to fill with site, and month specific data
lat_lon_site_variable_full <- data.table()

#load up hdf key with lat lon
hdf_key <- stack("/Users/kitchel/Downloads/cal_aco_3840_Latitude_Longitude.hdf")

#trim to study area
crop_ext <- extent(1617,2069,1727,2204)

hdf_key.c <- crop(hdf_key, crop_ext)

#convert into data table
hdf_key.xyz <- data.table(rasterToPoints(hdf_key.c))

colnames(hdf_key.xyz) <- c("x","y","latitude","longitude")

#download and process each kmz file and populate data table with chlorophyll and temp data
for(i in 1:nrow(full_link.dt)){
  temp <- tempdir()
  download.file(full_link.dt[i,link_long], file.path(temp, "temp.hdf"))
  hdf <- raster(file.path(temp, "temp.hdf"))
  
  #trim to study area
  crop_ext <- extent(1617,2069,1727,2204)
  
  hdf.c <- crop(hdf, crop_ext)
  
  #change 255 values to NA ()
  hdf.c.re <- reclassify(hdf.c, cbind(255, NA)) #"values of 0 and 255 are considered invalid"
  
  #convert into data table
  hdf.xyz <- data.table(rasterToPoints(hdf.c.re))
  
  colnames(hdf.xyz) <- c("x","y","value")

  #merge with key
  hdf_merge <- hdf.xyz[hdf_key.xyz, on = c("x","y")]
  
  #set up empty raster
  e <- extent(min(hdf_merge$longitude), max(hdf_merge$longitude), min(hdf_merge$latitude),max(hdf_merge$latitude))
  
  r <- raster(e,ncol = 452, nrow = 477)
  
  #then, rasterize xyz
  new_raster <- rasterize(x = hdf_merge[,.(longitude,latitude)], y = r, field = hdf_merge[,.(value)], fun = mean, na.rm = T)
  
  #set crs
  crs(new_raster) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
  
  #match crs of lat lon to crs of raster
  lat_lon_site.t <- st_transform(lat_lon_site.sf, crs(new_raster))
  
  value <- raster::extract(new_raster,lat_lon_site.t) #21 lat lon from VRG are missing! check these some time!
  
  #copy input lat lon 
  lat_lon_site_variable <- copy(lat_lon_site_fix)
  
  #fill new columns
  lat_lon_site_variable[,year:=full_link.dt[i,year]][,month:=full_link.dt[i,month]][,variable:= full_link.dt[i,variable]][,value:=  value]

  #add to full data table
  lat_lon_site_variable_full <- rbind(lat_lon_site_variable_full,   lat_lon_site_variable)
  
  print(paste0(i," out of ",nrow(full_link.dt)))
  
}

#Note that values use 1 byte per pixel with standard scaling. Linear scaling is used for SST and logarithmic scaling for Chl: https://spg-satdata.ucsd.edu/Readme.htm
#so, conversions
#SST (deg C) = 0.15 * PV - 3.0
lat_lon_site_variable_full[variable == "sst",value_adj := 0.15*value-3.0]

#Chl (mg m-3) = 10^(0.015 * PV - 2.0), i.e. 10 to the power of 0.015 * PV - 2.0
lat_lon_site_variable_full[variable == "chl",value_adj :=10^(0.015*value-2.0)]

#save as csv
fwrite(lat_lon_site_variable_full, file.path("data","enviro_predictors","lat_lon_site_variable_full.csv"))
#lat_lon_site_variable_full <- fread(file.path("data","enviro_predictors","lat_lon_site_variable_full.csv"))

#mean_value per site across all months in time series, also across DepthZones, because the resolution is 1km, so we wouldn't expect depth zones to be different
lat_lon_site_variable_full[, mean_value_adj := mean(value_adj, na.rm = T),.(Site, variable)]

#min and max value per site across all months in time series, also across DepthZones. Will first take min of each year, and then mean of minimums. Same for max
  #mean annual minimum
  lat_lon_site_variable_full[, min_annual_value_adj := min(value_adj, na.rm = T),.(Site, variable, year)][,min_annual_value_adj := ifelse(is.infinite(min_annual_value_adj), NA, min_annual_value_adj)]
  lat_lon_site_variable_full[, min_value_adj := mean(min_annual_value_adj, na.rm = T),.(Site, variable)]
  #mean annual maximum
  lat_lon_site_variable_full[, max_annual_value_adj := max(value_adj, na.rm = T),.(Site, variable, year)][,max_annual_value_adj := ifelse(is.infinite(max_annual_value_adj), NA, max_annual_value_adj)]
  lat_lon_site_variable_full[, max_value_adj := mean(max_annual_value_adj, na.rm = T),.(Site, variable)]

#reduce to one row per site
lat_lon_site_variable_full.r <- unique(lat_lon_site_variable_full[,.(Site, variable, mean_value_adj, min_value_adj, max_value_adj)])

#split into temperature and chlorophyll
sst_bysite <- lat_lon_site_variable_full.r[variable == "sst"][,mean_sst_C := mean_value_adj][,max_sst_C := max_value_adj][,min_sst_C := min_value_adj][,.(Site, mean_sst_C, max_sst_C, min_sst_C)]
chl_bysite <- lat_lon_site_variable_full.r[variable == "chl"][,mean_chl_mg_m3 := mean_value_adj][,max_chl_mg_m3 := max_value_adj][,min_chl_mg_m3 := min_value_adj][,.(Site, mean_chl_mg_m3, max_chl_mg_m3, min_chl_mg_m3)]

#######Distance to 200 m isobath

distance_200mbathy_bysite <- add_depth_columns(lat_lon_site_fix, ETOPO = F, CDFW = F, USGS_socal=F, dist_200m = T)

#remove lat lon columns
distance_200mbathy_bysite <- unique(distance_200mbathy_bysite[,.(Site, DepthZone, dist_200m_bath)])

###############################################
#merge all site level environmental variables
###############################################
#insitu
substrate_bysite
all_env_lat_lon <- substrate_bysite[lat_lon_site_fix, on = c("Site","DepthZone")]
relief_bysite
all_env_lat_lon <- relief_bysite[all_env_lat_lon, on = c("Site","DepthZone")]
macro_density_bysite
all_env_lat_lon <- macro_density_bysite[all_env_lat_lon, on = c("Site","DepthZone")]
stipe_density_bysite
all_env_lat_lon <- stipe_density_bysite[all_env_lat_lon, on = c("Site","DepthZone")]


#large scale
distance_200mbathy_bysite
all_env_lat_lon <- distance_200mbathy_bysite[all_env_lat_lon, on = c("Site","DepthZone")]
sst_bysite
all_env_lat_lon <- sst_bysite[all_env_lat_lon, on = c("Site")]
chl_bysite
all_env_lat_lon <- chl_bysite[all_env_lat_lon, on = c("Site")]

#save
fwrite(all_env_lat_lon,file.path("data","enviro_predictors","all_env_lat_lon.csv"))

###############################################
#plot by colors to make sure that these variables look right
###############################################

#mean value per site (get rid of depth zones)
all_env_lat_lon_nodepthzone <- all_env_lat_lon[,c(1:7,9:19)]
all_env_lat_lon.r <- all_env_lat_lon_nodepthzone[, lapply(.SD, mean, na.rm = T), by=c("Site")] #some values look into this 

library(ggplot2)
library(data.table)
library(dplyr)
library(vegan)
library(ggvegan)
library(labdsv)
library(cowplot)
library(marmap) #to pull depth data
library(rasterVis)
library(ggspatial)
library(ggmap) #background google map
library(sp)
library(RColorBrewer)

##############################################type_sum()
#Map
#######################
#set square from which to extract bathy data from NOAA server
bathy_VRG <- getNOAA.bathy(min(unique_lat_lon$avg_lon)-2, max(unique_lat_lon$avg_lon)+2, min(unique_lat_lon$avg_lat)-0.5, max(unique_lat_lon$avg_lat)+0.5, resolution = 0.000001) #bathymetry matrix

#map of bathymetry w kelp density
kelp_site_map <- autoplot.bathy(bathy_VRG, geom=c("tile"
                                             # ,"contour" exclude contour
), coast = F) +
  scale_fill_etopo(breaks = c(-12000,-6000, 0, 6000, 10000), labels = c(12, 6, 0, 6, 10), #from marmap, great way to visualize land and water instead of 'world' object
                   guide = NULL) +
  geom_point(data = all_env_lat_lon.r, aes(x = Longitude, y = Latitude, size = giantkelp_density_m2, color = giantkelp_density_m2), alpha = 0.8) +
  scale_size(range = c(1,3.5)) +
  scale_color_gradient(low = "white",high = "darkviolet") +
  labs(y = "Latitude", x = "Longitude", fill = "Elevation/Depth\n(1000s of m)") +
  # scale_x_continuous(breaks = c(-120:-117), labels = c("120˚W" ,"119˚W" ,"118˚W" ,"117˚W")) +
  scale_y_continuous(breaks = c(33:34),labels = c("33˚N" ,"34˚N" )) +
  coord_sf(xlim = c(-119.1,-116.75), ylim = c(32.6, 34.3), expand = F) +
  theme_classic() +
  theme(axis.title = element_blank())

ggsave(kelp_site_map, filename = "kelp_site_map.jpg", path = file.path("figures"), width = 10, height = 7, units = "in")

#map of bathymetry w temperature
temperature_site_map <- autoplot.bathy(bathy_VRG, geom=c("tile"
                                                  # ,"contour" exclude contour
), coast = F) +
  scale_fill_etopo(breaks = c(-12000,-6000, 0, 6000, 10000), labels = c(12, 6, 0, 6, 10), #from marmap, great way to visualize land and water instead of 'world' object
                   guide = NULL) +
  geom_point(data = all_env_lat_lon.r, aes(x = Longitude, y = Latitude, color = mean_sst_C), alpha = 0.8, size = 3) +
  scale_color_viridis(option = "C") +
  labs(y = "Latitude", x = "Longitude", fill = "Elevation/Depth\n(1000s of m)") +
  # scale_x_continuous(breaks = c(-120:-117), labels = c("120˚W" ,"119˚W" ,"118˚W" ,"117˚W")) +
  scale_y_continuous(breaks = c(33:34),labels = c("33˚N" ,"34˚N" )) +
  coord_sf(xlim = c(-119.1,-116.75), ylim = c(32.6, 34.3), expand = F) +
  theme_classic() +
  theme(axis.title = element_blank())

ggsave(temperature_site_map, filename = "temperature_site_map.jpg", path = file.path("figures"), width = 10, height = 7, units = "in")

#map of bathymetry w/ 200m distance
distance_200_site_map <- autoplot.bathy(bathy_VRG, geom=c("tile"
                                                  # ,"contour" exclude contour
), coast = F) +
  scale_fill_etopo(breaks = c(-12000,-6000, 0, 6000, 10000), labels = c(12, 6, 0, 6, 10), #from marmap, great way to visualize land and water instead of 'world' object
                   guide = NULL) +
  geom_point(data = all_env_lat_lon.r, aes(x = Longitude, y = Latitude, color = dist_200m_bath)) +
  scale_color_viridis(option = "turbo", direction = -1) +
  labs(y = "Latitude", x = "Longitude", fill = "Elevation/Depth\n(1000s of m)") +
  # scale_x_continuous(breaks = c(-120:-117), labels = c("120˚W" ,"119˚W" ,"118˚W" ,"117˚W")) +
  scale_y_continuous(breaks = c(33:34),labels = c("33˚N" ,"34˚N" )) +
  coord_sf(xlim = c(-119.1,-116.75), ylim = c(32.6, 34.3), expand = F) +
  theme_classic() +
  theme(axis.title = element_blank())

ggsave(distance_200_site_map, filename = "distance_200_site_map.jpg", path = file.path("figures"), width = 10, height = 7, units = "in")



