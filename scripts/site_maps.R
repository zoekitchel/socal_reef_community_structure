# CREATION DATE 25 March 2024

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Site Map

#############################
##Setup
#############################
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

#add type column
lat_lon_site_fix[,`Reef type` := factor(
  ifelse(DepthZone == "ARM", "Artificial reef","Natural reef"))]

#change col names
colnames(lat_lon_site_fix) <- c("Site","DepthZone","Latitude","Longitude")

#convert to spatial points
lat_lon_site.sf <- st_as_sf(lat_lon_site_fix,
                            coords = c("Longitude","Latitude"),
                            crs = 4326)

#avg values
#mean lat and lon
lat_lon_site_fix[,avg_lon := mean(Longitude,na.rm = T),Site][,avg_lat := mean(Latitude,na.rm = T),Site]
lat_lon_site_fix.r <- unique(lat_lon_site_fix[,.(Site, `Reef type`, avg_lon, avg_lat)])

##############################################type_sum()
#Map
#######################
#set square from which to extract bathy data from NOAA server
bathy_VRG <- getNOAA.bathy(min(lat_lon_site_fix.r$avg_lon)-2, max(lat_lon_site_fix.r$avg_lon)+2, min(lat_lon_site_fix.r$avg_lat)-0.5, max(lat_lon_site_fix.r$avg_lat)+0.5, resolution = 0.000001) #bathymetry matrix

#map of bathymetry
site_map <- autoplot.bathy(bathy_VRG, geom=c("raster"
                                 # ,"contour" exclude contour
), coast = F) +
  scale_fill_etopo(breaks = c(-12000,-6000, 0, 6000, 10000), labels = c(12, 6, 0, 6, 10), #from marmap, great way to visualize land and water instead of 'world' object
                   guide = NULL) +
  geom_point(data = lat_lon_site_fix.r, aes(x = avg_lon, y = avg_lat, color = `Reef type`), size = 2, alpha = 0.5) +
  scale_color_manual(values = c("darkturquoise","brown1"))+
  labs(y = "Latitude", x = "Longitude", fill = "Elevation/Depth\n(1000s of m)") +
 # scale_x_continuous(breaks = c(-120:-117), labels = c("120˚W" ,"119˚W" ,"118˚W" ,"117˚W")) +
  scale_y_continuous(breaks = c(33:34),labels = c("33˚N" ,"34˚N" )) +
  coord_sf(xlim = c(-119.1,-116.75), ylim = c(32.6, 34.3), expand = F) +
  theme_classic() +
  theme(axis.title = element_blank(),
        legend.position = c(0.15, 0.15),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.background = element_rect(colour = "black", fill = "transparent"))

ggsave(site_map, filename = "site_map.jpg", path = file.path("figures"), width = 10, height = 7, units = "in")
