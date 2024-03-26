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

########################
##Load data
########################
dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))

#avg lat lon for each site
dat_event.r[,avg_lat := mean(Latitude, na.rm = T),.(Site)][,avg_lon := mean(Longitude, na.rm = T),.(Site)][,`Reef type` := factor(
  ifelse(DepthZone == "ARM", "Artificial reef","Natural reef"))]

#unique values
unique_lat_lon <- unique(dat_event.r[,.(Site, `Reef type`, avg_lon, avg_lat)])

##############################################type_sum()
#Map
#######################
#set square from which to extract bathy data from NOAA server
bathy_VRG <- getNOAA.bathy(min(unique_lat_lon$avg_lon)-2, max(unique_lat_lon$avg_lon)+2, min(unique_lat_lon$avg_lat)-0.5, max(unique_lat_lon$avg_lat)+0.5, resolution = 0.000001) #bathymetry matrix

#map of bathymetry
site_map <- autoplot.bathy(bathy_VRG, geom=c("tile"
                                 # ,"contour" exclude contour
), coast = F) +
  scale_fill_etopo(breaks = c(-12000,-6000, 0, 6000, 10000), labels = c(12, 6, 0, 6, 10), #from marmap, great way to visualize land and water instead of 'world' object
                   guide = NULL) +
  geom_point(data = unique_lat_lon, aes(x = avg_lon, y = avg_lat, color = `Reef type`), size = 2, alpha = 0.5) +
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
