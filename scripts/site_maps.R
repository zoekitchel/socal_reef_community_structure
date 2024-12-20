# CREATION DATE 25 March 2024
# MODIFIED DATE 19 August 2024

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
library(sf)
library(stringr)
library(rnaturalearth)

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
lat_lon_site_fix <- unique(lat_lon_site_fix) #93 sites

#delete old lat lon columns from all sites
lat_lon_site_fix <- lat_lon_site_fix[,c(1,2,5,6)]

#add type column
lat_lon_site_fix[,`Reef type` := factor(
  ifelse(DepthZone == "ARM", "Artificial reef","Natural reef"))]

#change col names
colnames(lat_lon_site_fix) <- c("Site","DepthZone","Latitude","Longitude", "Reef type")

#set factor order 
lat_lon_site_fix[,DepthZone:= factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","ARM"),
                                     labels = c("Inner (5m)","Middle (10m)","Outer (15m)","Deep (25m)","AR"))]

#avg values
#mean lat and lon
lat_lon_site_fix[,avg_lon := mean(Longitude,na.rm = T),Site][,avg_lat := mean(Latitude,na.rm = T),Site]
lat_lon_site_fix.r <- unique(lat_lon_site_fix[,.(Site, `Reef type`, avg_lon, avg_lat)])

#Project points to simple feature
lat_lon_site_fix.sf <- st_as_sf(lat_lon_site_fix, coords = c("Longitude","Latitude"), crs = st_crs(4326)) #all depth zones
lat_lon_site_fix.r.sf <- st_as_sf(lat_lon_site_fix.r, coords = c("avg_lon","avg_lat"), crs = st_crs(4326)) #avg across depth zones for each site

##############################################type_sum()
#Map
#######################
#set square from which to extract bathy data from NOAA server
bathy_VRG <- getNOAA.bathy(min(lat_lon_site_fix.r$avg_lon)-2, max(lat_lon_site_fix.r$avg_lon)+2, min(lat_lon_site_fix.r$avg_lat)-0.5, max(lat_lon_site_fix.r$avg_lat)+0.5, resolution = 0.000001) #bathymetry matrix
bathy_VRG_xyz <- data.table(as.xyz(bathy_VRG))
#only 200m
bathy_VRG_xyz.200m <- bathy_VRG_xyz[V3 >= -200,]

#california
country <- ne_countries(scale = "small",returnclass = "sf")

#map of bathymetry
site_map_colored <- autoplot.bathy(bathy_VRG, geom=c("raster"
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

ggsave(site_map_colored, filename = "site_map_colored.jpg", path = file.path("figures"), width = 10, height = 7, units = "in")

#site map b&W
site_map_colored <- autoplot.bathy(bathy_VRG, geom=c("raster"
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

ggsave(site_map_colored, filename = "site_map_colored.jpg", path = file.path("figures"), width = 10, height = 7, units = "in")




# Plot using ggplot and sf
# Get US states and Mexico
usa <- rnaturalearth::ne_countries(country = "United States of America", returnclass = "sf", scale = "large")
mexico <- rnaturalearth::ne_states(country="Mexico", returnclass = "sf")

# Load the state boundaries for the USA
states <- ne_states(country = "united states of america", returnclass = "sf")

# Filter to get only California
california <- states[states$name == "California", ]

#Higher rez California
CA_Map <- sf::st_read("~/Dropbox/VRG Files/R Code/Mapping/CA_Map_Nov2023.shp")

#plot
site_map_basic <- ggplot() + 
 # geom_tile(data = bathy_VRG_xyz.200m, aes(x = V1, y = V2), fill = "lightblue", color = "lightblue") +
  geom_sf(data = usa) +
  geom_sf(data = mexico) +
  geom_point(data = lat_lon_site_fix.r, aes(x = avg_lon, y = avg_lat, fill = `Reef type`, shape = `Reef type`), size = 2, color = "black") +
  scale_fill_manual(values = c("darkturquoise","brown1"))+
  scale_shape_manual(values = c(21,24)) +
  coord_sf(xlim = c(-120.5,-116.75), ylim = c(32.6, 34.1), expand = F) +
  geom_rect(aes(xmin = -118.27,
                xmax = -118.48,
                ymin = 33.68,
                ymax = 33.8),
            fill = NA, color = "black", size = 0.5, linetype = "dashed") +  # Add red rectangle to highlight PVR
  theme_classic() +
  theme(axis.title = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(0.72,0.22), legend.direction = "vertical",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12)) +
  guides(
    fill = guide_legend(override.aes = list(size = 4)),
    shape = guide_legend(override.aes = list(size = 4))
  )

ggsave(site_map_basic, filename = "site_map_basic.jpg", path = file.path("figures"), width = 15, height = 10, units = "in", dpi = 300)

# Create an inset plot for California
california_inset <- ggplot() +
  geom_sf(data = california, lwd = 0.3) +
  geom_rect(aes(xmin = -120.5,
            xmax = -116.75,
            ymin = 32.6,
            ymax = 34.3),
            fill = NA, color = "black", size = 0.4) +  # Add red rectangle to highlight the area
  theme_void()

#Convert lat/lon for each site to match CA map crs
lat_lon_site_fix.sf <- st_transform(lat_lon_site_fix.sf, crs = st_crs(CA_Map))

#Extract legend
inset_leg <- get_legend(ggplot() +
  geom_point(data = lat_lon_site_fix, aes(x = Longitude, y = Latitude,color = DepthZone, shape = DepthZone), size = 2) +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100","black")) +
    scale_shape_manual(values = c(17,17,17,17,21)) +
    labs(color = "Depth zone/\nreef type", shape = "Depth zone/\nreef type") +
    theme_classic() +
    theme(legend.background = element_rect(fill = "NA")))

#Create inset plot for PVR
PVR_inset <- ggplot() +
  geom_sf(data = CA_Map) +
  geom_sf(data = lat_lon_site_fix.sf, aes(fill = DepthZone, shape = `Reef type`), size = 1.5, color = "black", stroke = 0.1) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100","white"))+
  scale_shape_manual(values = c(21,24)) +
 coord_sf(xlim = c(366453,381074.8), ylim = c(3729596, 3741414), expand = F) + #bounding box values from "PalosVerdes_DepthZones.shp"
  ggspatial::annotation_scale(location = 'bl') +
  theme_classic() +
  theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1, linetype = "dashed"), legend.position = "null", plot.background = element_rect(fill = NA))


# Add the CA outline inset to the main plot
site_map_with_inset <- site_map_basic +
  annotation_custom(
    grob = ggplotGrob(california_inset), 
    xmin = -117.7, 
    ymin = 33.1)

#Add legend, and small PVR plot to main plot
site_map_with_insets <- ggdraw(site_map_with_inset) +
  draw_plot(PVR_inset, x = -0.28, y = 0.059, height = 0.55) +
  draw_plot(inset_leg, x = -0.19, y = 0.36, height = 0.05) +
  geom_segment(aes(x = 0.55, y = 0.73, xend = 0.37, yend = 0.64),
               arrow = arrow()) +
  geom_text(aes(label = "Palos Verdes Peninsula",x = 0.22, y = 0.63), size = 5) +
  geom_text(aes(label = "Santa Monica Bay"), x = 0.46, y = 0.85, size = 5)

# Save the plot
ggsave(site_map_with_insets, filename = "site_map_basic_with_insets.jpg", path = file.path("figures"), width = 9, height =4.5, units = "in", dpi = 300)

#Identify sites overlapping with MPAs
#Bring in MPA polygons
CA_mpas <- st_read(file.path("data","California_Marine_Protected_Areas_shp","California_Marine_Protected_Areas.shp"))

#Adjust CRS
CA_mpas.t <- st_transform(CA_mpas, crs = st_crs(lat_lon_site_fix.r.sf))

#Find overlap
overlaps <- st_intersects(lat_lon_site_fix.r.sf, CA_mpas.t)

#Add column to df identifying whether it overlaps with MPA (lengths = length of list or vector elements)
lat_lon_site_fix.r.sf$MPA_overlap <- lengths(overlaps) > 0

#Extract just site name depth zone and MPA overlap status
MPA_site_key <- data.table(lat_lon_site_fix.r.sf)

MPA_site_key <- unique(MPA_site_key[,.(Site,MPA_overlap)])

#Save to keys folder
saveRDS(MPA_site_key, file.path("keys","MPA_site_key.rds"))


