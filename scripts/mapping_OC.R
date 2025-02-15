# CREATION DATE 13 Feb 2025
# MODIFIED DATE 13 Feb 2025

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Map of OC for Jayson Smith with timeline of sampling

#############################
##Setup
#############################
library(ggplot2)
library(ggnewscale)
library(ggrepel)
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
# Get correct lat lon values for OC sites ####
#############################

#Chelsea recommends using Lat/Lon from "2023 Dive Site Priority" instead of All Sites File
dive_site_priority_list <- fread("dive_site_priority_list.csv")

#extract depth zone and site from location column
dive_site_priority_list[,DepthZone := word(Location,-1)][,Site := word(Location, start = 1, end = -2)]

#change column names
colnames(dive_site_priority_list) <- c("Location","Latitude_fix","Longitude_fix","DepthZone","Site")

#limit to lat, lon, depthzone, site
dive_site_priority_list.r <- dive_site_priority_list[,.(Latitude_fix = mean(Latitude_fix,na.rm = T),
                                                        Longitude_fix = mean(Longitude_fix,na.rm = T)),
                                                   .(Site)]

#Load CRANE dat_event_malibu
dat_event_OC <- readRDS(file.path("data","processed_crane", "dat_event_OC.rds"))

#Cabrillo Monument - Zone 1 = Inner (use dat event lat lon for this) 32.67007	-117.24709 (Lat Lon from All Sites file)
#Cabrillo Monument - Zone 2 = Inner (use dat event lat lon for this) 32.66668	-117.24737 (Lat Lon from All Sites file)
  #Compare Cabrillo Monument - Zone 1, Cabrillo Monument - Zone 2, and Cabrillo Monument
  #Okay, so, there is NO lat/lon in the data, so, we will ignore
#No Shaw's Cove in dive priority list, will have to use lat lon from different source
  #No lat lon in dive data, Sheephead, only sampled in 2006, pulled lat lon from AllSites file (33.54347	-117.79940)
#No South Carlsbad in dive priority list, will have to use lat lon from different source (use dat event lat lon for this)
  #Only sampled in 2011 for MPA Baseline 2011, lat lon IS in dive data
#No San Elijo in dive priority list, will have to use lat lon from different source (use dat event lat lon for this)
  #Only sampled in 2011 for MPA Baseline 2011, lat lon IS in dive data
#No NBPL (not sure what this is) in dive priority list, will have to use lat lon from different source (use dat event lat lon for this)
  #Only sampled in NBPL 2018, lat lon in dive data, NBPL: Naval Base Point Loma
#Abalone Point
  #33.55218	-117.81995 (From All Sites Data)

#Reduce to unique Site
lat_lon_site_orig.OC <- unique(dat_event_OC[,.(Site)])

#merge event data with fixed lat and lon from dive site priority list
lat_lon_site_fix.OC <- dive_site_priority_list.r[lat_lon_site_orig.OC, on = c("Site")]

#Manually assign some
#Single side, both Inners: 
  #Cabrillo Monument - Zone 1: 32.67007	-117.24709
  #Cabrillo Monument - Zone 2: 32.66668	-117.24737
#San Elijo FROM Event Lat Lon
#South Carlsbad FROM Event Lat Lon
#NBPL FROM Event Lat Lon
#Abalone Point 33.55218	-117.81995
#Shaw's Cove 33.54347	-117.79940
lat_lon_site_fix.OC[,Site_label := Site]
lat_lon_site_fix.OC[Site == "Cabrillo Monument - Zone 1", Latitude_fix := 32.67007][Site == "Cabrillo Monument - Zone 1", Longitude_fix := -117.24709]
lat_lon_site_fix.OC[Site == "Cabrillo Monument - Zone 2", Latitude_fix := 32.66668][Site == "Cabrillo Monument - Zone 2", Longitude_fix := -117.24737]
#Give both the same site name
lat_lon_site_fix.OC[Site == "Cabrillo Monument - Zone 1", Site_label := "Cabrillo National Monument"]
lat_lon_site_fix.OC[Site == "Cabrillo Monument - Zone 2", Site_label := "Cabrillo National Monument"]
lat_lon_site_fix.OC[Site == "San Elijo",Latitude_fix := mean(dat_event_OC[Site == "San Elijo",Latitude])][Site == "San Elijo",Longitude_fix := mean(dat_event_OC[Site == "San Elijo",Longitude])]
lat_lon_site_fix.OC[Site == "South Carlsbad",Latitude_fix := mean(dat_event_OC[Site == "South Carlsbad",Latitude])][Site == "South Carlsbad",Longitude_fix := mean(dat_event_OC[Site == "South Carlsbad",Longitude])]
lat_lon_site_fix.OC[grepl("NBPL",Site)==T, Latitude_fix := mean(dat_event_OC[grepl("NBPL",Site)==T,Latitude])][grepl("NBPL",Site)==T, Longitude_fix := mean(dat_event_OC[grepl("NBPL",Site)==T,Longitude])]
#Give all same site name
lat_lon_site_fix.OC[grepl("NBPL",Site)==T, Site_label := "NBPL"]
lat_lon_site_fix.OC[Site == "Abalone Point",Latitude_fix := 33.55218][Site == "Abalone Point",Longitude_fix := -117.81995]
lat_lon_site_fix.OC[Site == "Shaw's Cove",Latitude_fix := 33.54347][Site == "Shaw's Cove",Longitude_fix := -117.79940]


#only keep unique values
lat_lon_site_fix.OC[,Latitude_fix := mean(Latitude_fix),.(Site_label)][,Longitude_fix := mean(Longitude_fix),.(Site_label)]

#change col names
colnames(lat_lon_site_fix.OC) <- c("Site","Latitude","Longitude","Site_label")

#Link back to dat_event_OC
dat_event_OC <- lat_lon_site_fix.OC[dat_event_OC, on = "Site"]

#Count depth zones sampled
dat_event_OC[,DepthZonesSampled := uniqueN(DepthZone),Site]


#Number of years of sampling at a site
years_sampled_site_key <- unique(dat_event_OC[,.(Site_label,SampleYear, DepthZonesSampled)])

#Max depth zones sampled
max_depth_zones <- years_sampled_site_key[,.(max_depth_zones_count = as.factor(max(DepthZonesSampled))),.(Site_label)]

#set order of max depth zones to 2,3
max_depth_zones[,max_depth_zones_count := factor(max_depth_zones_count, levels = c("2","3"))]

#Add column for labeling on map by longitude
setkey(lat_lon_site_fix.OC,Latitude)

#Label natural reefs as letters
lat_lon_site_fix.OC[,label := rev(LETTERS[1:nrow(lat_lon_site_fix.OC)])]

#Link up with # depth zones sampled
lat_lon_site_fix.OC <- max_depth_zones[lat_lon_site_fix.OC, on = "Site_label"]

#Project points to simple feature points to simple feature
lat_lon_site_fix.OC.sf <- st_as_sf(lat_lon_site_fix.OC, coords = c("Longitude","Latitude"), crs = st_crs(4326)) #all depth zones merged


##############################################
#Map
#######################
#Plot using ggmap (google map) remember to 
#add API here: https://console.cloud.google.com/google/maps-apis/credentials?utm_source=Docs_CreateAPIKey&utm_content=Docs_maps-backend&_gl=1*xzqbpi*_ga*ODI3MDgxMjQ1LjE3MzgyMDI3OTI.*_ga_NRWSTWS78N*MTczODIwMjc5Mi4xLjEuMTczODIwMzE3OS4wLjAuMA..&project=alpine-canto-410820

#Malibu map
OC_basemap <- get_googlemap(center=c(-117.49421906606435,33.32424647992126), zoom = 10, maptype = "satellite")


#Malibu fire map with polygons
OC_sitemap <- ggmap(OC_basemap) +
  geom_point(data = lat_lon_site_fix.OC, aes(x = Longitude, y = Latitude,color = max_depth_zones_count), size = 2) +
  geom_label_repel(data = lat_lon_site_fix.OC, 
                   aes(x = Longitude, y = Latitude,label = label), 
                   size = 3, label.size = 0.5, label.padding = unit(0.5,"mm")) +
  scale_color_manual(values = c("#FF00FF","#FF8000")) +
  coord_sf(expand = F, crs = 4326) +
  theme_classic() +
  labs(color = "Depth zones\nsampled") +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(0.2,0.25), legend.direction = "vertical",
        legend.title = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 12)) +
   annotation_scale(    location = "bl",   # "bl" (bottom-left) is closest to bottom-middle
                        pad_x = unit(0.3, "npc"),  # Moves scale bar towards the center
                        pad_y = unit(0.2, "cm"))    # Slight padding from the bottom)

# Load the state boundaries for the USA
states <- ne_states(country = "united states of america", returnclass = "sf")

# Filter to get only California
california <- states[states$name == "California", ]

# Create an inset plot for California
california_inset <- ggplot() +
  geom_sf(data = california, lwd = 0.3) +
  geom_rect(aes(xmin = -118,
                xmax = -117,
                ymin = 32.95,
                ymax = 33.7),
             color = "yellow",fill = "yellow", linewidth = 0.6, alpha = 0.3) +  # Add rectangle to highlight the area
  theme_void()

# Add the CA outline inset to the main plot
OC_site_map_with_inset <- OC_sitemap +
  annotation_custom(
    grob = ggplotGrob(california_inset), 
    xmin = -117.4, 
    ymin = 33.3)

# Save the plot
ggsave(OC_site_map_with_inset, filename = "OC_site_map_with_inset.jpg",
       path = file.path("figures"), width = 6, height =6, units = "in", dpi = 300)

####################
# Make tile plot with years that had observations ####
####################

#Unique site, year data table
OC_site_year <- unique(dat_event_OC[,.(Site_label,SampleYear)])

OC_site_year[,sampled := T]

#expand grid for all possibilities
OC_site_year_all <- data.table(expand.grid(SampleYear = seq(from = min(OC_site_year$SampleYear), to = 2024,by = 1), Site_label = OC_site_year$Site_label))

#link all possibilities with true observations
OC_site_year_all <- OC_site_year[OC_site_year_all, on = c("Site_label","SampleYear")]

#If NA, not sampled
OC_site_year_all[,sampled := ifelse(is.na(sampled),F,T)]

#Manually, what sites sampled in 2024?
OC_2024_sites <- c(
  "Buck Gully",
"Crystal Cove",
"Heisler Park",
"Laguna Beach",
"Dana Point",
"San Mateo Kelp",
"Leucadia",
"Swami's",
"Matlahuayl",
"Children's Pool",
"South La Jolla",
"Point Loma",
"Cabrillo National Monument"
  )

#Sampled in 2024
OC_site_year_all[Site_label %in% OC_2024_sites & SampleYear == 2024,sampled := T]

#Link lat lon with labels
OC_site_year_all <- lat_lon_site_fix.OC[OC_site_year_all, on = "Site_label"]

OC_site_year_all[,full_label := paste0(label,": ",Site_label)]

#Create tileplot
# Create the tile plot
OC_sampling_years <- ggplot(OC_site_year_all, aes(x = SampleYear, y = reorder(full_label, -Longitude), fill = factor(sampled))) + 
  geom_tile(color = "grey") +  # Tile plot with white borders
  scale_fill_manual(values = c("white", "turquoise"), name = "Sampled", labels = c("No", "Yes")) + 
  labs(x = "Year", y = "Site") + 
  scale_x_continuous(breaks = seq(min(OC_site_year_all$SampleYear), 2024, by = 1), expand = c(0,0)) +  # Set breaks for each year
  theme_classic() + 
  theme(
    legend.position = "top", 
    legend.direction = "horizontal", 
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 12),   # Increase axis text size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels by 45 degrees
    legend.text = element_text(size = 12),  # Increase legend text size
    legend.title = element_text(size = 14), # Increase legend title size
    plot.title = element_text(size = 16),   # Increase plot title size (if you add one)
    strip.text = element_text(size = 12)    # For facet labels if you're using facets
  )

ggsave(OC_sampling_years, filename = "OC_sampling_years.jpg", path = file.path("figures"), width = 12, height =6, units = "in", dpi = 300)
