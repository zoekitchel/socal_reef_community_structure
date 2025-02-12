# CREATION DATE 29 Jan 2025
# MODIFIED DATE 3 Feb 2025

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Map of Malibu with fire inset for 2 pager

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
# Get correct lat lon values for Malibut sites ####
#############################

#Chelsea recommends using Lat/Lon from "2023 Dive Site Priority" instead of All Sites File
dive_site_priority_list <- fread("dive_site_priority_list.csv")

#extract depth zone and site from location column
dive_site_priority_list[,DepthZone := word(Location,-1)][,Site := word(Location, start = 1, end = -2)]

#Additionally manual fixes to wonky names without depthzones
dive_site_priority_list[Location == "Malibu Bluffs Eelgrass",DepthZone := "Middle"][Location == "Malibu Bluffs Eelgrass",Site := "Malibu Bluffs Eelgrass"]
dive_site_priority_list[Location == "Santa Monica Jetty",DepthZone := "ARM"][Location == "Santa Monica Jetty",Site := "Santa Monica Jetty"]
dive_site_priority_list[Location == "Santa Monica Jetty - Interior",DepthZone := "ARM"][Location == "Santa Monica Jetty - Interior",Site := "Santa Monica Jetty"]
dive_site_priority_list[Location == "Marina del Rey Breakwater - Exterior",DepthZone := "ARM"][Location == "Marina del Rey Breakwater - Exterior",Site := "Marina del Rey Breakwater"]

#For SMB AR, MDR AR, and SoS AR, merge into only 3 to make them visible on the map
dive_site_priority_list[grepl("Santa Monica Bay AR", Site) == T, Site := "Santa Monica Bay AR"]
dive_site_priority_list[grepl("Santa Monica AR", Site) == T, Site := "Santa Monica AR"]
dive_site_priority_list[grepl("Marina del Rey AR", Site) == T, Site := "Marina del Rey AR"]

#change column names
colnames(dive_site_priority_list) <- c("Location","Latitude_fix","Longitude_fix","DepthZone","Site")

#limit to lat, lon, depthzone, site
dive_site_priority_list.r <- dive_site_priority_list[,.(Latitude_fix = mean(Latitude_fix,na.rm = T),
                                                        Longitude_fix = mean(Longitude_fix,na.rm = T)),
                                                   .(Site)]

#Load CRANE dat_event_malibu
dat_event_malibu <- readRDS(file.path("data","processed_crane", "dat_event_malibu.rds"))

#Again, merge SMB and MDR AR sites into 1
dat_event_malibu[grepl("Santa Monica AR", Site) == T, Site := "Santa Monica AR"]
dat_event_malibu[grepl("Santa Monica Bay AR", Site) == T, Site := "Santa Monica Bay AR"]
dat_event_malibu[grepl("Marina del Rey AR", Site) == T, Site := "Marina del Rey AR"]

#Add column for reef type
dat_event_malibu[,reef_type := ifelse(DepthZone == "ARM","Artificial_reef","Natural_reef")]

#Manually adjust for Marina del Rey Breakwater - Exterior
dat_event_malibu[Site == "Marina del Rey Breakwater - Exterior", reef_type := "Artificial_reef"]
dat_event_malibu[Site == "Santa Monica Jetty - Interior",Site := "Santa Monica Jetty"]
dat_event_malibu[Site == "Santa Monica Jetty - Exterior",Site := "Santa Monica Jetty"]
dat_event_malibu[Site == "Marina del Rey Breakwater - Exterior",Site := "Marina del Rey Breakwater"]

#Number of years of sampling at a site
dat_event_malibu[,years_sampled := uniqueN(SampleYear),Site]

#Reduce to unique Site
lat_lon_site_orig.malibu <- unique(dat_event_malibu[,.(Site,reef_type, years_sampled)])

#merge event data with fixed lat and lon from dive site priority list
lat_lon_site_fix.malibu <- lat_lon_site_orig.malibu[dive_site_priority_list.r, on = c("Site")]

#Delete any rows without values, and use this as key for Site, Lat and Long
lat_lon_site_fix.malibu <- lat_lon_site_fix.malibu[complete.cases(lat_lon_site_fix.malibu),]

#only keep unique values
lat_lon_site_fix.malibu <- unique(lat_lon_site_fix.malibu) #24 sites

#delete old lat lon columns from all sites
lat_lon_site_fix.malibu <- lat_lon_site_fix.malibu[,.(Site, reef_type, years_sampled, Latitude_fix, Longitude_fix)]

#Per Chelsea's advice, remove Malibu Bluffs Eelgrass and Santa Monica Bay ARs
lat_lon_site_fix.malibu <- lat_lon_site_fix.malibu[!(Site %in% c("Malibu Bluffs Eelgrass", "Santa Monica Bay AR")),]

#change col names
colnames(lat_lon_site_fix.malibu) <- c("Site","Reef_type","Years_sampled","Latitude","Longitude")

#Add column for labeling on map by longitude
setkey(lat_lon_site_fix.malibu, Reef_type, Longitude)

#Label natural reefs as letters, label artificial reefs as numbers
lat_lon_site_fix.malibu[,label := 
                            c(as.character(seq(1,nrow(lat_lon_site_fix.malibu[Reef_type != "Natural_reef"]),by=1)),
                              LETTERS[1:nrow(lat_lon_site_fix.malibu[Reef_type == "Natural_reef"])])]

#Project points to simple feature points to simple feature
lat_lon_site_fix.malibu.sf <- st_as_sf(lat_lon_site_fix.malibu, coords = c("Longitude","Latitude"), crs = st_crs(4326)) #all depth zones

#############################
# Pull in fire outlines for Woolsey, Franklin, and Palisades ####
#############################
#Shapefiles of Woolsey fire perimeter downloaded January 29, 2025 from https://hub-calfire-forestry.hub.arcgis.com/datasets/CALFIRE-Forestry::california-historical-fire-perimeters/explore?layer=2&location=37.251628%2C-119.269051%2C6.27
#Filtered to FIRENAME = Woolsey
#Filename = California_Fire_Perimeters_(1950%2B)/California_Fire_Perimeters_(1950%2B).shp

woolsey_perimeter <- sf::st_read(file.path("data","California_Fire_Perimeters_(1950%2B)","California_Fire_Perimeters_(1950%2B).shp"))

#transform
woolsey_perimeter.t <- st_transform(woolsey_perimeter, crs = 4326)

#select focal columns to rbind with other fires
woolsey_perimeter.t <- woolsey_perimeter.t |>
  select(FIRE_NAME, geometry) |>
  rename(IncidentNa = FIRE_NAME) |>
  mutate(IncidentNa = "2018 Woolsey")

#Shapefiles of Franklin fire perimeter not available as of today. If I find them, I can add, BUT we can argue this wasn't a super urban fire.

#franklin_perimeters <- 

#Current burn perimeter of Palisades fire downloaded January 29, 2025 from https://occidental.maps.arcgis.com/home/item.html?id=d957997ccee7408287a963600a77f61f
#I had to be logged in, but soon this should just be available through calfire as well
#Filename = USA_Current_Wildfires/Current_Perimeters.shp
current_fire_permimeters <- sf::st_read(file.path("data","USA_Current_Wildfires","Current_Perimeters.shp"))
palisades_perimeter <- current_fire_permimeters |> 
  filter(IncidentNa == "PALISADES")

#remove large full perimeter shapefile
rm(current_fire_permimeters)

palisades_perimeter.t <- st_transform(palisades_perimeter, crs = 4326)

palisades_perimeter.t <- palisades_perimeter.t |>
  select(IncidentNa, geometry) |>
  mutate(IncidentNa = "2025 Palisades")

#Merge into a single sf
fire_perimeters.t <- rbind(woolsey_perimeter.t, palisades_perimeter.t)

##############################################type_sum()
#Map
#######################
#Plot using ggmap (google map) remember to 
#add API here: https://console.cloud.google.com/google/maps-apis/credentials?utm_source=Docs_CreateAPIKey&utm_content=Docs_maps-backend&_gl=1*xzqbpi*_ga*ODI3MDgxMjQ1LjE3MzgyMDI3OTI.*_ga_NRWSTWS78N*MTczODIwMjc5Mi4xLjEuMTczODIwMzE3OS4wLjAuMA..&project=alpine-canto-410820

#Malibu map
malibu_basemap <- get_googlemap("Malibu, CA, USA", zoom = 10, maptype = "satellite")


#Malibu fire map with polygons
malibu_firemap <- ggmap(malibu_basemap) +
  geom_sf(data = fire_perimeters.t,inherit.aes = F,
          aes(color = IncidentNa, fill = IncidentNa), alpha = 0.5) + #inherit aes false allows me to plot polygon on top of ggmap without coordinates
  scale_color_manual(values = c("darkred","red")) +
  scale_fill_manual(values = c("darkred","red")) +
  geom_point(data = lat_lon_site_fix.malibu, aes(x = Longitude, y = Latitude), color = "white") +
  geom_text(data = lat_lon_site_fix.malibu, aes(x = Longitude, y = Latitude,label = label),size = 1.5) +
  coord_sf(xlim = c(-119.1,-118.4), ylim = c(33.9, 34.3), expand = F, crs = 4326) +
  theme_classic() +
  labs(color = "Fire incident",fill = "Fire incident") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(0.16,0.19), legend.direction = "vertical",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12)) +
  guides(
    fill = guide_legend(override.aes = list(size = 4)),
    shape = guide_legend(override.aes = list(size = 4))) +   
  geom_rect(aes(xmin = -119.1,
                      xmax = -118.88,
                      ymin = 33.902,
                      ymax = 34.01),
                  fill = "white", color = "white") +  # Add white rectangle to highlight scalebar
  annotation_scale(location = "bl")

#Malibu PALISADES fire map with polygons
malibu_palisades_firemap <- ggmap(malibu_basemap) +
  geom_sf(data = fire_perimeters.t %>% filter(IncidentNa == "2025 Palisades"),inherit.aes = F,
          aes(color = IncidentNa, fill = IncidentNa), alpha = 0.5) + #inherit aes false allows me to plot polygon on top of ggmap without coordinates
  scale_color_manual(values = c("red")) +
  scale_fill_manual(values = c("red")) +
  geom_point(data = lat_lon_site_fix.malibu, aes(x = Longitude, y = Latitude), color = "white") +
  geom_text(data = lat_lon_site_fix.malibu, aes(x = Longitude, y = Latitude,label = label),size = 1.5) +
  coord_sf(xlim = c(-119.0,-118.4), ylim = c(33.9, 34.25), expand = F, crs = 4326) +
  theme_classic() +
  labs(color = "Fire incident",fill = "Fire incident") +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 12),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = c(0.2,0.87), legend.direction = "vertical",
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12)) +
  guides(
    fill = guide_legend(override.aes = list(size = 4)),
    shape = guide_legend(override.aes = list(size = 4))) +   
 # geom_rect(aes(xmin = -119.1,
 #               xmax = -118.88,
 #               ymin = 33.902,
 #               ymax = 34.01),
 #           fill = "white", color = "white") +  # Add white rectangle to highlight scalebar
  annotation_scale(location = "bl")

# Load the state boundaries for the USA
states <- ne_states(country = "united states of america", returnclass = "sf")

# Filter to get only California
california <- states[states$name == "California", ]

# Create an inset plot for California
california_inset <- ggplot() +
  geom_sf(data = california, lwd = 0.3) +
  geom_rect(aes(xmin = -119.1,
                xmax = -118.4,
                ymin = 33.9,
                ymax = 34.3),
             color = "yellow",fill = "yellow", linewidth = 0.6, alpha = 0.3) +  # Add rectangle to highlight the area
  theme_void()

# Add the CA outline inset to the main plot
malibu_site_map_with_inset <- malibu_firemap +
  annotation_custom(
    grob = ggplotGrob(california_inset), 
    xmin = -118.58, 
    ymin = 34.12)

# Save the plot
ggsave(malibu_site_map_with_inset, filename = "malibu_site_map_with_inset.jpg", path = file.path("figures"), width = 6, height =4.5, units = "in", dpi = 300)

# Add the CA outline inset to the Palisades only map
malibu_site_map_with_inset_palisadesonly <- malibu_palisades_firemap +
  annotation_custom(
    grob = ggplotGrob(california_inset), 
    xmin = -118.58, 
    ymin = 34.12)

# Save the plot
ggsave(malibu_site_map_with_inset_palisadesonly, filename = "malibu_site_map_with_inset_palisadesonly.jpg", path = file.path("figures"), width = 6, height =4.5, units = "in", dpi = 300)

####################
# Make tile plot with years that had observations ####
####################

#Unique site, year data table
malibu_site_year <- unique(dat_event_malibu[,.(Site,SampleYear)])

malibu_site_year[,sampled := T]

#Per Chelsea's advice, remove Malibu Bluffs Eelgrass and Santa Monica Bay ARs
malibu_site_year <- malibu_site_year[!(Site %in% c("Malibu Bluffs Eelgrass", "Santa Monica Bay AR")),]

#expand grid for all possibilities
malibu_site_year_all <- data.table(expand.grid(SampleYear = seq(from = min(malibu_site_year$SampleYear), to = 2025,by = 1), Site = malibu_site_year$Site))

#link all possibilities with true observations
malibu_site_year_all <- malibu_site_year[malibu_site_year_all, on = c("Site","SampleYear")]

#If NA, not sampled
malibu_site_year_all[,sampled := ifelse(is.na(sampled),F,T)]

#Manually, what sites sampled in 2024?
malibu_2024_sites <- c(
  "Deep Hole East" ,
  "Deep Hole West" ,
  "Big Rock",
  "El Matador" ,
  "El Pescador",
  "El Sol" ,
  "Escondido East",
  "Escondido West" ,
  "La Piedra",
  "Lechuza" ,
  "Leo Carrillo",
  "Little Dume East" ,
  "Little Dume West",
  "Malibu Bluffs" ,
  "Nicholas Canyon East",
  "Nicholas Canyon West" ,
  "Point Dume",
  "Santa Monica AR",
  "Marina del Rey AR",
  "Santa Monica Jetty",
  "Marina del Rey Breakwater"
  )

#Sampled in 2024
malibu_site_year_all[Site %in% malibu_2024_sites & SampleYear == 2024,sampled := T]

#Link lat lon with labels
malibu_site_year_all <- lat_lon_site_fix.malibu[malibu_site_year_all, on = "Site"]

malibu_site_year_all[,full_label := paste0(label,": ",Site)]

#Create tileplot
# Create the tile plot
malibu_sampling_years <- ggplot(malibu_site_year_all, aes(x = SampleYear, y = reorder(full_label, -Longitude), fill = factor(sampled))) + 
  geom_tile(color = "grey") +  # Tile plot with white borders
  scale_fill_manual(values = c("white", "turquoise"), name = "Sampled", labels = c("No", "Yes")) + 
  labs(x = "Year", y = "Site") + 
  scale_x_continuous(breaks = seq(min(malibu_site_year_all$SampleYear), 2025, by = 1), expand = c(0,0)) +  # Set breaks for each year
  #geom_vline(xintercept = 2018, color = "darkred", size = 1) +  # Add vertical red lines at 2017 for Woolsey
  geom_vline(xintercept = 2024.5, color = "red", size = 1, linetype = "dashed", alpha = 0.8) +  # Add vertical red lines at 2025 for Palisades
  #annotate("text", x = 2018, y = 5, label = " Woolsey", angle = 90, hjust = -3,vjust = -0.5, size = 5, color = "darkred") +  # Add text rotated at 2017
  annotate("text", x = 2025, y = 5, label = "Palisades Fire", angle = 90, hjust = -0.2,size = 5, color = "red") +  # Add text rotated at 2025
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

ggsave(malibu_sampling_years, filename = "malibu_sampling_years.jpg", path = file.path("figures"), width = 12, height =6, units = "in", dpi = 300)

#Abbreviated version
malibu_sampling_years_2019_2024 <- ggplot(malibu_site_year_all, aes(x = SampleYear, y = reorder(full_label, -Longitude), fill = factor(sampled))) + 
  geom_tile(color = "grey") +  # Tile plot with white borders
  scale_fill_manual(values = c("white", "turquoise"), name = "Sampled", labels = c("No", "Yes")) + 
  labs(x = "Year", y = "Site") + 
  scale_x_continuous(breaks = seq(2019, 2025, by = 1), expand = c(0,0), limits = c(2018.5,2025.5)) +  # Set breaks for each year
  #geom_vline(xintercept = 2018, color = "darkred", size = 1) +  # Add vertical red lines at 2017 for Woolsey
  geom_vline(xintercept = 2024.5, color = "red", size = 1, linetype = "dashed", alpha = 0.8) +  # Add vertical red lines at 2025 for Palisades
  #annotate("text", x = 2018, y = 5, label = " Woolsey", angle = 90, hjust = -3,vjust = -0.5, size = 5, color = "darkred") +  # Add text rotated at 2017
  annotate("text", x = 2025, y = 5, label = "Palisades Fire", angle = 90, hjust = -0.2,size = 5, color = "red") +  # Add text rotated at 2025
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


ggsave(malibu_sampling_years_2019_2024, filename = "malibu_sampling_years_2019_2024.jpg", path = file.path("figures"), width = 6, height =6, units = "in", dpi = 300)

#Plot Panels left and right
merged_malibu <- cowplot::plot_grid(malibu_site_map_with_inset,
                   malibu_sampling_years,
                   labels = c("a.","b."),
                   ncol = 1, rel_widths = c(1,1.5))

