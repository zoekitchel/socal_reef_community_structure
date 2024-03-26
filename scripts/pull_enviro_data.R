# CREATION DATE 14 Jan 2023

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

#############
#Prep in situ habitat characteristics
#############

#######KELP
#first, only macrocystis
dat_kelp_site_averages <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages.rds"))

#kelp
dat_kelp_long_density <- dat_kelp_site_averages[,.(BenthicReefSpecies,Region,Site,DepthZone, mean_density_m2)]
dat_kelp_long_density[,mean_density_m2 := mean(mean_density_m2),.(BenthicReefSpecies,Region,Site,DepthZone)]
dat_kelp_long_density <- unique(dat_kelp_long_density[,.(BenthicReefSpecies,Region,Site,DepthZone,mean_density_m2)])

macro_density_bysite <- dat_kelp_long_density[BenthicReefSpecies == "Macrocystis pyrifera"][,.(Site, DepthZone, mean_density_m2)]

colnames(macro_density_bysite) <- c("Site","DepthZone","macro_mean_density_m2")

#save
saveRDS(macro_density_bysite, file.path("data","enviro_predictors","macro_density_bysite.rds"))

###########relief


###########substrate


##########
#Prep large scale habitat characteristics
##########
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

#######Distance to 200 m isobath

dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))

df.with_depth <- add_depth_columns(dat_event.r, ETOPO = F, CDFW = F, USGS_socal=F, dist_200m = T)

df.with_depth[,dist_200m_bath := mean(dist_200m_bath),.(Site,DepthZone)]

distance_200mbathy_bysite <-  unique(df.with_depth[,.(Site, DepthZone, dist_200m_bath)])


#save
saveRDS(distance_200mbathy_bysite, file.path("data","enviro_predictors","distance_200mbathy_bysite.rds"))
