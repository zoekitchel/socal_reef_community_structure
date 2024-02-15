# CREATION DATE 14 Jan 2023

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Pull environmental data, both in situ and external

#############################
##Setup
#############################
library(ggplot2)
library(data.table)
library(dplyr)

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

#redice
