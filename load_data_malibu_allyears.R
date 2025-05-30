# CREATION DATE 24 Jan 2025
# MODIFIED DATE 24 Jan 2025

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Pull in Malibu data for all years for Wildfire~rocky reef proposal

#############################
##Setup
#############################

library(data.table)
library(vegan)
library(dplyr)

#############################
##Load data
#############################

# Read in reference tables (conversions, station names, species names, codes for UPC, ISC, etc.)
source("~/Dropbox/VRG Files/R Code/DataFiles/READ_Sites_Fish_BRS.r")

#Loads up and processes CRANE data
source("~/Dropbox/VRG Files/R Code/Integrated Dive General/CRANE_data_prep.R")

CRANE_data_prep(SELECT = "Malibu_sites",
                add_0s = "sp_0s", #summed to by Species with 0s added
                fish_prod = F, #check what this does, don't think I need it now
                Rem_YOYs_dens = T) #talk to Jeremy about this, see FISH_PREP2.R for details on what this does

########################
##Split into fish, inverts, and kelp
########################
#EVENT data
#dat_event

#FISH (both density (density_m2) and biomass (wt_density_g_m2))
#summed across rows with same Transect & Species, zeros added when not present!
#dat_fish_t

#MACROINVERTS (SWATH)
#swath_melt_T_sp, and then %in% "nudibranchs", "anemones", "sea cucumbers", "sea urchins",
#"sea stars", "molluscs - mobile", "crustaceans", "sponges", "molluscs - sessile", "gorgonians",
#"tunicates",  "hydrocorals", "sea pens" 
dat_macroinvert <- swath_melt_T_sp %>%
  filter(SpeciesGroupF %in% c("nudibranchs", "anemones", "sea cucumbers", "sea urchins", "sea stars", "molluscs - mobile",
                              "crustaceans", "sponges", "molluscs - sessile", "gorgonians", "tunicates",  "hydrocorals", "sea pens"))

#KELP (SWATH)
#swath_melt_T_sp, and then %in% "kelp - understory", "kelp - canopy" ,"Sargassum", "green algae", "giant kelp stipes"
#      "giant kelp stipes" included
dat_kelp <- swath_melt_T_sp %>%
  filter(SpeciesGroupF %in% c("kelp - understory", "kelp - canopy" ,"Sargassum", "green algae","giant kelp stipes"))

########################
##Add Level Orders to Region, Site, DepthZone (mostly from Jeremy's code)
########################
source("~/Dropbox/VRG Files/R Code/Integrated Dive General/CRANE_level_orders.R")

dat_fish_t <- CRANE_level_orders(
  dat_fish_t,
  PVR_Monitor_Cat = T, #need this so PVR_Reefing_Area works below
  AR_Region = T,
  AR_Complex = T,
  # reassign relevant SiteType & DepthZone (Halo) for Reefing Area
  PVR_Reefing_Area = T #after this can filter out DepthZone != "Halo"
) |> 
  filter(DepthZone != "Halo") |> #remove halo observations
  droplevels()

dat_macroinvert <- CRANE_level_orders(dat_macroinvert, AR_Region = T, AR_Complex = T) |> 
  droplevels()

dat_kelp <- CRANE_level_orders(dat_kelp, AR_Region = T, AR_Complex = T) |> 
  droplevels()

########################
##Filtering (mostly from Jeremy's code)
########################

# Filter out PVR Pre-Construction & sampling in 2020 right after construction
dat_fish_t <- dat_fish_t |>
  filter(!(str_detect(Site, "PVR") & SampleYear < 2021 )) |> 
  droplevels()

dat_fish_t <- dat_fish_t |>
  filter(!Site %in% c("Leucadia", "Old 18th")) |> #some depth zones were way out there in nMDS (probably very few species?)
  droplevels() 

dat_macroinvert <- dat_macroinvert |>
  filter(!(str_detect(Site, "PVR") & SampleYear < 2021 )) |> 
  droplevels()

dat_UPC_VRG <- dat_UPC_VRG |>
  filter(!(str_detect(Site, "PVR") & SampleYear < 2021 )) |> 
  droplevels()

dat_macroinvert <- dat_macroinvert |>
  filter(!Site %in% c("Leucadia", "Old 18th")) |> #some depth zones were way out there in nMDS (probably very few species?)
  droplevels() 

dat_kelp <- dat_kelp |>
  filter(!(str_detect(Site, "PVR") & SampleYear < 2021 )) |> 
  droplevels()

dat_kelp <- dat_kelp |>
  filter(!Site %in% c("Leucadia", "Old 18th")) |> #some depth zones were way out there in nMDS (probably very few species?)
  droplevels() 

dat_UPC_VRG <- dat_UPC_VRG |>
  filter(!Site %in% c("Leucadia", "Old 18th")) |> #some depth zones were way out there in nMDS (probably very few species?)
  droplevels() 


########################
##Restrict to well sampled sites for averages (sampled in at least X years)
########################

#Switch to data table

dat_event <- data.table(dat_event)
dat_fish_t <- data.table(dat_fish_t)
dat_macroinvert <- data.table(dat_macroinvert)
dat_kelp <- data.table(dat_kelp)

  #Malibu only
  
  #reduce all data tables to sites that are in Malibut, or are above 33.8977 lat (above manhattan beach)
  dat_event.malibu <- dat_event[Region == "Malibu" | Latitude > 33.8977362]
  
  #Site names
  sites_keep <- unique(dat_event.malibu$Site)
  
  dat_fish_t.malibu <- dat_fish_t[Site %in% sites_keep,]
  dat_macroinvert.malibu <- dat_macroinvert[Site %in% sites_keep,]
  dat_kelp.malibu <- dat_kelp[Site %in% sites_keep,]
  
  ########################
  ##Remove TBF 2016 Project (ASK WHAT THIS IS AND IF I SHOULD LEAVE IT IN)
  ########################
  dat_event.malibu <- dat_event.malibu[Project != "TBF2016",]
  dat_fish_t.malibu <- dat_fish_t.malibu[Project != "TBF2016",]
  dat_macroinvert.malibu <- dat_macroinvert.malibu[Project != "TBF2016",]
  dat_kelp.malibu <- dat_kelp.malibu[Project != "TBF2016",]
  
  ########################
  ##Do not keep early years, not totally trust worthy
  ########################
  dat_event.malibu <- dat_event.malibu[Project != "TBF2016",][SampleYear > 2010]
  dat_fish_t.malibu <- dat_fish_t.malibu[Project != "TBF2016",][SampleYear > 2010]
  dat_macroinvert.malibu <- dat_macroinvert.malibu[Project != "TBF2016",][SampleYear > 2010]
  dat_kelp.malibu <- dat_kelp.malibu[Project != "TBF2016",][SampleYear > 2010]
  
  ########################
  ##Save full UPC output (so we can calculate UPC relief/substrate metrics)
  ########################
  UPC_complete_malibu <- UPC_complete
  
  saveRDS(UPC_complete_malibu, file.path("data","full_crane","UPC_complete_malibu.rds"))
  
  ########################
  ##Take average values across all years of sampling (avg density and avg biomass of each species at site)
  ########################
  
  dat_fish_site_averages_malibu <- dat_fish_t.malibu[,.(mean_density_m2=mean(density_m2), mean_wt_density_g_m2=mean(wt_density_g_m2)),
                                         .(Species, Region, AR_Complex, Site, DepthZone)]
  dat_macroinvert_site_averages_malibu <- dat_macroinvert.malibu[,.(mean_density_m2=mean(Abundance/area.m2)),
                                                     .(BenthicReefSpecies, SpeciesGroupF, Region, AR_Complex, Site, DepthZone)]
  dat_kelp_site_averages_malibu <- dat_kelp.malibu[,.(mean_density_m2=mean(Abundance/area.m2)),
                                       .(BenthicReefSpecies, SpeciesGroupF, Region, AR_Complex, Site, DepthZone)]
  
  #Which regions are retained?
  unique(dat_event.malibu$Region)
  
  #What years do we have data from?
  unique(dat_event.malibu$SampleYear) #2008, 2004, 2007, 2009, 2019, 2020, 2021, 2022, 2023
  
  ########################
  ##Save average output
  ########################
  saveRDS(dat_event.malibu, file.path("data","processed_crane", "dat_event_malibu.rds"))
  saveRDS(dat_fish_site_averages_malibu, file.path("data","processed_crane", "dat_fish_site_averages_malibu.rds"))
  saveRDS(dat_macroinvert_site_averages_malibu, file.path("data","processed_crane", "dat_macroinvert_site_averages_malibu.rds"))
  saveRDS(dat_kelp_site_averages_malibu, file.path("data","processed_crane", "dat_kelp_site_averages_malibu.rds"))


#Identify inside outside MPA
MPA_site_key <- readRDS(file.path("keys","MPA_site_key.rds"))
dat_event.malibu <- MPA_site_key[dat_event.malibu, on = "Site"]

#Add lat lon from dive priority list
lat_lon_site_fix.r <- fread(file.path("keys","lat_lon_site_fix.r.csv"))
dat_event.malibu <- lat_lon_site_fix.r[dat_event.malibu, on = "Site"]


