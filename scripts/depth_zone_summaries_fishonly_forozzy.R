# CREATION DATE 30 Jul 2024
# MODIFIED DATE 30 Jul 2024

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Pull fish data for Ozzy to build models (most common and hefty fish)

#############################
##Setup
#############################
library(ggplot2)
library(data.table)
library(cowplot)
library(ggpattern) #patterned bars
library(RColorBrewer)

########################
##Load data
########################
dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))
dat_fish_site_averages <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages.rds"))

#pull in spp taxonomy info
species_key <- fread(file.path("keys","species_key.csv"))
#capitalize California
species_key[, common_name_final := gsub("california sheephead", "California sheephead", common_name_final)]

#find and replace Semicossyphus pulcher with Bodianus pulcher in VRG data
dat_fish_site_averages[, Species := gsub("Semicossyphus pulcher", "Bodianus pulcher", Species)]

#find and replace Hermosilla azurea with Kyphosus azureus in VRG data
dat_fish_site_averages[, Species := gsub("Hermosilla azurea","Kyphosus azureus", Species)]

#link site averaged data with species key
dat_fish_site_averages <- species_key[dat_fish_site_averages, on = c("taxa" = "Species")]

#number of sites per depth zone
dat_event.r[,number_sites_depthzone := uniqueN(Site),.(DepthZone)]

number_sites_depthzone <- unique(dat_event.r[,.(DepthZone, number_sites_depthzone)])

########################
##Averaged across all sites, top species per depth zone
########################

#fish
dat_fish_site_averages[,mean_depthzone_density_m2 := mean(mean_density_m2),.(taxa, DepthZone)] 
dat_fish_site_averages[,mean_depthzone_wt_density_g_m2 := mean(mean_wt_density_g_m2),.(taxa, DepthZone)] 

col_keep_fish <- colnames(dat_fish_site_averages[,c(1:10,14,17,18)])

dat_fish_averages <- unique(dat_fish_site_averages[,..col_keep_fish])

#Add column for whether species is in top 10 species (by density) for that zone
#fish density
dat_fish_averages[, Species_top10 := ifelse(frank(-mean_depthzone_density_m2)<=10,taxa,"Other"), .(DepthZone)]
dat_fish_averages[, common_top10 := ifelse(frank(-mean_depthzone_density_m2)<=10,common_name_final,""), .(DepthZone)]
dat_fish_averages[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top10)]
dat_fishdensity_averages_deep_ar.u <- unique(dat_fish_averages[,.(Species_top10,common_top10, DepthZone, summed_mean_depthzone_density_m2)])
dat_fishdensity_averages_deep_ar.u[,full_label := ifelse(Species_top10 == "Other","Other",paste0(Species_top10,"\n", common_top10))]

#top species by biomass all sites
top_spp_density <- unique(dat_fishdensity_averages_deep_ar.u[Species_top10 != "Other",.(Species_top10,common_top10)])

fwrite(top_spp_density, file.path("output","top_spp_lists","top_spp_density.csv"))

#Add column for whether species is in top 5 species (by biomass) for that zone, and then sum biomass for all others
#fish biomass
dat_fish_averages[, Species_top10_biomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=10,taxa,"Other"), .(DepthZone)]
dat_fish_averages[, common_top10_biomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=10,common_name_final,""), .(DepthZone)]
dat_fish_averages[, summed_mean_depthzone_biomass_m2 := sum(mean_depthzone_wt_density_g_m2), .(DepthZone, Species_top10_biomass)]
dat_fishbiomass_averages_deep_ar.u <- unique(dat_fish_averages[,.(Species_top10_biomass,common_top10_biomass, DepthZone, summed_mean_depthzone_biomass_m2)])
dat_fishbiomass_averages_deep_ar.u[,full_label := ifelse(Species_top10_biomass == "Other","Other",paste0(Species_top10_biomass,"\n", common_top10_biomass))]

#top species by biomass all sites
top_spp_biomass <- unique(dat_fishbiomass_averages_deep_ar.u[Species_top10_biomass != "Other",.(Species_top10_biomass,common_top10_biomass)])

fwrite(top_spp_biomass, file.path("output","top_spp_lists","top_spp_biomass.csv"))

#merge top biomass and top density

top_spp <- unique(rbind(top_spp_biomass, top_spp_density, use.names = F)[,.(Species_top10_biomass,common_top10_biomass)])

colnames(top_spp) <- c("spp","common")

fwrite(top_spp, file.path("output","top_spp_lists","top_spp.csv"))

############################################################
#Now, average biomass and density by Island vs ARM vs Coast
############################################################

#New Column Identifying ARM vs Island vs Natural Coast
dat_fish_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]

#average by type and depth zone
dat_fish_averages_sitetype <- dat_fish_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2),
                                                        mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                                     .(taxa, common_name_final, DepthZone, type)] 

#Identify top 5 species, sum other into 'other' category
###FISH DENSITY#####

dat_fish_averages_sitetype[, Species_top10 := ifelse(frank(-mean_depthzone_density_m2)<=10,taxa,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, common_top10 := ifelse(frank(-mean_depthzone_density_m2)<=10,common_name_final,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, summed_mean_depthzone_sitetype_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, type, Species_top10, common_top10)]
dat_fishdensity_averages_sitetype.u <- unique(dat_fish_averages_sitetype[,.(Species_top10, type, DepthZone, summed_mean_depthzone_sitetype_density_m2, common_top10)])
dat_fishdensity_averages_sitetype.u[,full_label := ifelse(Species_top10 == "Other","Other",paste0(Species_top10,"\n", common_top10))]


fish_spp_density_type_island_unique <- unique(dat_fishdensity_averages_sitetype.u[Species_top10 != "Other" & type == "Island",.(Species_top10, common_top10)])
fwrite(fish_spp_density_type_island_unique, file.path("output","top_spp_lists","fish_spp_density_type_island_unique.csv"))
fish_spp_density_type_mainland_unique <- unique(dat_fishdensity_averages_sitetype.u[Species_top10 != "Other" &type == "Natural mainland",.(Species_top10, common_top10)])
fwrite(fish_spp_density_type_mainland_unique, file.path("output","top_spp_lists","fish_spp_density_type_mainland_unique.csv"))
fish_spp_density_type_AR_unique <- unique(dat_fishdensity_averages_sitetype.u[Species_top10 != "Other"&type == "ARM",.(Species_top10, common_top10)])
fwrite(fish_spp_density_type_AR_unique, file.path("output","top_spp_lists","fish_spp_density_type_AR_unique.csv"))
fish_spp_density_type_all_unique <- unique(dat_fishdensity_averages_sitetype.u[Species_top10 != "Other",.(Species_top10, common_top10)])
fwrite(fish_spp_density_type_all_unique, file.path("output","top_spp_lists","fish_spp_density_type_all_unique.csv"))


######FISHBIOMASS######
dat_fish_averages_sitetype[, Species_top10_fishbiomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,taxa,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, common_top10_fishbiomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,common_name_final,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, summed_mean_depthzone_sitetype_biomass_m2 := sum(mean_depthzone_wt_density_g_m2), .(DepthZone, type, Species_top10_fishbiomass, common_top10_fishbiomass)]
dat_fishbiomass_averages_sitetype.u <- unique(dat_fish_averages_sitetype[,.(Species_top10_fishbiomass, type, DepthZone, summed_mean_depthzone_sitetype_biomass_m2, common_top10_fishbiomass)])
dat_fishbiomass_averages_sitetype.u[,full_label := ifelse(Species_top10_fishbiomass == "Other","Other",paste0(Species_top10_fishbiomass,"\n", common_top10_fishbiomass))]


fish_spp_biomass_type_island_unique <- unique(dat_fishbiomass_averages_sitetype.u[Species_top10_fishbiomass != "Other"&type == "Island",.(Species_top10_fishbiomass, common_top10_fishbiomass)])
fwrite(fish_spp_biomass_type_island_unique, file.path("output","top_spp_lists","fish_spp_biomass_type_island_unique.csv"))
fish_spp_biomass_type_mainland_unique <- unique(dat_fishbiomass_averages_sitetype.u[Species_top10_fishbiomass != "Other"&type == "Natural mainland",.(Species_top10_fishbiomass, common_top10_fishbiomass)])
fwrite(fish_spp_biomass_type_mainland_unique, file.path("output","top_spp_lists","fish_spp_biomass_type_mainland_unique.csv"))
fish_spp_biomass_type_AR_unique <- unique(dat_fishbiomass_averages_sitetype.u[Species_top10_fishbiomass != "Other"&type == "ARM",.(Species_top10_fishbiomass, common_top10_fishbiomass)])
fwrite(fish_spp_biomass_type_AR_unique, file.path("output","top_spp_lists","fish_spp_biomass_type_AR_unique.csv"))
fish_spp_biomass_type_all_unique <- unique(dat_fishbiomass_averages_sitetype.u[Species_top10_fishbiomass != "Other",.(Species_top10_fishbiomass, common_top10_fishbiomass)])
fwrite(fish_spp_biomass_type_all_unique, file.path("output","top_spp_lists","fish_spp_biomass_type_all_unique.csv"))

