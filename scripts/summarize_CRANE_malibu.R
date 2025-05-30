# CREATION DATE 7 July 2024
# MODIFIED DATE 19 December 2024

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Abundance and diversity with depth

#############################
##Setup
#############################
library(ggplot2)
library(data.table)
library(cowplot)
library(ggpattern) #patterned bars
library(RColorBrewer)
library(nlme)
library(MuMIn)
library(ggpattern)
library(multcompView)
library(grid)
library(MASS)
library(patchwork)
library(dunn.test)
#library(rcompanion)

########################
##Load data
########################
dat_event_malibu <- readRDS(file.path("data","processed_crane", "dat_event_malibu.rds"))
dat_fish_site_averages_malibu <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages_malibu.rds"))
dat_macroinvert_site_averages_malibu <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages_malibu.rds"))
dat_kelp_site_averages_malibu <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages_malibu.rds"))

#included giant kelp stipes for in situ habitat data, but DELETE for community analyses
dat_kelp_site_averages_malibu <- dat_kelp_site_averages_malibu[BenthicReefSpecies != "Macrocystis pyrifera stipes",] 

#Just two sites that overlap fire area and watershed area, and only since 2016
dat_event_malibu.r <- dat_event_malibu[Site %in% c("Big Rock","Malibu Bluffs")][SampleYear > 2016]
dat_fish_site_averages_malibu.r <- dat_fish_site_averages_malibu[Site %in% c("Big Rock","Malibu Bluffs")]
dat_macroinvert_site_averages_malibu.r <- dat_macroinvert_site_averages_malibu[Site %in% c("Big Rock","Malibu Bluffs")]
dat_kelp_site_averages_malibu.r <- dat_kelp_site_averages_malibu[Site %in% c("Big Rock","Malibu Bluffs")]

#pull in spp taxonomy info
species_key <- fread(file.path("keys","species_key.csv"))
#capitalize California
species_key[, common_name_final := gsub("california sheephead", "California sheephead", common_name_final)]

#find and replace Semicossyphus pulcher with Bodianus pulcher in VRG data
dat_fish_site_averages_malibu[, Species := gsub("Semicossyphus pulcher", "Bodianus pulcher", Species)]

#link site averaged data with species key
dat_fish_site_averages_malibu.r <- species_key[dat_fish_site_averages_malibu.r, on = c("taxa" = "Species")]
dat_macroinvert_site_averages_malibu.r <- species_key[dat_macroinvert_site_averages_malibu.r, on = c("taxa" = "BenthicReefSpecies")]
dat_kelp_site_averages_malibu.r <- species_key[dat_kelp_site_averages_malibu.r, on = c("taxa" = "BenthicReefSpecies")]

########################
##Total abundance (summed density of all taxa, summed biomass of all taxa) by depth zone
########################

#abundance for fish
dat_fish_site_averages_malibu.r[,total_abundance_depthzone_site := sum(mean_density_m2)*100,.(Site, DepthZone)] #multiply by 100 to get # per 100m^2
dat_fish_site_averages_malibu.r[,total_biomass_depthzone_site := sum(mean_wt_density_g_m2)*100/1000,.(Site, DepthZone)] #total fish biomass in kg per 100m^2 (multiply by 100, divide by 1000)

dat_fish_total_abundances <- unique(dat_fish_site_averages_malibu.r[,.(Site, DepthZone, total_abundance_depthzone_site, total_biomass_depthzone_site)])

#Set factor order for depth zone
dat_fish_total_abundances[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep"))]

#abundance for macroinverts
dat_macroinvert_site_averages_malibu.r[,total_abundance_depthzone_site := sum(mean_density_m2)*100,.(Site, DepthZone)] #multiply by 100 to get # per 100m^2
dat_macroinvert_total_abundances <- unique(dat_macroinvert_site_averages_malibu.r[,.(Site, DepthZone, total_abundance_depthzone_site)])

#Set factor order for depth zone
dat_macroinvert_total_abundances[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep"))]

#abundance for kelp
dat_kelp_site_averages_malibu.r[,total_abundance_depthzone_site := sum(mean_density_m2)*100,.(Site, DepthZone)] #multiply by 100 to get # per 100m^2
dat_kelp_total_abundances <- unique(dat_kelp_site_averages_malibu.r[,.(Site, DepthZone,  total_abundance_depthzone_site)])
#Set factor order for depth zone
dat_kelp_total_abundances[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep"))]


#Visualize all fish density ####
fish_abundance_depthzone_malibu <- ggplot(dat_fish_total_abundances) +
  geom_boxplot(aes(x = Site, y = total_abundance_depthzone_site, color = DepthZone), position = position_dodge2(preserve = "single")) +
  labs(x = "", y = bquote("   Density\n(count per 100 m"^2*")") , color = "Depth") +
  theme_classic() +
  ggtitle("Fish")+
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave(fish_abundance_depthzone_malibu, path = file.path("figures"), filename ="fish_abundance_depthzone_malibu.jpg", height = 4.5, width = 6, units = "in")

#Tables of species sorted by abundance
malibu_fish_table <- dat_fish_site_averages_malibu.r[mean_density_m2>0,.(taxa,common_name_final, Site, DepthZone, mean_density_m2, mean_wt_density_g_m2)]

#Sort by Site and then depth zone, and then descending by frequency
#Just do in excel
fwrite(malibu_fish_table, file = file.path("output","malibu_fish_table.csv"))

#Tables of species sorted by abundance
malibu_macroinvert_table <- dat_macroinvert_site_averages_malibu.r[mean_density_m2>0,.(taxa,common_name_final, Site, DepthZone, mean_density_m2)]

#Sort by Site and then depth zone, and then descending by frequency
#Just do in excel
fwrite(malibu_macroinvert_table, file = file.path("output","malibu_macroinvert_table.csv"))

#Tables of species sorted by abundance
malibu_kelp_table <- dat_kelp_site_averages_malibu.r[mean_density_m2>0,.(taxa,common_name_final, Site, DepthZone, mean_density_m2)]

#Sort by Site and then depth zone, and then descending by frequency
#Just do in excel
fwrite(malibu_kelp_table, file = file.path("output","malibu_kelp_table.csv"))
