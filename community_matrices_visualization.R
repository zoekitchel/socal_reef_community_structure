# CREATION DATE 28 Jan 2024
# UPDATED DATE 17 April 2024

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Create depth zone summaries

#NOTE: SOME SITES SURVEYED MULTIPLE TIMES A YEAR FOR ONE PROJECT, NOT SURE HOW TO DEAL WITH THIS, TAKING MEAN FOR NOW

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
library(ggnewscale)
library(ggrepel)
library(gllvm) #latent variable modeling

source(file.path("functions","return_spptaxonomy_function.R"))

########################
##Load data
########################
#biotic
dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))
dat_fish_site_averages <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages.rds"))
dat_macroinvert_site_averages <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages.rds"))
dat_kelp_site_averages <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages.rds"))

#environmental
all_env_lat_lon <- fread(file.path("data","enviro_predictors","all_env_lat_lon.csv"))


##################################################
#Long data to wide data for vegan analyses
##################################################
#fish
dat_fish_long_density <- dat_fish_site_averages[,.(Species,Region,Site,DepthZone, mean_density_m2)]
dat_fish_long_density[,mean_density_m2 := mean(mean_density_m2),.(Species,Region,Site,DepthZone)]
dat_fish_long_density <- unique(dat_fish_long_density[,.(Species,Region,Site,DepthZone,mean_density_m2)])

dat_fish_long_biomass <- dat_fish_site_averages[,.(Species,Region,Site,DepthZone, mean_wt_density_g_m2)]
dat_fish_long_biomass[,mean_wt_density_g_m2 := mean(mean_wt_density_g_m2),.(Species,Region,Site,DepthZone)]
dat_fish_long_biomass <- unique(dat_fish_long_biomass[,.(Species,Region,Site,DepthZone,mean_wt_density_g_m2)])

#macroinvert
dat_macroinvert_long_density <- dat_macroinvert_site_averages[,.(BenthicReefSpecies,Region,Site,DepthZone, mean_density_m2)]
dat_macroinvert_long_density[,mean_density_m2 := mean(mean_density_m2),.(BenthicReefSpecies,Region,Site,DepthZone)]
dat_macroinvert_long_density <- unique(dat_macroinvert_long_density[,.(BenthicReefSpecies,Region,Site,DepthZone,mean_density_m2)])

#kelp
dat_kelp_long_density <- dat_kelp_site_averages[,.(BenthicReefSpecies,Region,Site,DepthZone, mean_density_m2)]
dat_kelp_long_density[,mean_density_m2 := mean(mean_density_m2),.(BenthicReefSpecies,Region,Site,DepthZone)]
dat_kelp_long_density <- unique(dat_kelp_long_density[,.(BenthicReefSpecies,Region,Site,DepthZone,mean_density_m2)])

#melt long to wide
dat_fish_wide_density <- dcast(dat_fish_long_density, Region + Site + DepthZone ~ Species, value.var = "mean_density_m2", fun = mean)

dat_fish_wide_biomass <- dcast(dat_fish_long_biomass, Region + Site + DepthZone ~ Species, value.var = "mean_wt_density_g_m2", fun = mean)

dat_macroinvert_wide_density <- dcast(dat_macroinvert_long_density, Region + Site + DepthZone ~ BenthicReefSpecies, value.var = "mean_density_m2", fun = mean)

dat_kelp_wide_density <- dcast(dat_kelp_long_density, Region + Site + DepthZone ~ BenthicReefSpecies, value.var = "mean_density_m2", fun = mean)

#spp per depthzone
dat_fish_long_density_removezeros <- dat_fish_long_density[mean_density_m2>0,]
dat_fish_long_biomass_removezeros <- dat_fish_long_biomass[mean_wt_density_g_m2>0,]
dat_macroinvert_long_density_removezeros <- dat_macroinvert_long_density[mean_density_m2>0,]
dat_kelp_long_density_removezeros <- dat_kelp_long_density[mean_density_m2>0,]


#number of species per depth zone?
fish_depthzone_site_richness_abun <- dat_fish_long_density_removezeros[,.(count_spp = .N,sum_abun = sum(mean_density_m2)),.(Site, DepthZone, Region)]
fish_depthzone_site_richness_abun[,category:="fish"]
fish_depthzone_site_biomass <- dat_fish_long_biomass_removezeros[,.(sum_biomass = sum(mean_wt_density_g_m2)),.(Site, DepthZone, Region)]
fish_depthzone_site_richness_abun_biomass <- fish_depthzone_site_richness_abun[fish_depthzone_site_biomass, on = c("Site","DepthZone","Region")]

macroinvert_depthzone_site_richness <- dat_macroinvert_long_density_removezeros[,.(count_spp = .N, sum_abun = sum(mean_density_m2)),.(Site, DepthZone, Region)]
macroinvert_depthzone_site_richness[,category:="macroinvert"]
kelp_depthzone_site_richness <- dat_kelp_long_density_removezeros[,.(count_spp = .N, sum_abun = sum(mean_density_m2)),.(Site, DepthZone, Region)]
kelp_depthzone_site_richness[,category:="kelp"]

#rbind these to make a single figure
depthzone_site_richness <- rbind(fish_depthzone_site_richness_abun_biomass, macroinvert_depthzone_site_richness, kelp_depthzone_site_richness, fill = TRUE)

##################################################
#Diversity metric visualizations
##################################################

depthzone_site_richness[,category := factor(category, labels = c("Fish","Macroalgae","Macroinvertebrates"))]

#boxplot
pal_5 <- c("lightsalmon1", "gold1", "palegreen4","mediumpurple1","skyblue")

depthzone_site_richness_boxplot <- ggplot(depthzone_site_richness, aes(x = DepthZone, y = count_spp, fill = DepthZone)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_5) +
  scale_x_discrete(labels = c("Inner \n (n = 64)", "Middle \n (n = 63)", "Outer \n (n = 55)", "Deep \n (n = 25)", "ARM \n (n = 25)")) +
  labs(x = "Depth Zone",
       y = "No taxa") +
  facet_grid(~category)+
  theme_classic()+
  theme(legend.position = "none")

ggsave(depthzone_site_richness_boxplot, path = "figures", filename = "depthzone_site_richness_boxplot.jpg", height = 6, width = 8, unit = "in")

#alternatively, box plot also split by ARM, island, mainland
depthzone_site_richness[,type := ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")]

depthzone_site_type_richness_boxplot <- ggplot(depthzone_site_richness, aes(x = DepthZone, y = count_spp, fill = type)) +
  geom_boxplot(outlier.size = 1, position = position_dodge(preserve = "single")) + #the median, two hinges and two whiskers, and all "outlying" points individually
  scale_fill_manual(values = c("grey","white"), ) +
  scale_x_discrete(labels = c("Inner", "Middle", "Outer", "Deep", "ARM")) +
  labs(x = "Depth Zone",
       y = "No taxa",
       fill = "Site location") +
  facet_grid(~category)+
  theme_classic()+
  theme()

ggsave(depthzone_site_type_richness_boxplot, path = "figures", filename = "depthzone_site_type_richness_boxplot.jpg", height = 5, width = 8, unit = "in")

#biomass and abundance boxplots
#abundance macro kelp
depthzone_site_type_abundance_macro_kelp_boxplot <- ggplot(depthzone_site_richness[category != "Fish"], aes(x = DepthZone, y = sum_abun, fill = type)) +
  geom_boxplot(outlier.size = 1, position = position_dodge(preserve = "single")) + #the median, two hinges and two whiskers, and all "outlying" points individually
  scale_fill_manual(values = c("grey","white"), ) +
  scale_x_discrete(labels = c("Inner", "Middle", "Outer", "Deep", "ARM")) +
  labs(x = "Depth Zone",
       y = "Abundance summed across all taxa\n(count per m^2)",
       fill = "Site location") +
  facet_wrap(~category, scales = "free_y")+
  theme_classic()+
  theme(legend.position = "bottom", legend.justification = "center", legend.direction = "horizontal", panel.spacing = unit(3, "lines"), axis.line = element_line())  +
  scale_y_continuous(limits=c(0,6.5))

#abundance fish
depthzone_site_type_abundance_fish_boxplot <- ggplot(depthzone_site_richness[category == "Fish"], aes(x = DepthZone, y = sum_abun, fill = type)) +
  geom_boxplot(outlier.size = 1, position = position_dodge(preserve = "single")) + #the median, two hinges and two whiskers, and all "outlying" points individually
  scale_fill_manual(values = c("grey","white"), ) +
  scale_x_discrete(labels = c("Inner", "Middle", "Outer", "Deep", "ARM")) +
  labs(x = "Depth Zone",
       y = "Abundance summed across all taxa\n(count per m^2)",
       fill = "Site location") +
  facet_grid(~category)+
  theme_classic()+
  theme(legend.position = "null")

#biomass
depthzone_site_type_biomass_boxplot <- ggplot(depthzone_site_richness[category == "Fish"], aes(x = DepthZone, y = sum_biomass/1000, fill = type)) +
  geom_boxplot(outlier.size = 1, position = position_dodge(preserve = "single")) + #the median, two hinges and two whiskers, and all "outlying" points individually
  scale_fill_manual(values = c("grey","white"), ) +
  scale_x_discrete(labels = c("Inner", "Middle", "Outer", "Deep", "ARM")) +
  labs(x = "Depth Zone",
       y = "Biomass summed across all taxa\n(kg per m^2)",
       fill = "Site location") +
  facet_grid(~category)+
  theme_classic()+
  theme(legend.position = "null")

depthzone_site_type_fish_biomass_abun_boxplot_merge <- plot_grid(depthzone_site_type_abundance_fish_boxplot,
                                                                 depthzone_site_type_biomass_boxplot, ncol = 2)

depthzone_site_type_biomass_abun_boxplot_merge <- plot_grid(depthzone_site_type_fish_biomass_abun_boxplot_merge,
                                                            depthzone_site_type_abundance_macro_kelp_boxplot , ncol = 1)

ggsave(depthzone_site_type_biomass_abun_boxplot_merge, path = "figures", filename = "depthzone_site_type_biomass_abun_boxplot_merge.jpg", height = 8, width = 9, unit = "in")

##################################################################################################
#Long data to wide data for multivariate community visualizations and analyses
##################################################################################################

#melt long to wide
dat_fish_averages_bysite.wide <- dcast(dat_fish_long_density_removezeros, Region + Site + DepthZone ~ Species, value.var = "mean_density_m2", fun = mean, fill = 0)

dat_fish_biomass_averages_bysite.wide <- dcast(dat_fish_long_biomass_removezeros, Region + Site + DepthZone ~ Species, value.var = "mean_wt_density_g_m2", fun = mean, fill = 0)

dat_macroinvert_averages_bysite.wide <- dcast(dat_macroinvert_long_density_removezeros, Region + Site + DepthZone ~ BenthicReefSpecies, value.var = "mean_density_m2", fun = mean, fill = 0)

dat_kelp_averages_bysite.wide <- dcast(dat_kelp_long_density_removezeros, Region + Site + DepthZone ~ BenthicReefSpecies, value.var = "mean_density_m2", fun = mean, fill = 0)

#Merge for all density data tables for species visualization (biomass of fish kept separate)
dat_averages_bysite.wide <- dat_fish_averages_bysite.wide[dat_macroinvert_averages_bysite.wide, on = c("Region","Site","DepthZone")]
dat_averages_bysite.wide <- dat_kelp_averages_bysite.wide[dat_averages_bysite.wide, on = c("Region","Site","DepthZone")]

#there are some sites that have no kelp, and these come up as NAs, change to 0 
dat_averages_bysite.wide[is.na(dat_averages_bysite.wide)] <- 0

####Add variables to single data table
#for all
dat_averages_bysite.wide.envir <- all_env_lat_lon[dat_averages_bysite.wide, on = c("Site","DepthZone")]
#for fish abun only
dat_averages_fishdensity_bysite.wide.envir <- all_env_lat_lon[dat_fish_averages_bysite.wide, on = c("Site","DepthZone")]
#for fish biomass only
dat_averages_fishbiomass_bysite.wide.envir <- all_env_lat_lon[dat_fish_biomass_averages_bysite.wide, on = c("Site","DepthZone")]
#for kelp only
dat_averages_kelpdensity_bysite.wide.envir <- all_env_lat_lon[dat_kelp_averages_bysite.wide, on = c("Site","DepthZone")]
#for macroinvert only
dat_averages_macroinvertdensity_bysite.wide.envir <- all_env_lat_lon[dat_macroinvert_averages_bysite.wide, on = c("Site","DepthZone")]

#ALTERNATIVELY, to match gllvm package example
#species only
site_spp <- dat_averages_bysite.wide.envir[,c(21:234)] #check these #s
site_env <- dat_averages_bysite.wide.envir[,c(2:7,9:17)] #check these #s, AND NOTE I EXCLUDED DEPTH ZONE
site_env.s <- scale(site_env)
#merge with depth zone
site_env.f <- cbind(site_env.s, site_depth)
site_depth <- dat_averages_bysite.wide.envir[,8]