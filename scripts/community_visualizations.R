# CREATION DATE 28 Jan 2023

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

source(file.path("functions","return_spptaxonomy_function.R"))

########################
##Load data
########################
dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))
dat_fish_site_averages <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages.rds"))
dat_macroinvert_site_averages <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages.rds"))
dat_kelp_site_averages <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages.rds"))

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
dat_macroinvert_long_density_removezeros <- dat_macroinvert_long_density[mean_density_m2>0,]
dat_kelp_long_density_removezeros <- dat_kelp_long_density[mean_density_m2>0,]


#number of species per depth zone?
fish_depthzone_site_richness <- dat_fish_long_density_removezeros[,.(count_spp = .N),.(Site, DepthZone)]
fish_depthzone_site_richness[,category:="fish"]
macroinvert_depthzone_site_richness <- dat_macroinvert_long_density_removezeros[,.(count_spp = .N),.(Site, DepthZone)]
macroinvert_depthzone_site_richness[,category:="macroinvert"]
kelp_depthzone_site_richness <- dat_kelp_long_density_removezeros[,.(count_spp = .N),.(Site, DepthZone)]
kelp_depthzone_site_richness[,category:="kelp"]

#rbind these to make a single figure
depthzone_site_richness <- rbind(fish_depthzone_site_richness, macroinvert_depthzone_site_richness, kelp_depthzone_site_richness)


#boxplot
pal_5 <- c("lightsalmon1", "gold1", "palegreen4","mediumpurple1","skyblue")

depthzone_site_richness_boxplot <- ggplot(depthzone_site_richness, aes(x = DepthZone, y = count_spp, fill = DepthZone)) +
  geom_boxplot() +
  scale_fill_manual(values = pal_5) +
  scale_x_discrete(labels = c("Inner \n (n = 64)", "Middle \n (n = 63)", "Outer \n (n = 55)", "Deep \n (n = 25)", "ARM \n (n = 25)")) +
  labs(x = "Depth Zone",
       y = "Number of species per site",
       title = "Species richness") +
  facet_grid(~category)+
  theme_classic()+
  theme(legend.position = "none")

ggsave(depthzone_site_richness_boxplot, path = "figures", filename = "depthzone_site_richness_boxplot.jpg", height = 6, width = 8, unit = "in")

#Now, simple NMDS

#Permutational Multivariate Analysis of Variance (perMANOVA)
dat_fish_wide_density.trim <- dat_fish_wide_density[,4:ncol(dat_fish_wide_density)]
dat_fish_wide_density.env <- dat_fish_wide_density[,1:3]

fish_perm <- adonis2(dat_fish_wide_density.trim ~ DepthZone, dat_fish_wide_density.env)
fish_perm

#Fish communities differ by depthzone

dat_macroinvert_wide_density.trim <- dat_macroinvert_wide_density[,4:ncol(dat_macroinvert_wide_density)]
dat_macroinvert_wide_density.env <- dat_macroinvert_wide_density[,1:3]

macroinvert_perm <- adonis2(dat_macroinvert_wide_density.trim ~ DepthZone, dat_macroinvert_wide_density.env)
macroinvert_perm

#sometimes, no kelp at all, need to delete these rows
kelp_row_sums <- rowSums(dat_kelp_wide_density[,4:ncol(dat_kelp_wide_density)])

dat_kelp_wide_density[,rowSums := kelp_row_sums]

dat_kelp_wide_density.r <- dat_kelp_wide_density[rowSums>0]

dat_kelp_wide_density.trim <- dat_kelp_wide_density.r[,4:21]
dat_kelp_wide_density.env <- dat_kelp_wide_density.r[,1:3]

kelp_perm <- adonis2(dat_kelp_wide_density.trim ~ DepthZone, dat_kelp_wide_density.env)
kelp_perm

#Now, PCA
fishdensityPCA <- rda(dat_fish_wide_density.trim)
fishdensityPCA

fishdensityPCAscores <- scores(fishdensityPCA, display = "sites") %>% 
  cbind(dat_fish_wide_density.env) %>%
  as.data.table()

fishdensityPCAvect <- scores(fishdensityPCA, display = "species") %>% 
  as.data.table() %>%
  

plot_fishdensity_PCA <- ggplot() +
  geom_point(data = fishdensityPCAscores, aes(x = PC1, y = PC2, color = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = fishdensityPCAvect, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = fishdensityPCAvect, aes(x = PC1, y = PC2, label = rownames(fishdensityPCAvect))) +
  labs(x = "PC1 (47.12%)",
       y = "PC2 (11.96%)",
       title = "Principal Components Analysis") +
  theme_bw()
plot_fishdensity_PCA

#this clearly isn't looking right, maybe transform data before starting?

#To try next
##use species presence/absence instead of abundance
##transform, like Jeremy recommended
##use Jeremy's code
