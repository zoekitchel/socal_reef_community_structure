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
library(labdsv)

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
#FISH ONLY
#Permutational Multivariate Analysis of Variance (perMANOVA)
dat_fish_wide_density.trim <- dat_fish_wide_density[,4:ncol(dat_fish_wide_density)]
dat_fish_wide_density.env <- dat_fish_wide_density[,1:3]

#Hellinger transformation
###a method of pre-transforming species composition data for the use in linear ordination methods, 
###resulting in transformation-based ordination analysis (tb-PCA, tb-RDA).
###Calculates a sample total standardization (all values in a row are divided by the row sum),
###and then takes the square root of the values.

dat_fish_wide_density.hellin <- hellinger(dat_fish_wide_density.trim)

fishdensity_perm <- adonis2(dat_fish_wide_density.hellin ~ DepthZone, dat_fish_wide_density.env)
fishdensity_perm

#Fish communities differ by depthzone

#Now, PCA
fishdensityPCA <- rda(dat_fish_wide_density.hellin)
fishdensityPCA

fishdensityPCAscores <- vegan::scores(fishdensityPCA, display = "sites") %>% 
  cbind(dat_fish_wide_density.env) %>%
  as.data.table()

fishdensityPCAvect <- vegan::scores(fishdensityPCA, display = "species") %>% 
  as.data.table()

fishdensityPCAvect[,"Species" := colnames(dat_fish_wide_density.trim)]

fishdensity_eig1 <- round(fishdensityPCA$CA$eig[1]*100,1)
fishdensity_eig2 <- round(fishdensityPCA$CA$eig[2]*100,1)

#no labels, and excluding ARMs
plot_fishdensity_PCA_nolabel_noarm <- ggplot() +
  geom_point(data = fishdensityPCAscores[DepthZone != "ARM"], aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  labs(x = paste0("PC1 (",fishdensity_eig1,"%)"),
       y = paste0("PC2 (",fishdensity_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-0.4,0.45) +
  ylim(-.8,0.25) +
  theme_bw()
plot_fishdensity_PCA_nolabel_noarm

ggsave(plot_fishdensity_PCA_nolabel_noarm, path = "figures", filename = "plot_fishdensity_PCA_nolabel_noarm.jpg", height = 4, width = 5)

#no labels
plot_fishdensity_PCA_nolabel <- ggplot() +
  geom_point(data = fishdensityPCAscores, aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  labs(x = paste0("PC1 (",fishdensity_eig1,"%)"),
       y = paste0("PC2 (",fishdensity_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-0.4,0.45) +
  ylim(-.8,0.25) +
  theme_bw()
plot_fishdensity_PCA_nolabel

ggsave(plot_fishdensity_PCA_nolabel, path = "figures", filename = "plot_fishdensity_PCA_nolabel.jpg", height = 4, width = 5)

#vectors and all data
plot_fishdensity_PCA <- ggplot() +
  geom_point(data = fishdensityPCAscores, aes(x = PC1, y = PC2, color = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = fishdensityPCAvect, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = fishdensityPCAvect, aes(x = PC1, y = PC2, label = Species), size = 3, nudge_y = 0.05) +
  labs(x = paste0("PC1 (",fishdensity_eig1,"%)"),
       y = paste0("PC2 (",fishdensity_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-.9,1.7) +
  theme_bw()
plot_fishdensity_PCA

ggsave(plot_fishdensity_PCA, path = "figures", filename = "plot_fishdensity_PCA.jpg", height = 4, width = 5)

View(fishdensityPCAvect[,abs_PC1:= abs(PC1)][,abs_PC2:= abs(PC2)])

#only show top species
top_fishdensity_spp <- c("Brachyistius frenatus",
                  "Hypsypops rubicundus",
                  "Lythrypnus dalli",
                  "Chromis punctipinnis",
                  "Oxyjulis californica",
                  "Halichoeres semicinctus",
                  "Semicossyphus pulcher")

plot_fishdensity_PCA_top5label <- ggplot() +
  geom_point(data = fishdensityPCAscores, aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone), size = 1) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = fishdensityPCAvect[Species %in% top_fishdensity_spp,], aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = fishdensityPCAvect[Species %in% top_fishdensity_spp,], aes(x = PC1, y = PC2, label = Species), size = 3, nudge_y = 0.05) +
  labs(x = paste0("PC1 (",fishdensity_eig1,"%)"),
       y = paste0("PC2 (",fishdensity_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-.9,1.7) +
  theme_bw()
plot_fishdensity_PCA_top5label

ggsave(plot_fishdensity_PCA_top5label, path = "figures", filename = "plot_fishdensity_PCA_top5label.jpg", height = 4, width = 5)

################
#FISH BIOMASS
################
#FISH ONLY
#Permutational Multivariate Analysis of Variance (perMANOVA)
dat_fish_wide_biomass.trim <- dat_fish_wide_biomass[,4:ncol(dat_fish_wide_biomass)]
dat_fish_wide_biomass.env <- dat_fish_wide_biomass[,1:3]

#Hellinger transformation

dat_fish_wide_biomass.hellin <- hellinger(dat_fish_wide_biomass.trim)

fishbiomass_perm <- adonis2(dat_fish_wide_biomass.hellin ~ DepthZone, dat_fish_wide_biomass.env)
fishbiomass_perm

#Fish communities differ by depthzone

#Now, PCA
fishbiomassPCA <- rda(dat_fish_wide_biomass.hellin)
fishbiomassPCA

fishbiomassPCAscores <- vegan::scores(fishbiomassPCA, display = "sites") %>% 
  cbind(dat_fish_wide_biomass.env) %>%
  as.data.table()

fishbiomassPCAvect <- vegan::scores(fishbiomassPCA, display = "species") %>% 
  as.data.table()

fishbiomassPCAvect[,"Species" := colnames(dat_fish_wide_biomass.trim)]

fishbiomass_eig1 <- round(fishbiomassPCA$CA$eig[1]*100,1)
fishbiomass_eig2 <- round(fishbiomassPCA$CA$eig[2]*100,1)

#no labels, and excluding ARMs
plot_fishbiomass_PCA_nolabel_noarm <- ggplot() +
  geom_point(data = fishbiomassPCAscores[DepthZone != "ARM"], aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  labs(x = paste0("PC1 (",fishbiomass_eig1,"%)"),
       y = paste0("PC2 (",fishbiomass_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-0.43,0.40) +
  ylim(-.51,0.72) +
  theme_bw()
plot_fishbiomass_PCA_nolabel_noarm

ggsave(plot_fishbiomass_PCA_nolabel_noarm, path = "figures", filename = "plot_fishbiomass_PCA_nolabel_noarm.jpg", height = 4, width = 5)

#no labels
plot_fishbiomass_PCA_nolabel <- ggplot() +
  geom_point(data = fishbiomassPCAscores, aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  labs(x = paste0("PC1 (",fishbiomass_eig1,"%)"),
       y = paste0("PC2 (",fishbiomass_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-0.43,0.40) +
  ylim(-.51,0.72) +
  theme_bw()
plot_fishbiomass_PCA_nolabel

ggsave(plot_fishbiomass_PCA_nolabel, path = "figures", filename = "plot_fishbiomass_PCA_nolabel.jpg", height = 4, width = 5)

#vectors and all data
plot_fishbiomass_PCA <- ggplot() +
  geom_point(data = fishbiomassPCAscores, aes(x = PC1, y = PC2, color = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = fishbiomassPCAvect, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = fishbiomassPCAvect, aes(x = PC1, y = PC2, label = Species), size = 3, nudge_y = 0.05) +
  labs(x = paste0("PC1 (",fishbiomass_eig1,"%)"),
       y = paste0("PC2 (",fishbiomass_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-.55,.9) +
  theme_bw()
plot_fishbiomass_PCA

ggsave(plot_fishbiomass_PCA, path = "figures", filename = "plot_fishbiomass_PCA.jpg", height = 4, width = 5)

View(fishbiomassPCAvect[,abs_PC1:= abs(PC1)][,abs_PC2:= abs(PC2)])

#only show top species
top_fishbiomass_spp <- c("Paralabrax nebulifer",
                         "Chromis punctipinnis",
                         "Stereolepis gigas",
                         "Girella nigricans",
                         "Hypsypops rubicundus",
                         "Paralabrax clathratus")

plot_fishbiomass_PCA_top5label <- ggplot() +
  geom_point(data = fishbiomassPCAscores, aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone), size = 1) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = fishbiomassPCAvect[Species %in% top_fishbiomass_spp,], aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = fishbiomassPCAvect[Species %in% top_fishbiomass_spp,], aes(x = PC1, y = PC2, label = Species), size = 3, nudge_y = 0.05) +
  labs(x = paste0("PC1 (",fishbiomass_eig1,"%)"),
       y = paste0("PC2 (",fishbiomass_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-.7,1.1) +
  theme_bw()
plot_fishbiomass_PCA_top5label

ggsave(plot_fishbiomass_PCA_top5label, path = "figures", filename = "plot_fishbiomass_PCA_top5label.jpg", height = 4, width = 5)

#################
#MACROINVERT ONLY
#################
dat_macroinvert_wide_density.trim <- dat_macroinvert_wide_density[,4:ncol(dat_macroinvert_wide_density)]
dat_macroinvert_wide_density.env <- dat_macroinvert_wide_density[,1:3]

#Hellinger transformation
dat_macroinvert_wide_density.hellin <- hellinger(dat_macroinvert_wide_density.trim)

macroinvert_perm <- adonis2(dat_macroinvert_wide_density.hellin ~ DepthZone, dat_macroinvert_wide_density.env)
macroinvert_perm

#significantly different

#Now, PCA
macroinvertdensityPCA <- rda(dat_macroinvert_wide_density.hellin)
macroinvertdensityPCA

macroinvertdensityPCAscores <- vegan::scores(macroinvertdensityPCA, display = "sites") %>% 
  cbind(dat_macroinvert_wide_density.env) %>%
  as.data.table()

macroinvertdensityPCAvect <- vegan::scores(macroinvertdensityPCA, display = "species") %>% 
  as.data.table()

macroinvertdensityPCAvect[,"Species" := colnames(dat_macroinvert_wide_density.trim)]

macroinvert_eig1 <- round(macroinvertdensityPCA$CA$eig[1]*100,1)
macroinvert_eig2 <- round(macroinvertdensityPCA$CA$eig[2]*100,1)

#no labels, and excluding ARMs
plot_macroinvertdensity_PCA_nolabel_noarm <- ggplot() +
  geom_point(data = macroinvertdensityPCAscores[DepthZone != "ARM"], aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  labs(x = paste0("PC1 (",macroinvert_eig1,"%)"),
       y = paste0("PC2 (",macroinvert_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-0.5,0.4) +
  theme_bw()
plot_macroinvertdensity_PCA_nolabel_noarm

ggsave(plot_macroinvertdensity_PCA_nolabel_noarm, path = "figures", filename = "plot_macroinvertdensity_PCA_nolabel_noarm.jpg", height = 4, width = 5)

#no labels
plot_macroinvertdensity_PCA_nolabel <- ggplot() +
  geom_point(data = macroinvertdensityPCAscores, aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  labs(x = paste0("PC1 (",macroinvert_eig1,"%)"),
       y = paste0("PC2 (",macroinvert_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-0.5,0.4) +
  theme_bw()
plot_macroinvertdensity_PCA_nolabel

ggsave(plot_macroinvertdensity_PCA_nolabel, path = "figures", filename = "plot_macroinvertdensity_PCA_nolabel.jpg", height = 4, width = 5)
  
#vectors and all data
  plot_macroinvertdensity_PCA <- ggplot() +
  geom_point(data = macroinvertdensityPCAscores, aes(x = PC1, y = PC2, color = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = macroinvertdensityPCAvect, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = macroinvertdensityPCAvect, aes(x = PC1, y = PC2, label = Species), size = 3, nudge_y = 0.05) +
  labs(x = paste0("PC1 (",macroinvert_eig1,"%)"),
       y = paste0("PC2 (",macroinvert_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
    xlim(-1,1.2) +
  theme_bw()
plot_macroinvertdensity_PCA

ggsave(plot_macroinvertdensity_PCA, path = "figures", filename = "plot_macroinvertdensity_PCA.jpg", height = 4, width = 5)

View(macroinvertdensityPCAvect[,abs_PC1:= abs(PC1)][,abs_PC2:= abs(PC2)])

#only show top species
top_macroinvert_spp <- c("Strongylocentrotus purpuratus", "Megastraea undosa", "Centrostephanus coronatus", "Leptogorgia chilensis",
                         "Mesocentrotus franciscanus", "Kelletia kelletii", "Panulirus interruptus")

plot_macroinvertdensity_PCA_top5label <- ggplot() +
  geom_point(data = macroinvertdensityPCAscores, aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone), size = 1) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = macroinvertdensityPCAvect[Species %in% top_macroinvert_spp,], aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = macroinvertdensityPCAvect[Species %in% top_macroinvert_spp,], aes(x = PC1, y = PC2, label = Species), size = 3, nudge_y = 0.05) +
  labs(x = paste0("PC1 (",macroinvert_eig1,"%)"),
       y = paste0("PC2 (",macroinvert_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-1,1.5) +
  theme_bw()
plot_macroinvertdensity_PCA_top5label

ggsave(plot_macroinvertdensity_PCA_top5label, path = "figures", filename = "plot_macroinvertdensity_PCA_top5label.jpg", height = 4, width = 5)

#################
#KELP ONLY
#################
#sometimes, no kelp at all, need to delete these rows
kelp_row_sums <- rowSums(dat_kelp_wide_density[,4:ncol(dat_kelp_wide_density)])

dat_kelp_wide_density[,rowSums := kelp_row_sums]

dat_kelp_wide_density.r <- dat_kelp_wide_density[rowSums>0]

dat_kelp_wide_density.trim <- dat_kelp_wide_density.r[,4:21]
dat_kelp_wide_density.env <- dat_kelp_wide_density.r[,1:3]

#Hellinger transformation
dat_kelp_wide_density.hellin <- hellinger(dat_kelp_wide_density.trim)

kelp_perm <- adonis2(dat_kelp_wide_density.trim ~ DepthZone, dat_kelp_wide_density.env)
kelp_perm

#Now, PCA
kelpdensityPCA <- rda(dat_kelp_wide_density.hellin)
kelpdensityPCA

kelpdensityPCAscores <- vegan::scores(kelpdensityPCA, display = "sites") %>% 
  cbind(dat_kelp_wide_density.env) %>%
  as.data.table()

kelpdensityPCAvect <- vegan::scores(kelpdensityPCA, display = "species") %>% 
  as.data.table()

kelpdensityPCAvect[,"Species" := colnames(dat_kelp_wide_density.trim)]

kelp_eig1 <- round(kelpdensityPCA$CA$eig[1]*100,1)
kelp_eig2 <- round(kelpdensityPCA$CA$eig[2]*100,1)

#no labels, and excluding ARMs
plot_kelpdensity_PCA_nolabel_noarm <- ggplot() +
  geom_point(data = kelpdensityPCAscores[DepthZone != "ARM"], aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  labs(x = paste0("PC1 (",kelp_eig1,"%)"),
       y = paste0("PC2 (",kelp_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-0.35,0.5) +
  theme_bw()
plot_kelpdensity_PCA_nolabel_noarm

ggsave(plot_kelpdensity_PCA_nolabel_noarm, path = "figures", filename = "plot_kelpdensity_PCA_nolabel_noarm.jpg", height = 4, width = 5)

#no labels
plot_kelpdensity_PCA_nolabel <- ggplot() +
  geom_point(data = kelpdensityPCAscores, aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  labs(x = paste0("PC1 (",kelp_eig1,"%)"),
       y = paste0("PC2 (",kelp_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-0.35,0.5) +
  theme_bw()
plot_kelpdensity_PCA_nolabel

ggsave(plot_kelpdensity_PCA_nolabel, path = "figures", filename = "plot_kelpdensity_PCA_nolabel.jpg", height = 4, width = 5)

#vectors and all data
plot_kelpdensity_PCA <- ggplot() +
  geom_point(data = kelpdensityPCAscores, aes(x = PC1, y = PC2, color = DepthZone)) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = kelpdensityPCAvect, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = kelpdensityPCAvect, aes(x = PC1, y = PC2, label = Species), size = 3, nudge_y = 0.05) +
  labs(x = paste0("PC1 (",kelp_eig1,"%)"),
       y = paste0("PC2 (",kelp_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-1,1.9) +
  theme_bw()
plot_kelpdensity_PCA

ggsave(plot_kelpdensity_PCA, path = "figures", filename = "plot_kelpdensity_PCA.jpg", height = 4, width = 5)

View(kelpdensityPCAvect[,abs_PC1:= abs(PC1)][,abs_PC2:= abs(PC2)])

#only show top species
top_kelp_spp <- c("Laminaria farlowii", "Sargassum horneri", "Eisenia arborea", "Stephanocystis spp.,
                  Macrocystis pyrifera", "Pterygophora californica", "Sargassum palmeri")

plot_kelpdensity_PCA_top5label <- ggplot() +
  geom_point(data = kelpdensityPCAscores, aes(x = PC1, y = PC2, color = DepthZone, shape = DepthZone), size = 1) +
  scale_color_manual(values = pal_5) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  geom_segment(data = kelpdensityPCAvect[Species %in% top_kelp_spp,], aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(length = unit(0.2, "cm"))) +
  geom_text(data = kelpdensityPCAvect[Species %in% top_kelp_spp,], aes(x = PC1, y = PC2, label = Species), size = 3, nudge_y = 0.05) +
  labs(x = paste0("PC1 (",kelp_eig1,"%)"),
       y = paste0("PC2 (",kelp_eig2,"%)"),
       title = "PCA: Hellinger transformation") +
  xlim(-1,1.9) +
  theme_bw()
plot_kelpdensity_PCA_top5label

ggsave(plot_kelpdensity_PCA_top5label, path = "figures", filename = "plot_kelpdensity_PCA_top5label.jpg", height = 4, width = 5)

#################
#MERGE ALL SPECIES
#################
dat_averages_bysite <- rbind(dat_fish_averages_bysite, dat_macroinvert_averages_bysite, dat_kelp_averages_bysite, use.names = FALSE)

#spp per depthzone
dat_fish_long_density_removezeros <- dat_fish_long_density[mean_density_m2>0,]
dat_macroinvert_long_density_removezeros <- dat_macroinvert_long_density[mean_density_m2>0,]
dat_kelp_long_density_removezeros <- dat_kelp_long_density[mean_density_m2>0,]

##################################################
#Long data to wide data for vegan analyses
##################################################

#merge


#melt long to wide
dat_fish_averages_bysite.wide <- dcast(dat_fish_averages_bysite, Region + Site + DepthZone + `Site type` ~ Species, value.var = "mean_depthzone_density_m2", fun = mean)

dat_macroinvert_averages_bysite.wide <- dcast(dat_macroinvert_averages_bysite, Region + Site + DepthZone + `Site type` ~ BenthicReefSpecies, value.var = "mean_depthzone_density_m2", fun = mean)

dat_kelp_averages_bysite.wide <- dcast(dat_kelp_averages_bysite, Region + Site + DepthZone + `Site type` ~ BenthicReefSpecies, value.var = "mean_depthzone_density_m2", fun = mean)

dat_averages_bysite.wide <- dcast(dat_averages_bysite, Region + Site + DepthZone + `Site type` ~ Species, value.var = "mean_depthzone_density_m2", fun = mean)

#####################
#PERMANOVA, Permutational Multivariate Analysis of Variance (perMANOVA)
#####################
#FULL COMMUNITY

#Start with root transformation (high abundance spp less pull in ordination)
dat_averages_bysite.wide.rt <- sqrt(dat_averages_bysite.wide[,c(5:ncol(dat_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_averages_bysite.wide.rt <- cbind(dat_averages_bysite.wide[,c(1:4), with = FALSE], dat_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_averages_bysite.wide.rt.trim <- dat_averages_bysite.wide.rt[,c(5:ncol(dat_averages_bysite.wide.rt)), with = FALSE]

permanova_allspp <- adonis2(
  dat_averages_bysite.wide.rt.trim ~ dat_averages_bysite.wide.rt $DepthZone,
  method = "bray"
)

permanova_allspp

#whether or not a site is an artificial reef affects community composition, but
#accounts for relatively small amount of variation (12.3% of variation in composition)

#ONLY FISH

#Start with root transformation (high abundance spp less pull in ordination)
dat_fish_averages_bysite.wide.rt <- sqrt(dat_fish_averages_bysite.wide[,c(5:ncol(dat_fish_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_fish_averages_bysite.wide.rt <- cbind(dat_fish_averages_bysite.wide[,c(1:4), with = FALSE], dat_fish_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_fish_averages_bysite.wide.rt.trim <- dat_fish_averages_bysite.wide.rt[,c(5:ncol(dat_fish_averages_bysite.wide.rt)), with = FALSE]

permanova_fish <- adonis2(
  dat_fish_averages_bysite.wide.rt.trim ~ dat_fish_averages_bysite.wide.rt $DepthZone,
  method = "bray"
)

permanova_fish

#whether or not a site is an artificial reef affects fish community composition, accounting for 14.5% of variation in composition

#ONLY MACRO

#Start with root transformation (high abundance spp less pull in ordination)
dat_macroinvert_averages_bysite.wide.rt <- sqrt(dat_macroinvert_averages_bysite.wide[,c(5:ncol(dat_macroinvert_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_macroinvert_averages_bysite.wide.rt <- cbind(dat_macroinvert_averages_bysite.wide[,c(1:4), with = FALSE], dat_macroinvert_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_macroinvert_averages_bysite.wide.rt.trim <- dat_macroinvert_averages_bysite.wide.rt[,c(5:ncol(dat_macroinvert_averages_bysite.wide.rt)), with = FALSE]

permanova_macroinvert <- adonis2(
  dat_macroinvert_averages_bysite.wide.rt.trim ~ dat_macroinvert_averages_bysite.wide.rt $DepthZone,
  method = "bray"
)

permanova_macroinvert

#whether or not a site is an artificial reef affects macroinvert community composition, but only accounting for 10.5% of variation in composition


#ONLY KELP

#sometimes, no kelp at all, need to delete these rows
kelp_row_sums <- rowSums(dat_kelp_averages_bysite.wide[,5:ncol(dat_kelp_wide_density)])

dat_kelp_averages_bysite.wide[,rowSums := kelp_row_sums]

dat_kelp_averages_bysite.wide.r <- dat_kelp_averages_bysite.wide[rowSums>0]

#Start with root transformation (high abundance spp less pull in ordination)
dat_kelp_averages_bysite.wide.rt <- sqrt(dat_kelp_averages_bysite.wide.r[,c(5:ncol(dat_kelp_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_kelp_averages_bysite.wide.rt <- cbind(dat_kelp_averages_bysite.wide.r[,c(1:4), with = FALSE], dat_kelp_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_kelp_averages_bysite.wide.rt.trim <- dat_kelp_averages_bysite.wide.rt[,c(5:ncol(dat_kelp_averages_bysite.wide.rt)), with = FALSE]

permanova_kelp <- adonis2(
  dat_kelp_averages_bysite.wide.rt.trim ~ dat_kelp_averages_bysite.wide.rt $DepthZone,
  method = "bray"
)

permanova_kelp

#whether or not a site is an artificial reef affects fish community composition, accounting for 14.2% of variation in composition

#Take away:
#Artificial vs natural accounts for the most variation for fish and kelp communities, and the least relative variation for macroinvert communities


###################
#PLOT NMDS
###################
#NMDS “preserves the rank order of the inter-point dissimilarities as well as possible within the constraints of a small number of dimensions” (Anderson et al. 2008, p. 105).
#It has the least restrictive assumptions (any distance measure, no assumptions about shape of species response curves).

#all species

#convert to distance matrix
dat_averages_bysite.dist <- vegdist(dat_averages_bysite.wide[,5:ncol(dat_averages_bysite.wide)], distance = "bray")

full_nmds <- metaMDS(dat_averages_bysite.wide[,5:ncol(dat_averages_bysite.wide)], #metaMDS automatically does this for us
                     autotransform = FALSE,
                     distance = "bray",
                     engine = "monoMDS",
                     k = 4, #need at least twice as many sample units as dimensions, will impact results
                     weakties = TRUE,#Set to true if you commonly have sample units that don't share species 
                     #TRUE allows observations to be positioned at different coordinates within the ordination space even though they are the same distance apart in the matrix
                     #FALSE forces observations that are the same distance apart in the distance matrix to also be the same distance apart in ordination space
                     model = "global",
                     maxit = 300,
                     try = 40,
                     trymax = 100)

print(full_nmds)


#The number of dimensions is often selected by conducting NMDS ordinations with different 
#numbers of dimensions (e.g., k = 1 to 5) and then plotting stress as a function of k. 
#The intent here is to identify the point where adding the complexity of an additional dimension contributes relatively little to the reduction in stress.

#Stress: 
#         k = 2, stress = .168
#         k = 3, stress = .112
#         k = 4, stress = 0.07, good ordination with no real risk of drawing false inferences (Clark 1993 pg 126)
#         k = 5, stress = 0.06
#         k = 6, stress = 0.05
#         k = 3, stress = 0.04

#The fit of a NMDS ordination can be assessed by plotting the original dissimilarities (z$diss) against the (Euclidean) ordination distances (z$dist) 

stressplot(full_nmds)

#link NMDS points to site details
dat_averages_bysite.nmds <- cbind(dat_fish_averages_bysite.wide[,1:4], full_nmds$points)

#set apart SoS
dat_averages_bysite.nmds[Site == "Star of Scotland",`Site type` := "Star of Scotland\nn = 1"]

#For some reason, Marina Del Rey coded as Region = Artificial reef, change to Santa Monica Bay

dat_averages_bysite.nmds[Region == "Artificial Reef",Region := "Santa Monica Bay"]

full_nmds_plot_bysitetype <- ggplot(data = dat_averages_bysite.nmds) +
  geom_point(aes(x = MDS1, y = MDS2, fill = `Site type`, shape = `Site type`), size = 3, color = "transparent") +
  scale_fill_manual(values = c("#BF6992", "#54AEAA", "#EA3323")) +
  scale_shape_manual(values = c(21,24,24)) +
  lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = c(.15, .9))


ggsave(full_nmds_plot_bysitetype, path = "figures", filename = "full_nmds_plot_bysitetype.jpg", height = 6, width = 6)

full_nmds_plot_byregion <- ggplot(data = dat_averages_bysite.nmds) +
  geom_point(aes(x = MDS1, y = MDS2, color = Region, shape = `Site type`), size = 3) +
  scale_color_manual(values = c("#83752D", "#EAC4AA", "#D38EB6", "#F6D798","#76A6B8", "#2F624E")) +
  scale_shape_manual(values = c(19,17,17), guide = NULL) +
  lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = c(.2, .83),
        legend.background = element_blank(), legend.key = element_blank())


ggsave(full_nmds_plot_byregion, path = "figures", filename = "full_nmds_plot_byregion.jpg", height = 6, width = 6)


#SPECIES SCORES
species_scores <- data.table(vegan::scores(full_nmds, "species"))
spp_list <- colnames(dat_averages_bysite.wide[,5:ncol(dat_averages_bysite.wide)])
species_scores[,Species := spp_list]

#top and bottom three for NMDS1 & 2
top_all_spp <- c("Stereolepis gigas","Pachycerianthus fimbriatus","Crossata californica",
                 "Paralabrax maculatofasciatus","Sebastes constellatus","Heterostichus rostratus",
                 "Lythrypnus dalli","Haliotis fulgens","Centrostephanus coronatus",
                 "Berthella californica","Pomaulax gibberosus","Acanthodoris rhodoceras")

species_scores_top <- species_scores[Species %in% top_all_spp]

species_scores_top[,Species := factor(Species,
                                      labels = 
                                        c("Stereolepis gigas\ngiant sea bass","Pachycerianthus fimbriatus\ntube anenome","Crossata californica\nfrogsnail",
                                          "Paralabrax maculatofasciatus\nspotted sand bass","Sebastes constellatus\nstarry rockfish","Heterostichus rostratus\ngiant kelpfish",
                                          "Lythrypnus dalli\nblue-banded goby","Haliotis fulgens\ngreen abalone","Centrostephanus coronatus\ncrowned sea urchin",
                                          "Berthella californica\nsea slug spp.","Pomaulax gibberosus\nred turban snail","Acanthodoris rhodoceras\nsea slug spp."))]

#add to simple NMDS plot

full_nmds_plot_speciesscores <- ggplot(data = dat_averages_bysite.nmds) +
  geom_point(aes(x = MDS1, y = MDS2,shape = `Site type`), size = 3, alpha = 0.3) +
  geom_text_repel(data = species_scores_top, aes(x = NMDS1, y = NMDS2, label = Species), size = 2) +
  scale_shape_manual(values = c(19,17,17)) +
  lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = c(.15, .9))


ggsave(full_nmds_plot_speciesscores, path = "figures", filename = "full_nmds_plot_speciesscores.jpg", height = 5, width = 10)

#plot together
full_nmds_plot <- plot_grid(full_nmds_plot_bysitetype, full_nmds_plot_byregion,full_nmds_plot_speciesscores + theme(legend.position = "none"), ncol = 3)

ggsave(full_nmds_plot, path = "figures", filename = "full_nmds_plot.jpg", height = 5, width = 15)
ggsave(full_nmds_plot, path = "figures", filename = "full_nmds_plot.pdf", height = 5, width = 15)

#We can create ellipses around groups

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(factor(dat_averages_bysite.nmds$DepthZone))){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(dat_averages_bysite.nmds[dat_averages_bysite.nmds$DepthZone==g,],
                                                   veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

full_nmds_plot_bysitetype +
  geom_path(data=df_ell, aes(x=MDS1, y=MDS2,color=group), size=1, linetype=2) +
  scale_color_manual(values = c("#BF6992", "#54AEAA"), guide = NULL)

#NOT IMMEDIATELY SURE WHAT THIS ELLIPSE IS TELLING US



######################
#To try next
##repeat for fish, fish biomass, macro, kelp, and all species clumped together (maybe do this one by presence absence only)
##Integrate environmental data (distance to 200m isobath, insitu: depth, temp, rugosity or whatever)




##################################################################################
#Visual summaries for OSM poster
##################################################################################

########################
##Load data (deep and AR, 2022-2023)
########################
dat_event_OSM.r <- readRDS(file.path("data","processed_crane", "dat_event_OSM.r.rds"))
dat_fish_site_averages_OSM <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages_OSM.rds"))
dat_macroinvert_site_averages_OSM <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages_OSM.rds"))
dat_kelp_site_averages_OSM <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages_OSM.rds"))

########################
##Deep only
########################
dat_event_OSM.r.deep_ar <- dat_event_OSM.r[DepthZone %in% c("ARM","Deep")]
dat_fish_site_averages_OSM.deep_ar <- dat_fish_site_averages_OSM[DepthZone %in% c("ARM","Deep")]
dat_macroinvert_site_averages_OSM.deep_ar <- dat_macroinvert_site_averages_OSM[DepthZone %in% c("ARM","Deep")]
dat_kelp_site_averages_OSM.deep_ar <- dat_kelp_site_averages_OSM[DepthZone %in% c("ARM","Deep")]

#######################
##Counts of each site type
######################
count_natural <- length(unique(dat_event_OSM.r.deep_ar[DepthZone == "Deep",Site])) #27
count_artificial <- length(unique(dat_event_OSM.r.deep_ar[DepthZone == "ARM",Site])) #28

########################
##Average across years for each site
########################

dat_fish_averages_bysite <- dat_fish_site_averages_OSM.deep_ar[,.(mean_depthzone_density_m2 = mean(mean_density_m2),
                                                                   mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                                                .(Site, Region,DepthZone, Species)] 
dat_macroinvert_averages_bysite <- dat_macroinvert_site_averages_OSM.deep_ar[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                                              .(Site, Region, DepthZone, BenthicReefSpecies)]  
dat_kelp_averages_bysite <- dat_kelp_site_averages_OSM.deep_ar[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                                .(Site, Region,DepthZone, BenthicReefSpecies)]  

#site type column
dat_fish_averages_bysite[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]
dat_macroinvert_averages_bysite[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]
dat_kelp_averages_bysite[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]

#reorder these factors
dat_fish_averages_bysite[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
                                                          labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]
dat_macroinvert_averages_bysite[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
                                                labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]
dat_kelp_averages_bysite[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
                                                labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]

#add column for species type
dat_fish_averages_bysite[,spp_class := "fish"]
dat_macroinvert_averages_bysite[,spp_class := "macroinvert"]
dat_kelp_averages_bysite[,spp_class := "kelp"]

#reorder columns
dat_fish_averages_bysite <- dat_fish_averages_bysite[,.(Region, Site, DepthZone, `Site type`,spp_class, Species, mean_depthzone_density_m2)] #note, I'm excluding biomass for now!!
dat_macroinvert_averages_bysite <- dat_macroinvert_averages_bysite[,.(Region, Site, DepthZone, `Site type`,spp_class, BenthicReefSpecies, mean_depthzone_density_m2)]
dat_kelp_averages_bysite <- dat_kelp_averages_bysite[,.(Region, Site, DepthZone, `Site type`,spp_class, BenthicReefSpecies, mean_depthzone_density_m2)]

#merge into single data table to allow for analyses with all species

dat_averages_bysite <- rbind(dat_fish_averages_bysite, dat_macroinvert_averages_bysite, dat_kelp_averages_bysite, use.names = FALSE)

##################################################
#Long data to wide data for vegan analyses
##################################################


#melt long to wide
dat_fish_averages_bysite.wide <- dcast(dat_fish_averages_bysite, Region + Site + DepthZone + `Site type` ~ Species, value.var = "mean_depthzone_density_m2", fun = mean)

dat_macroinvert_averages_bysite.wide <- dcast(dat_macroinvert_averages_bysite, Region + Site + DepthZone + `Site type` ~ BenthicReefSpecies, value.var = "mean_depthzone_density_m2", fun = mean)

dat_kelp_averages_bysite.wide <- dcast(dat_kelp_averages_bysite, Region + Site + DepthZone + `Site type` ~ BenthicReefSpecies, value.var = "mean_depthzone_density_m2", fun = mean)

dat_averages_bysite.wide <- dcast(dat_averages_bysite, Region + Site + DepthZone + `Site type` ~ Species, value.var = "mean_depthzone_density_m2", fun = mean)

#####################
#PERMANOVA, Permutational Multivariate Analysis of Variance (perMANOVA)
#####################
#FULL COMMUNITY

#Start with root transformation (high abundance spp less pull in ordination)
dat_averages_bysite.wide.rt <- sqrt(dat_averages_bysite.wide[,c(5:ncol(dat_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_averages_bysite.wide.rt <- cbind(dat_averages_bysite.wide[,c(1:4), with = FALSE], dat_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_averages_bysite.wide.rt.trim <- dat_averages_bysite.wide.rt[,c(5:ncol(dat_averages_bysite.wide.rt)), with = FALSE]

permanova_allspp <- adonis2(
  dat_averages_bysite.wide.rt.trim ~ dat_averages_bysite.wide.rt $DepthZone,
  method = "bray"
)

permanova_allspp

#whether or not a site is an artificial reef affects community composition, but
#accounts for relatively small amount of variation (12.3% of variation in composition)

#ONLY FISH

#Start with root transformation (high abundance spp less pull in ordination)
dat_fish_averages_bysite.wide.rt <- sqrt(dat_fish_averages_bysite.wide[,c(5:ncol(dat_fish_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_fish_averages_bysite.wide.rt <- cbind(dat_fish_averages_bysite.wide[,c(1:4), with = FALSE], dat_fish_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_fish_averages_bysite.wide.rt.trim <- dat_fish_averages_bysite.wide.rt[,c(5:ncol(dat_fish_averages_bysite.wide.rt)), with = FALSE]

permanova_fish <- adonis2(
  dat_fish_averages_bysite.wide.rt.trim ~ dat_fish_averages_bysite.wide.rt $DepthZone,
  method = "bray"
)

permanova_fish

#whether or not a site is an artificial reef affects fish community composition, accounting for 14.5% of variation in composition

#ONLY MACRO

#Start with root transformation (high abundance spp less pull in ordination)
dat_macroinvert_averages_bysite.wide.rt <- sqrt(dat_macroinvert_averages_bysite.wide[,c(5:ncol(dat_macroinvert_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_macroinvert_averages_bysite.wide.rt <- cbind(dat_macroinvert_averages_bysite.wide[,c(1:4), with = FALSE], dat_macroinvert_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_macroinvert_averages_bysite.wide.rt.trim <- dat_macroinvert_averages_bysite.wide.rt[,c(5:ncol(dat_macroinvert_averages_bysite.wide.rt)), with = FALSE]

permanova_macroinvert <- adonis2(
  dat_macroinvert_averages_bysite.wide.rt.trim ~ dat_macroinvert_averages_bysite.wide.rt $DepthZone,
  method = "bray"
)

permanova_macroinvert

#whether or not a site is an artificial reef affects macroinvert community composition, but only accounting for 10.5% of variation in composition


#ONLY KELP

#sometimes, no kelp at all, need to delete these rows
kelp_row_sums <- rowSums(dat_kelp_averages_bysite.wide[,5:ncol(dat_kelp_wide_density)])

dat_kelp_averages_bysite.wide[,rowSums := kelp_row_sums]

dat_kelp_averages_bysite.wide.r <- dat_kelp_averages_bysite.wide[rowSums>0]

#Start with root transformation (high abundance spp less pull in ordination)
dat_kelp_averages_bysite.wide.rt <- sqrt(dat_kelp_averages_bysite.wide.r[,c(5:ncol(dat_kelp_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_kelp_averages_bysite.wide.rt <- cbind(dat_kelp_averages_bysite.wide.r[,c(1:4), with = FALSE], dat_kelp_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_kelp_averages_bysite.wide.rt.trim <- dat_kelp_averages_bysite.wide.rt[,c(5:ncol(dat_kelp_averages_bysite.wide.rt)), with = FALSE]

permanova_kelp <- adonis2(
  dat_kelp_averages_bysite.wide.rt.trim ~ dat_kelp_averages_bysite.wide.rt $DepthZone,
  method = "bray"
)

permanova_kelp

#whether or not a site is an artificial reef affects fish community composition, accounting for 14.2% of variation in composition

#Take away:
#Artificial vs natural accounts for the most variation for fish and kelp communities, and the least relative variation for macroinvert communities


###################
#PLOT NMDS
###################
#NMDS “preserves the rank order of the inter-point dissimilarities as well as possible within the constraints of a small number of dimensions” (Anderson et al. 2008, p. 105).
#It has the least restrictive assumptions (any distance measure, no assumptions about shape of species response curves).

#all species

#convert to distance matrix
dat_averages_bysite.dist <- vegdist(dat_averages_bysite.wide[,5:ncol(dat_averages_bysite.wide)], distance = "bray")

full_nmds <- metaMDS(dat_averages_bysite.wide[,5:ncol(dat_averages_bysite.wide)], #metaMDS automatically does this for us
             autotransform = FALSE,
             distance = "bray",
             engine = "monoMDS",
             k = 4, #need at least twice as many sample units as dimensions, will impact results
             weakties = TRUE,#Set to true if you commonly have sample units that don't share species 
                              #TRUE allows observations to be positioned at different coordinates within the ordination space even though they are the same distance apart in the matrix
                              #FALSE forces observations that are the same distance apart in the distance matrix to also be the same distance apart in ordination space
             model = "global",
             maxit = 300,
             try = 40,
             trymax = 100)

print(full_nmds)


#The number of dimensions is often selected by conducting NMDS ordinations with different 
#numbers of dimensions (e.g., k = 1 to 5) and then plotting stress as a function of k. 
#The intent here is to identify the point where adding the complexity of an additional dimension contributes relatively little to the reduction in stress.

#Stress: 
#         k = 2, stress = .168
#         k = 3, stress = .112
#         k = 4, stress = 0.07, good ordination with no real risk of drawing false inferences (Clark 1993 pg 126)
#         k = 5, stress = 0.06
#         k = 6, stress = 0.05
#         k = 3, stress = 0.04

#The fit of a NMDS ordination can be assessed by plotting the original dissimilarities (z$diss) against the (Euclidean) ordination distances (z$dist) 

stressplot(full_nmds)

#link NMDS points to site details
dat_averages_bysite.nmds <- cbind(dat_fish_averages_bysite.wide[,1:4], full_nmds$points)

#set apart SoS
dat_averages_bysite.nmds[Site == "Star of Scotland",`Site type` := "Star of Scotland\nn = 1"]

#For some reason, Marina Del Rey coded as Region = Artificial reef, change to Santa Monica Bay

dat_averages_bysite.nmds[Region == "Artificial Reef",Region := "Santa Monica Bay"]

full_nmds_plot_bysitetype <- ggplot(data = dat_averages_bysite.nmds) +
  geom_point(aes(x = MDS1, y = MDS2, fill = `Site type`, shape = `Site type`), size = 3, color = "transparent") +
  scale_fill_manual(values = c("#BF6992", "#54AEAA", "#EA3323")) +
  scale_shape_manual(values = c(21,24,24)) +
  lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = c(.15, .9))


ggsave(full_nmds_plot_bysitetype, path = "figures", filename = "full_nmds_plot_bysitetype.jpg", height = 6, width = 6)

full_nmds_plot_byregion <- ggplot(data = dat_averages_bysite.nmds) +
  geom_point(aes(x = MDS1, y = MDS2, color = Region, shape = `Site type`), size = 3) +
  scale_color_manual(values = c("#83752D", "#EAC4AA", "#D38EB6", "#F6D798","#76A6B8", "#2F624E")) +
  scale_shape_manual(values = c(19,17,17), guide = NULL) +
  lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = c(.2, .83),
        legend.background = element_blank(), legend.key = element_blank())


ggsave(full_nmds_plot_byregion, path = "figures", filename = "full_nmds_plot_byregion.jpg", height = 6, width = 6)


#SPECIES SCORES
species_scores <- data.table(vegan::scores(full_nmds, "species"))
spp_list <- colnames(dat_averages_bysite.wide[,5:ncol(dat_averages_bysite.wide)])
species_scores[,Species := spp_list]

#top and bottom three for NMDS1 & 2
top_all_spp <- c("Stereolepis gigas","Pachycerianthus fimbriatus","Crossata californica",
                 "Paralabrax maculatofasciatus","Sebastes constellatus","Heterostichus rostratus",
                 "Lythrypnus dalli","Haliotis fulgens","Centrostephanus coronatus",
                 "Berthella californica","Pomaulax gibberosus","Acanthodoris rhodoceras")

species_scores_top <- species_scores[Species %in% top_all_spp]

species_scores_top[,Species := factor(Species,
                                      labels = 
                                        c("Stereolepis gigas\ngiant sea bass","Pachycerianthus fimbriatus\ntube anenome","Crossata californica\nfrogsnail",
                                                 "Paralabrax maculatofasciatus\nspotted sand bass","Sebastes constellatus\nstarry rockfish","Heterostichus rostratus\ngiant kelpfish",
                                                 "Lythrypnus dalli\nblue-banded goby","Haliotis fulgens\ngreen abalone","Centrostephanus coronatus\ncrowned sea urchin",
                                                 "Berthella californica\nsea slug spp.","Pomaulax gibberosus\nred turban snail","Acanthodoris rhodoceras\nsea slug spp."))]

#add to simple NMDS plot

full_nmds_plot_speciesscores <- ggplot(data = dat_averages_bysite.nmds) +
  geom_point(aes(x = MDS1, y = MDS2,shape = `Site type`), size = 3, alpha = 0.3) +
  geom_text_repel(data = species_scores_top, aes(x = NMDS1, y = NMDS2, label = Species), size = 2) +
  scale_shape_manual(values = c(19,17,17)) +
  lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = c(.15, .9))


ggsave(full_nmds_plot_speciesscores, path = "figures", filename = "full_nmds_plot_speciesscores.jpg", height = 5, width = 10)

#plot together
full_nmds_plot <- plot_grid(full_nmds_plot_bysitetype, full_nmds_plot_byregion,full_nmds_plot_speciesscores + theme(legend.position = "none"), ncol = 3)

ggsave(full_nmds_plot, path = "figures", filename = "full_nmds_plot.jpg", height = 5, width = 15)
ggsave(full_nmds_plot, path = "figures", filename = "full_nmds_plot.pdf", height = 5, width = 15)

#We can create ellipses around groups

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(factor(dat_averages_bysite.nmds$DepthZone))){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(dat_averages_bysite.nmds[dat_averages_bysite.nmds$DepthZone==g,],
                                                   veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                                ,group=g))
}

full_nmds_plot_bysitetype +
geom_path(data=df_ell, aes(x=MDS1, y=MDS2,color=group), size=1, linetype=2) +
  scale_color_manual(values = c("#BF6992", "#54AEAA"), guide = NULL)

#NOT IMMEDIATELY SURE WHAT THIS ELLIPSE IS TELLING US



######################
#To try next
##repeat for fish, fish biomass, macro, kelp, and all species clumped together (maybe do this one by presence absence only)
##Integrate environmental data (distance to 200m isobath, insitu: depth, temp, rugosity or whatever)





