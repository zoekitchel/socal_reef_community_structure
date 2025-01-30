# CREATION DATE 28 Jan 2024
# MODIFICATION DATE 13 Dec 2024

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Community visualizations and analyses

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
dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))
dat_fish_site_averages <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages.rds"))
dat_macroinvert_site_averages <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages.rds"))
dat_kelp_site_averages <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages.rds"))

#Small data adjustments
  #included giant kelp stipes for in situ habitat data, but DELETE for community analyses
  dat_kelp_site_averages <- dat_kelp_site_averages[BenthicReefSpecies != "Macrocystis pyrifera stipes",] 
  
  #pull in spp taxonomy info
  species_key <- fread(file.path("keys","species_key.csv"))
  
  #capitalize California
  species_key[, common_name_final := gsub("california sheephead", "California sheephead", common_name_final)]
  
  #find and replace Semicossyphus pulcher with Bodianus pulcher in VRG data
  dat_fish_site_averages[, Species := gsub("Semicossyphus pulcher", "Bodianus pulcher", Species)]
  
  #link site averaged data with species key
  dat_fish_site_averages <- species_key[dat_fish_site_averages, on = c("taxa" = "Species")]
  dat_macroinvert_site_averages <- species_key[dat_macroinvert_site_averages, on = c("taxa" = "BenthicReefSpecies")]
  dat_kelp_site_averages <- species_key[dat_kelp_site_averages, on = c("taxa" = "BenthicReefSpecies")]
  
  #New Column Identifying ARM vs Island vs Natural Coast
  dat_fish_site_averages[,type := ifelse(DepthZone %in% c("ARM","Module"),"ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]
  dat_macroinvert_site_averages[,type := ifelse(DepthZone %in% c("ARM","Module"),"ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]
  dat_kelp_site_averages[,type := ifelse(DepthZone %in% c("ARM","Module"),"ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]
  
  #Edit so that Depth Zone = Module -> Depth Zone ARM
  dat_fish_site_averages[,DepthZone := ifelse(DepthZone %in% c("ARM","Module"),"ARM",as.character(DepthZone))]
  dat_macroinvert_site_averages[,DepthZone := ifelse(DepthZone %in% c("ARM","Module"),"ARM",as.character(DepthZone))]
  dat_kelp_site_averages[,DepthZone := ifelse(DepthZone %in% c("ARM","Module"),"ARM",as.character(DepthZone))]
  
  
##################################################
#Manipulating data, extract richness per site and depth zone
##################################################

#In rare cases when more than one observation per site per year, take average, and then reduce to required columns
#Fish density
dat_fish_long_density <- dat_fish_site_averages[,.(taxa,Region,Site,DepthZone, mean_density_m2)]
dat_fish_long_density[,mean_density_m2 := mean(mean_density_m2),.(taxa,Region,Site,DepthZone)]
dat_fish_long_density <- unique(dat_fish_long_density[,.(taxa,Region,Site,DepthZone,mean_density_m2)])

#Fish biomass
dat_fish_long_biomass <- dat_fish_site_averages[,.(taxa,Region,Site,DepthZone, mean_wt_density_g_m2)]
dat_fish_long_biomass[,mean_wt_density_g_m2 := mean(mean_wt_density_g_m2),.(taxa,Region,Site,DepthZone)]
dat_fish_long_biomass <- unique(dat_fish_long_biomass[,.(taxa,Region,Site,DepthZone,mean_wt_density_g_m2)])

#Macroinvert density
dat_macroinvert_long_density <- dat_macroinvert_site_averages[,.(taxa,Region,Site,DepthZone, mean_density_m2)]
dat_macroinvert_long_density[,mean_density_m2 := mean(mean_density_m2),.(taxa,Region,Site,DepthZone)]
dat_macroinvert_long_density <- unique(dat_macroinvert_long_density[,.(taxa,Region,Site,DepthZone,mean_density_m2)])

#Kelp density
dat_kelp_long_density <- dat_kelp_site_averages[,.(taxa,Region,Site,DepthZone, mean_density_m2)]
dat_kelp_long_density[,mean_density_m2 := mean(mean_density_m2),.(taxa,Region,Site,DepthZone)]
dat_kelp_long_density <- unique(dat_kelp_long_density[,.(taxa,Region,Site,DepthZone,mean_density_m2)])

#Remove species with 0 observations
dat_fish_long_density_removezeros <- dat_fish_long_density[mean_density_m2>0,]
dat_fish_long_biomass_removezeros <- dat_fish_long_biomass[mean_wt_density_g_m2>0,]
dat_macroinvert_long_density_removezeros <- dat_macroinvert_long_density[mean_density_m2>0,]
dat_kelp_long_density_removezeros <- dat_kelp_long_density[mean_density_m2>0,]

#Number of taxa per depth zone?
fish_depthzone_site_richness_abun <- dat_fish_long_density_removezeros[,.(count_spp = .N,sum_abun = sum(mean_density_m2)),.(Site, DepthZone, Region)]
fish_depthzone_site_richness_abun[,category:="fish"]
fish_depthzone_site_biomass <- dat_fish_long_biomass_removezeros[,.(sum_biomass = sum(mean_wt_density_g_m2)),.(Site, DepthZone, Region)]
fish_depthzone_site_richness_abun_biomass <- fish_depthzone_site_richness_abun[fish_depthzone_site_biomass, on = c("Site","DepthZone","Region")]

macroinvert_depthzone_site_richness <- dat_macroinvert_long_density_removezeros[,.(count_spp = .N, sum_abun = sum(mean_density_m2)),.(Site, DepthZone, Region)]
macroinvert_depthzone_site_richness[,category:="macroinvert"]

kelp_depthzone_site_richness <- dat_kelp_long_density_removezeros[,.(count_spp = .N, sum_abun = sum(mean_density_m2)),.(Site, DepthZone, Region)]
kelp_depthzone_site_richness[,category:="kelp"]

#Merge across taxa
depthzone_site_richness <- rbind(fish_depthzone_site_richness_abun_biomass, macroinvert_depthzone_site_richness, kelp_depthzone_site_richness, fill = TRUE)

##################################################################################################
#Long data to wide data for multivariate community visualizations and analyses (without zeros!)
##################################################################################################

#melt long to wide
dat_fish_averages_bysite.wide <- dcast(dat_fish_long_density_removezeros, Region + Site + DepthZone ~ taxa, value.var = "mean_density_m2", fun = mean, fill = 0)

dat_fish_biomass_averages_bysite.wide <- dcast(dat_fish_long_biomass_removezeros, Region + Site + DepthZone ~ taxa, value.var = "mean_wt_density_g_m2", fun = mean, fill = 0)

dat_macroinvert_averages_bysite.wide <- dcast(dat_macroinvert_long_density_removezeros, Region + Site + DepthZone ~ taxa, value.var = "mean_density_m2", fun = mean, fill = 0)

dat_kelp_averages_bysite.wide <- dcast(dat_kelp_long_density_removezeros, Region + Site + DepthZone ~ taxa, value.var = "mean_density_m2", fun = mean, fill = 0)

#Merge for all density data tables for species visualization (biomass of fish kept separate)
dat_averages_bysite.wide <- dat_fish_averages_bysite.wide[dat_macroinvert_averages_bysite.wide, on = c("Region","Site","DepthZone")]
dat_averages_bysite.wide <- dat_kelp_averages_bysite.wide[dat_averages_bysite.wide, on = c("Region","Site","DepthZone")]

#there are some sites that have no kelp, and these come up as NAs, change to 0 as vegan can't deal
dat_averages_bysite.wide[is.na(dat_averages_bysite.wide)] <- 0


##################################################################################################
#Visualize full PcoA
##################################################################################################
MPA_site_key <- readRDS(file.path("keys","MPA_site_key.rds"))

#CHECK FOR NUMBERS IN 165
stopifnot(colnames(dat_averages_bysite.wide[,c(1:3,24,36)]) == c("Region","Site","DepthZone","Alloclinus holderi","Embiotoca lateralis"))
          
#species only
dat_averages_bysite.wide.spp <- dat_averages_bysite.wide[,c(4:199)]

#Using square root for PERMANOVA, so I should use square root here too
dat_averages_bysite.wide.spp.l <- sqrt(dat_averages_bysite.wide.spp)

#transformed but with attributes attached
dat_averages_bysite.wide.spp.attributes <- cbind(dat_averages_bysite.wide[,1:3], dat_averages_bysite.wide.spp.l)
dat_averages_bysite.wide.spp.attributes <- dat_averages_bysite.wide.spp.attributes[,DepthZone_wAR := ifelse(DepthZone == "ARM"& grepl("PVR",Site) == T,"AR_PVR",ifelse(DepthZone == "ARM" & grepl("PVR",Site) == F, "AR_SM",as.character(DepthZone)))]
dat_averages_bysite.wide.spp.attributes <- MPA_site_key[dat_averages_bysite.wide.spp.attributes, on = "Site"]
dat_averages_bysite.wide.spp.attributes <- dat_averages_bysite.wide.spp.attributes[,type := ifelse(DepthZone %in% c("ARM","Module"),"ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]

#distance matrix (distance = bray curtis)
dat_averages_bysite.dist <- vegdist(dat_averages_bysite.wide.spp.l, distance = "bray")

#PCoA
allspp_PCoA <- wcmdscale(dat_averages_bysite.dist, eig = TRUE)

allspp_PCoA.pnts <- data.table(allspp_PCoA$points)

#Add group variables
PCoA_sqrt_grouping <- cbind(dat_averages_bysite.wide[,c(1:3)],allspp_PCoA.pnts[,1:2])

#New column splitting SMB and PVR ARs
PCoA_sqrt_grouping[,DepthZone_wAR := ifelse(DepthZone == "ARM"& grepl("PVR",Site) == T,"AR_PVR",ifelse(DepthZone == "ARM" & grepl("PVR",Site) == F, "AR_SM",as.character(DepthZone)))]

#Adjust factor order
PCoA_sqrt_grouping[,DepthZone := factor(DepthZone,
                                       levels = c("Inner","Middle","Outer","Deep","ARM"),
                                       labels = c("Inner","Middle","Outer","Deep","AR"))]

PCoA_sqrt_grouping[,DepthZone_wAR := factor(DepthZone_wAR,
                                   levels = c("Inner","Middle","Outer","Deep","AR_SM","AR_PVR"),
                                   labels = c("Inner","Middle","Outer","Deep","Santa Monica Bay","Palos Verdes"))]

#Add mainland versus island designation
PCoA_sqrt_grouping[,type := factor(ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")))]


#Visualize without ellipses
PCoA_allspp_allsite_points <- ggplot(PCoA_sqrt_grouping) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100","black"), guide = guide_legend(reverse = TRUE)) +
  geom_point(aes(Dim1, Dim2, color = DepthZone, shape = DepthZone), size = 4, fill = "#FE6100") +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100","black"), guide = guide_legend(reverse = TRUE)) +
  # lims(x = c(-0.5,0.45),y = c(-0.5,0.35)) +
  scale_shape_manual(values = c(15,17,19,23,7), guide = guide_legend(reverse = TRUE)) +
  theme_classic()

#Natural sites colored by depth zone
PCoA_allspp_allsite_points_natural_only <- ggplot(PCoA_sqrt_grouping[DepthZone != "AR"]) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(aes(Dim1, Dim2, color = DepthZone, shape = DepthZone), size = 4, fill = "#FE6100") +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  lims(x = c(-0.50,0.6),y = c(-0.5,0.75)) +
  labs(x = paste0("PCoA dim 1 (",round(allspp_PCoA$eig[1],1),"%)"),y = paste0("PCoA dim 2 (",round(allspp_PCoA$eig[2],1),"%)"), shape = "Depth zone",color = "Depth zone") +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  annotate(geom = "text",label ="Square-root transformation\nDistance: Bray-Curtis dissimilarity", 
           x = 0.55, y = -0.5, size = 4, hjust = 1) +
  theme_classic() +
  theme(legend.position = c(0.5,0.95), legend.direction = "horizontal",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    color = guide_legend(override.aes = list(size = 6), title.position = "top", title.hjust = 0.5, reverse = T),
    shape = guide_legend(override.aes = list(size = 6), title.position = "top", title.hjust = 0.5, reverse = T),
  )

ggsave(PCoA_allspp_allsite_points_natural_only, path = "figures", filename = "PCoA_allspp_allsite_points_natural_only.jpg", height = 8, width = 12, unit = "in", dpi = 300)


#Mainland vs island
PCoA_allspp_allsite_points_natural_only_island_mainland <- ggplot(PCoA_sqrt_grouping[DepthZone != "AR"]) +
  geom_point(aes(Dim1, Dim2, fill = type), shape = 24, color = "black", size = 3) +
  scale_fill_manual(values = c("black","white")) +
  lims(x = c(-0.50,0.6),y = c(-0.5,0.75)) +
  labs(x = "PCoA dim 1",y = "PCoA dim 2", fill = "Reef type") +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(legend.direction = "horizontal",
        legend.position = c(.55, .95),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    fill = guide_legend(override.aes = list(size = 6), title.position = "top",title.hjust = 0.5)
  )

ggsave(PCoA_allspp_allsite_points_natural_only_island_mainland, path = "figures", filename = "PCoA_allspp_allsite_points_natural_only_island_mainland.jpg", height = 8, width = 12, unit = "in", dpi = 300)


##just natural sites
PCoA_allspp_allsite_points_natural_only_highlight_mainland <- ggplot(PCoA_sqrt_grouping[DepthZone != "AR"]) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(aes(Dim1, Dim2, color = DepthZone, shape = DepthZone, alpha = type), size = 4, fill = "#FE6100") +
  labs(x = "PCoA Dim 1",y = "PCoA Dim2") +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  scale_alpha_manual(values = c(0.3,1), guide = "none") +
  lims(x = c(-0.50,0.6),y = c(-0.5,0.75)) +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  labs(x = "PCoA Dim 1",y = "PCoA Dim2") +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    color = guide_legend(override.aes = list(size = 6)),
    shape = guide_legend(override.aes = list(size = 6)),
  )

ggsave(PCoA_allspp_allsite_points_natural_only_highlight_mainland, path = "figures", filename = "PCoA_allspp_allsite_points_natural_only_highlight_mainland.jpg", height = 8, width = 12, unit = "in", dpi = 300)


#
PCoA_allspp_allsite_points_natural_only_highlight_island <- ggplot(PCoA_sqrt_grouping[DepthZone != "AR"]) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(aes(Dim1, Dim2, color = DepthZone, shape = DepthZone, alpha = type), size = 4, fill = "#FE6100") +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  scale_alpha_manual(values = c(1,0.3), guide = "none") +
  lims(x = c(-0.50,0.6),y = c(-0.5,0.75)) +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  theme_classic()+
  labs(x = "PCoA Dim 1",y = "PCoA Dim2") +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    color = guide_legend(override.aes = list(size = 6)),
    shape = guide_legend(override.aes = list(size = 6))
  )

ggsave(PCoA_allspp_allsite_points_natural_only_highlight_island, path = "figures", filename = "PCoA_allspp_allsite_points_natural_only_highlight_island.jpg", height = 8, width = 12, unit = "in", dpi = 300)

#mark if it is PVR
PCoA_sqrt_grouping[,PVR := ifelse(grepl("PVR",Site)==T,TRUE,FALSE)]

#Add in Old and New AR reefs
PCoA_allspp_allsite_points_withAR <- ggplot() +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(data = PCoA_sqrt_grouping[DepthZone != "AR"], aes(Dim1, Dim2, color = DepthZone, shape = DepthZone), alpha = 0.6, size = 4, fill = "#FE6100") +
  geom_point(data = PCoA_sqrt_grouping[DepthZone == "AR" & PVR == TRUE], aes(Dim1, Dim2), shape = 12, color = "black", size = 5) + #PVR only
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  lims(x = c(-0.50,0.6),y = c(-0.5,0.75)) +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  theme_classic()+
  labs(x = paste0("PCoA Dim 1 (",round(allspp_PCoA$eig[1],1)," %)"),y = paste0("PCoA Dim 2 (",round(allspp_PCoA$eig[2],1)," %)")) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    color = guide_legend(override.aes = list(size = 6)),
    shape = guide_legend(override.aes = list(size = 6))
  )

ggsave(PCoA_allspp_allsite_points_withAR, path = "figures", filename = "PCoA_allspp_allsite_points_withAR.jpg", height = 8, width = 12, unit = "in", dpi = 300)

PCoA_allspp_allsite_points_withallAR <- PCoA_allspp_allsite_points_withAR +
  geom_point(data = PCoA_sqrt_grouping[DepthZone == "AR" & PVR == FALSE], aes(Dim1, Dim2), shape = 8, color = "black", size = 5) #ADD SMB
  
ggsave(PCoA_allspp_allsite_points_withallAR, path = "figures", filename = "PCoA_allspp_allsite_points_withallAR.jpg", height = 8, width = 12, unit = "in", dpi = 300)


#visualize
PCoA_allspp_allsite <- ggplot(PCoA_sqrt_grouping) +
  stat_ellipse(geom = "polygon", aes(Dim1, Dim2, fill = DepthZone), color = NA,alpha = 0.4) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100","black"), guide = guide_legend(reverse = TRUE)) +
  geom_point(aes(Dim1, Dim2, color = DepthZone, shape = DepthZone), size = 4, fill = "#FE6100") +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100","black"), guide = guide_legend(reverse = TRUE)) +
 # lims(x = c(-0.5,0.45),y = c(-0.5,0.35)) +
  scale_shape_manual(values = c(15,17,19,23,7), guide = guide_legend(reverse = TRUE)) +
  theme_classic()

#Central point for each
PCoA_sqrt_centroid <- PCoA_sqrt_grouping[, lapply(.SD, mean, na.rm=TRUE), by=c("DepthZone"), .SDcols = c("Dim1","Dim2")] 

PCoA_allspp_centroids <- ggplot() +
  stat_ellipse(data = PCoA_sqrt_grouping[DepthZone != "AR"], geom = "polygon", aes(Dim1, Dim2, fill = DepthZone), color = "white",alpha = 0.2) +
  stat_ellipse(data = PCoA_sqrt_grouping[DepthZone != "AR"], geom = "polygon", aes(Dim1, Dim2,linetype = type), color = "gray32",fill = NA) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(data = PCoA_sqrt_grouping[DepthZone == "AR"], aes(Dim1, Dim2, shape = DepthZone_wAR),color = "black", size = 4) +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values = c(9,10)) +
  annotate(geom = "text",label = "D", x = 0.4, y = 0.25, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "O", x = 0.05, y = 0.39, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "M", x = -0.06, y = 0.38, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "I", x = -0.31, y = 0.25, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "Mainland", x = -0.3, y = -0.27, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "Island", x = -0.4, y = 0.5, size = 5, fontface = 'bold') +
  guides(fill = "none",linetype = "none",shape = guide_legend(title.position = "top",title.hjust = 0.5)) +
  lims(x = c(-0.50,0.6),y = c(-0.5,0.75)) +
  labs(x = "PCoA dim 1",y = "PCoA dim 2", shape = "Artificial reef",color = "Depth zone", fill = "Depth zone", linetype = "Reef type") +
  theme_classic() +
  theme(legend.justification = "center",
        legend.direction = "horizontal",
        legend.position = c(0.5,0.95),
        legend.text = element_text(size = 14),   # Increase legend text size
        legend.key.size = unit(1.5, "lines"),
        legend.background = element_blank(),
              legend.title = element_text(size = 16),
              axis.text = element_text(size = 16),
              axis.title = element_text(size = 16))     # Increase legend key size

#plot natural reefs next to artificial reef highlight plot
PCoA_allspp_allsite_AR_natural_merge <- plot_grid(PCoA_allspp_allsite_points_natural_only,
                                                  PCoA_allspp_allsite_points_natural_only_island_mainland + theme(axis.title.y = element_blank(), axis.text.y = element_blank()),
                                                  PCoA_allspp_centroids+theme(axis.title.y = element_blank(), axis.text.y = element_blank()), ncol = 3, align = "hv", labels = c("a.","b.","c."), label_size = 20)


ggsave(PCoA_allspp_allsite_AR_natural_merge, path = "figures", filename = "PCoA_allspp_allsite_AR_natural_merge.jpg", width =18, height = 6)
  



##################################################################################################
#Visualize PcoA for each individual taxonomic group ####
##################################################################################################
#FISH DENSITY ####
#CHECK FOR NUMBERS IN 165
stopifnot(colnames(dat_fish_averages_bysite.wide[,c(1:4)]) == c("Region","Site","DepthZone","Alloclinus holderi"))

#fish species only
dat_fish_averages_bysite.wide.spp <- dat_fish_averages_bysite.wide[,c(5:74)]

#Using square root for permanova, so I should use square root here too
dat_fish_averages_bysite.wide.spp.l <- sqrt(dat_fish_averages_bysite.wide.spp)

#transformed but with depth zone attached
dat_fish_averages_bysite.wide.spp.depthzone <- cbind(dat_fish_averages_bysite.wide[,3], dat_fish_averages_bysite.wide.spp.l)

#distance matrix (distance = bray curtis)
dat_fish_averages_bysite.wide.spp.dist <- vegdist(dat_fish_averages_bysite.wide.spp.l, distance = "bray")

#PCoA
fish_PCoA <- wcmdscale(dat_fish_averages_bysite.wide.spp.dist, eig = TRUE)

fish_PCoA.pnts <- data.table(fish_PCoA$points)

#add group variables
PCoA_fish_sqrt_grouping <- cbind(dat_fish_averages_bysite.wide[,c(1:3)],fish_PCoA.pnts[,1:2])

#New column splitting SMB and PVR ARs
PCoA_fish_sqrt_grouping[,DepthZone_wAR := ifelse(DepthZone == "ARM"& grepl("PVR",Site) == T,"AR_PVR",ifelse(DepthZone == "ARM" & grepl("PVR",Site) == F, "AR_SM",as.character(DepthZone)))]

#adjust factor order
PCoA_fish_sqrt_grouping[,DepthZone := factor(DepthZone,
                                   levels = c("Inner","Middle","Outer","Deep","ARM"),
                                   labels = c("Inner","Middle","Outer","Deep","AR"))]

PCoA_fish_sqrt_grouping[,DepthZone_wAR := factor(DepthZone_wAR,
                                       levels = c("Inner","Middle","Outer","Deep","AR_SM","AR_PVR"),
                                       labels = c("Inner","Middle","Outer","Deep","Santa Monica Bay","Palos Verdes"))]

#add mainland versus island designation
PCoA_fish_sqrt_grouping[,type := factor(ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")))]

#visualize
#just natural sites
PCoA_fish_allsite_points_natural_only <- ggplot(PCoA_fish_sqrt_grouping[DepthZone != "AR"]) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(aes(Dim1, Dim2, color = DepthZone, shape = DepthZone), size = 4, fill = "#FE6100") +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  lims(x = c(-0.5,0.6),y = c(-0.65,0.5)) +
  labs(x = paste0("PCoA dim 1 (",round(fish_PCoA$eig[1],1),"%)"),y = paste0("Fish\nPCoA dim 2 (",round(fish_PCoA$eig[2],1),"%)"), shape = "Depth zone",color = "Depth zone") +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  annotate(geom = "text",label ="Square-root transformation\nDistance: Bray-Curtis dissimilarity", 
           x = 0.6, y = -0.6, size = 4, hjust = 1) +
  theme_classic() +
  theme(legend.position = c(0.5,0.95), legend.direction = "horizontal",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    color = guide_legend(override.aes = list(size = 6), title.position = "top", title.hjust = 0.5, reverse = T),
    shape = guide_legend(override.aes = list(size = 6), title.position = "top", title.hjust = 0.5, reverse = T),
  )

ggsave(PCoA_fish_allsite_points_natural_only, path = "figures", filename = "PCoA_fish_allsite_points_natural_only.jpg", height = 8, width = 12, unit = "in", dpi = 300)


#mainland vs island
#just natural sites
PCoA_fish_allsite_points_natural_only_island_mainland <- ggplot(PCoA_fish_sqrt_grouping[DepthZone != "AR"]) +
  geom_point(aes(Dim1, Dim2, fill = type), shape = 24, color = "black", size = 3) +
  scale_fill_manual(values = c("black","white")) +
  lims(x = c(-0.5,0.6),y = c(-0.65,0.5)) +
  labs(x = "PCoA dim 1",y = "PCoA dim 2", fill = "Reef type") +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(legend.direction = "horizontal",
        legend.position = c(.55, .95),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    fill = guide_legend(override.aes = list(size = 6), title.position = "top",title.hjust = 0.5)
  )

ggsave(PCoA_fish_allsite_points_natural_only_island_mainland, path = "figures", filename = "PCoA_fish_allsite_points_natural_only_island_mainland.jpg", height = 8, width = 12, unit = "in", dpi = 300)


#merge plot
PCoA_fish_centroids <- ggplot() +
  stat_ellipse(data = PCoA_fish_sqrt_grouping[DepthZone != "AR"], geom = "polygon", aes(Dim1, Dim2, fill = DepthZone), color = "white",alpha = 0.2) +
  stat_ellipse(data = PCoA_fish_sqrt_grouping[DepthZone != "AR"], geom = "polygon", aes(Dim1, Dim2,linetype = type), color = "gray32",fill = NA) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(data = PCoA_fish_sqrt_grouping[DepthZone == "AR"], aes(Dim1, Dim2, shape = DepthZone_wAR),color = "black", size = 4) +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values = c(9,10)) +
  annotate(geom = "text",label = "I", x = -0.45, y = -0.08, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "M", x = -0.4, y = -0.15, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "O", x = -0.29, y = -0.25, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "D", x = -0.1, y = -0.5, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "Island", x = 0.3, y = 0.3, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "Mainland", x = -0.15, y = -0.22, size = 5, fontface = 'bold') +
  guides(fill = "none",linetype = "none",shape = guide_legend(title.position = "top",title.hjust = 0.5)) +
  lims(x = c(-0.5,0.6),y = c(-0.65,0.5)) +
  labs(x = "PCoA dim 1",y = "PCoA dim 2", shape = "Depth zone",color = "Depth zone", fill = "Depth zone", linetype = "Reef type") +
  theme_classic() +
  theme(legend.justification = "center",legend.position = c(0.5,0.95), legend.direction = "horizontal",
        legend.text = element_text(size = 14),   # Increase legend text size
        legend.key.size = unit(1.5, "lines"),
        legend.title = element_text(size = 16),
        legend.background = element_blank(),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))     # Increase legend key size

#plot natural reefs next to artificial reef highlight plot
PCoA_fish_allsite_AR_natural_merge <- plot_grid(PCoA_fish_allsite_points_natural_only,
                                                  PCoA_fish_allsite_points_natural_only_island_mainland + theme(axis.title.y = element_blank(), axis.text.y = element_blank()),
                                                  PCoA_fish_centroids+theme(axis.title.y = element_blank(), axis.text.y = element_blank()), ncol = 3, align = "hv", labels = c("a.","b.","c."), label_size = 20)


ggsave(PCoA_fish_allsite_AR_natural_merge, path = "figures", filename = "PCoA_fish_allsite_AR_natural_merge.jpg", width =18, height = 6)

#MACROALGAE DENSITY ####
#CHECK FOR NUMBERS IN 165
stopifnot(colnames(dat_kelp_averages_bysite.wide[,c(1:4)]) == c("Region","Site","DepthZone","Agarum fimbriatum"))

#kelp species only
dat_kelp_averages_bysite.wide.spp <- dat_kelp_averages_bysite.wide[,c(4:23)]

#Using square root for permanova, so I should use square root here too
dat_kelp_averages_bysite.wide.spp.l <- sqrt(dat_kelp_averages_bysite.wide.spp)

#transformed but with depth zone attached
dat_kelp_averages_bysite.wide.spp.depthzone <- cbind(dat_kelp_averages_bysite.wide[,3], dat_kelp_averages_bysite.wide.spp.l)

#distance matrix (distance = bray curtis)
dat_kelp_averages_bysite.wide.spp.dist <- vegdist(dat_kelp_averages_bysite.wide.spp.l, distance = "bray")

#PCoA
kelp_PCoA <- wcmdscale(dat_kelp_averages_bysite.wide.spp.dist, eig = TRUE)

kelp_PCoA.pnts <- data.table(kelp_PCoA$points)

#add env and group variables
PCoA_kelp_sqrt_grouping <- cbind(dat_kelp_averages_bysite.wide[,c(1:3)],kelp_PCoA.pnts[,1:2])

#New column splitting SMB and PVR ARs
PCoA_kelp_sqrt_grouping[,DepthZone_wAR := ifelse(DepthZone == "ARM"& grepl("PVR",Site) == T,"AR_PVR",ifelse(DepthZone == "ARM" & grepl("PVR",Site) == F, "AR_SM",as.character(DepthZone)))]

#adjust factor order
PCoA_kelp_sqrt_grouping[,DepthZone := factor(DepthZone,
                                        levels = c("Inner","Middle","Outer","Deep","ARM"),
                                        labels = c("Inner","Middle","Outer","Deep","AR"))]

PCoA_kelp_sqrt_grouping[,DepthZone_wAR := factor(DepthZone_wAR,
                                            levels = c("Inner","Middle","Outer","Deep","AR_SM","AR_PVR"),
                                            labels = c("Inner","Middle","Outer","Deep","Santa Monica Bay","Palos Verdes"))]


#add mainland versus island designation
PCoA_kelp_sqrt_grouping[,type := factor(ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")))]

#visualize
#just natural sites
PCoA_kelp_allsite_points_natural_only <- ggplot(PCoA_kelp_sqrt_grouping[DepthZone != "AR"]) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(aes(Dim1, Dim2, color = DepthZone, shape = DepthZone), size = 4, fill = "#FE6100") +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  lims(x = c(-0.7,0.6),y = c(-0.45,0.65)) +
  labs(x = paste0("PCoA dim 1 (",round(kelp_PCoA$eig[1],1),"%)"),y = paste0("Macroalgae\nPCoA dim 2 (",round(kelp_PCoA$eig[2],1),"%)"), shape = "Depth zone",color = "Depth zone") +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  annotate(geom = "text",label ="Square-root transformation\nDistance: Bray-Curtis dissimilarity", 
           x = 0.6, y = -0.45, size = 4, hjust = 1) +
  theme_classic() +
  theme(legend.position = c(0.5,0.95), legend.direction = "horizontal",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    color = guide_legend(override.aes = list(size = 6), title.position = "top", title.hjust = 0.5, reverse = T),
    shape = guide_legend(override.aes = list(size = 6), title.position = "top", title.hjust = 0.5, reverse = T),
  )

ggsave(PCoA_kelp_allsite_points_natural_only, path = "figures", filename = "PCoA_kelp_allsite_points_natural_only.jpg", height = 8, width = 12, unit = "in", dpi = 300)


#mainland vs island
#just natural sites
PCoA_kelp_allsite_points_natural_only_island_mainland <- ggplot(PCoA_kelp_sqrt_grouping[DepthZone != "AR"]) +
  geom_point(aes(Dim1, Dim2, fill = type), shape = 24, color = "black", size = 3) +
  scale_fill_manual(values = c("black","white")) +
  lims(x = c(-0.7,0.6),y = c(-0.45,0.65)) +
  labs(x = "PCoA dim 1",y = "PCoA dim 2", fill = "Reef type") +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(legend.direction = "horizontal",
        legend.position = c(.55, .95),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    fill = guide_legend(override.aes = list(size = 6), title.position = "top",title.hjust = 0.5)
  )

ggsave(PCoA_kelp_allsite_points_natural_only_island_mainland, path = "figures", filename = "PCoA_kelp_allsite_points_natural_only_island_mainland.jpg", height = 8, width = 12, unit = "in", dpi = 300)


#merge plot
PCoA_kelp_centroids <- ggplot() +
  stat_ellipse(data = PCoA_kelp_sqrt_grouping[DepthZone != "AR"], geom = "polygon", aes(Dim1, Dim2, fill = DepthZone), color = "white",alpha = 0.2) +
  stat_ellipse(data = PCoA_kelp_sqrt_grouping[DepthZone != "AR"], geom = "polygon", aes(Dim1, Dim2,linetype = type), color = "gray32",fill = NA) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(data = PCoA_kelp_sqrt_grouping[DepthZone == "AR"], aes(Dim1, Dim2, shape = DepthZone_wAR),color = "black", size = 4) +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values = c(9,10)) +
  annotate(geom = "text",label = "I", x = 0.17, y = -0.44, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "M", x = 0.39, y = -0.35, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "O", x = 0.5, y = -0.05, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "D", x = 0.5, y = 0.26, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "Mainland", x = 0.1, y = 0, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "Island", x = -0.45, y = 0.27, size = 5, fontface = 'bold') +
  lims(x = c(-0.7,0.6),y = c(-0.45,0.65)) +
  guides(fill = "none",linetype = "none",shape = guide_legend(title.position = "top",title.hjust = 0.5)) +
  labs(x = "PCoA dim 1",y = "PCoA dim 2", shape = "Depth zone",color = "Depth zone", fill = "Depth zone", linetype = "Reef type") +
  theme_classic() +
  theme(legend.position = c(0.5,0.95), legend.justification = "center",legend.background = element_blank(),legend.direction = "horizontal",
        legend.text = element_text(size = 14),   # Increase legend text size
        legend.key.size = unit(1.5, "lines"),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))     # Increase legend key size

#plot natural reefs next to artificial reef highlight plot
PCoA_kelp_allsite_AR_natural_merge <- plot_grid(PCoA_kelp_allsite_points_natural_only + theme(legend.position = "null"),
                                                PCoA_kelp_allsite_points_natural_only_island_mainland + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "null"),
                                                PCoA_kelp_centroids+theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "null"), ncol = 3, align = "hv", labels = c("d.","e.","f."), label_size = 20)


ggsave(PCoA_kelp_allsite_AR_natural_merge, path = "figures", filename = "PCoA_kelp_allsite_AR_natural_merge.jpg", width =18, height = 6)


#MACROINVERT DENSITY#CHECK FOR NUMBERS IN 165
stopifnot(colnames(dat_macroinvert_averages_bysite.wide[,c(1:4)]) == c("Region","Site","DepthZone","Acanthodoris lutea"))

#macroinvert species only
dat_macroinvert_averages_bysite.wide.spp <- dat_macroinvert_averages_bysite.wide[,c(4:106)]

#Using square root for permanova, so I should use square root here too
dat_macroinvert_averages_bysite.wide.spp.l <- sqrt(dat_macroinvert_averages_bysite.wide.spp)

#transformed but with depth zone attached
dat_macroinvert_averages_bysite.wide.spp.depthzone <- cbind(dat_macroinvert_averages_bysite.wide[,3], dat_macroinvert_averages_bysite.wide.spp.l)

#distance matrix (distance = bray curtis)
dat_macroinvert_averages_bysite.wide.spp.dist <- vegdist(dat_macroinvert_averages_bysite.wide.spp.l, distance = "bray")

#PCoA
macroinvert_PCoA <- wcmdscale(dat_macroinvert_averages_bysite.wide.spp.dist, eig = TRUE)

macroinvert_PCoA.pnts <- data.table(macroinvert_PCoA$points)

#add env and group variables
PCoA_macroinvert_sqrt_grouping <- cbind(dat_macroinvert_averages_bysite.wide[,c(1:3)],macroinvert_PCoA.pnts[,1:2])

#New column splitting SMB and PVR ARs
PCoA_macroinvert_sqrt_grouping[,DepthZone_wAR := ifelse(DepthZone == "ARM"& grepl("PVR",Site) == T,"AR_PVR",ifelse(DepthZone == "ARM" & grepl("PVR",Site) == F, "AR_SM",as.character(DepthZone)))]

#adjust factor order
PCoA_macroinvert_sqrt_grouping[,DepthZone := factor(DepthZone,
                                        levels = c("Inner","Middle","Outer","Deep","ARM"),
                                        labels = c("Inner","Middle","Outer","Deep","AR"))]

PCoA_macroinvert_sqrt_grouping[,DepthZone_wAR := factor(DepthZone_wAR,
                                            levels = c("Inner","Middle","Outer","Deep","AR_SM","AR_PVR"),
                                            labels = c("Inner","Middle","Outer","Deep","Santa Monica Bay","Palos Verdes"))]


#add mainland versus island designation
PCoA_macroinvert_sqrt_grouping[,type := factor(ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")))]

#visualize
#just natural sites
PCoA_macroinvert_allsite_points_natural_only <- ggplot(PCoA_macroinvert_sqrt_grouping[DepthZone != "AR"]) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(aes(Dim1, Dim2, color = DepthZone, shape = DepthZone), size = 4, fill = "#FE6100") +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  lims(x = c(-0.7,0.6),y = c(-0.6,0.6)) +
  labs(x = paste0("PCoA dim 1 (",round(macroinvert_PCoA$eig[1],1),"%)"),y = paste0("Macroinvertebrate\nPCoA dim 2 (",round(macroinvert_PCoA$eig[2],1),"%)"), shape = "Depth zone",color = "Depth zone") +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  annotate(geom = "text",label ="Square-root transformation\nDistance: Bray-Curtis dissimilarity", 
           x = 0.6, y = -0.6, size = 4, hjust = 1) +
  theme_classic() +
  theme(legend.position = c(0.5,0.95), legend.direction = "horizontal",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    color = guide_legend(override.aes = list(size = 6), title.position = "top", title.hjust = 0.5, reverse = T),
    shape = guide_legend(override.aes = list(size = 6), title.position = "top", title.hjust = 0.5, reverse = T),
  )

ggsave(PCoA_macroinvert_allsite_points_natural_only, path = "figures", filename = "PCoA_macroinvert_allsite_points_natural_only.jpg", height = 8, width = 12, unit = "in", dpi = 300)


#mainland vs island
#just natural sites
PCoA_macroinvert_allsite_points_natural_only_island_mainland <- ggplot(PCoA_macroinvert_sqrt_grouping[DepthZone != "AR"]) +
  geom_point(aes(Dim1, Dim2, fill = type), shape = 24, color = "black", size = 3) +
  scale_fill_manual(values = c("black","white")) +
  lims(x = c(-0.7,0.6),y = c(-0.6,0.6)) +
  labs(x = "PCoA dim 1",y = "PCoA dim 2", fill = "Reef type") +
  scale_shape_manual(values = c(15,17,19,23), guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(legend.direction = "horizontal",
        legend.position = c(.55, .95),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  guides(
    fill = guide_legend(override.aes = list(size = 6), title.position = "top",title.hjust = 0.5)
  )

ggsave(PCoA_macroinvert_allsite_points_natural_only_island_mainland, path = "figures", filename = "PCoA_macroinvert_allsite_points_natural_only_island_mainland.jpg", height = 8, width = 12, unit = "in", dpi = 300)


#merge plot
PCoA_macroinvert_centroids <- ggplot() +
  stat_ellipse(data = PCoA_macroinvert_sqrt_grouping[DepthZone != "AR"], geom = "polygon", aes(Dim1, Dim2, fill = DepthZone), color = "white",alpha = 0.2) +
  stat_ellipse(data = PCoA_macroinvert_sqrt_grouping[DepthZone != "AR"], geom = "polygon", aes(Dim1, Dim2,linetype = type), color = "gray32",fill = NA) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  geom_point(data = PCoA_macroinvert_sqrt_grouping[DepthZone == "AR"], aes(Dim1, Dim2, shape = DepthZone_wAR),color = "black", size = 4) +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100"), guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values = c(9,10)) +
  annotate(geom = "text",label = "I", x = -0.5, y = 0.3, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "M", x = -0.3, y = 0.31, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "O", x = 0, y = 0.38, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "D", x = 0.3, y = 0.4, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "Island", x = -0.2, y = 0.45, size = 5, fontface = 'bold') +
  annotate(geom = "text",label = "Mainland", x = 0.45, y = 0.15, size = 5, fontface = 'bold') +
  guides(fill = "none",linetype = "none",shape = guide_legend(title.position = "top",title.hjust = 0.5)) +
  lims(x = c(-0.7,0.6),y = c(-0.6,0.6)) +
  labs(x = "PCoA dim 1",y = "PCoA dim 2", shape = "Depth zone",color = "Depth zone", fill = "Depth zone", linetype = "Reef type") +
  theme_classic() +
  theme(legend.position = c(0.5,0.95), legend.justification = "center",legend.background = element_blank(), legend.direction = "horizontal",
        legend.text = element_text(size = 14),   # Increase legend text size
        legend.key.size = unit(1.5, "lines"),
        legend.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))     # Increase legend key size

#plot natural reefs next to artificial reef highlight plot
PCoA_macroinvert_allsite_AR_natural_merge <- plot_grid(PCoA_macroinvert_allsite_points_natural_only + theme(legend.position = "null"),
                                                PCoA_macroinvert_allsite_points_natural_only_island_mainland + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "null"),
                                                PCoA_macroinvert_centroids+theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "null"), ncol = 3, align = "hv", labels = c("g.","h.","i."), label_size = 20)


ggsave(PCoA_macroinvert_allsite_AR_natural_merge, path = "figures", filename = "PCoA_macroinvert_allsite_AR_natural_merge.jpg", width =18, height = 6)

#full merge
PCoA_allsite_AR_natural_merge_taxagroups <- plot_grid(PCoA_fish_allsite_AR_natural_merge,
          PCoA_kelp_allsite_AR_natural_merge,
          PCoA_macroinvert_allsite_AR_natural_merge, ncol = 1, align = "hv")

ggsave(PCoA_allsite_AR_natural_merge_taxagroups, path = "figures", filename = "PCoA_allsite_AR_natural_merge_taxagroups.jpg", width =18, height = 18)


#PLOT BY REGION FOR SUPPLEMENT
PCoA_sqrt_grouping[,Region := ifelse(grepl("Santa Monica",Site) == T, "Santa Monica Bay",as.character(Region))]

#Set order by latitude and mainland/island
PCoA_sqrt_grouping[,Region := factor(Region, levels = c( "Malibu"    , 
                                                    "Santa Monica Bay",
                                                    "Palos Verdes"           ,
                                                    "Orange and North County",
                                                    "La Jolla and Point Loma",
                                                    "San Clemente Island"    ,
                                                    "Santa Barbara Island"   ,
                                                    "Santa Catalina Island"))]

#AR or NO AR column
PCoA_sqrt_grouping[,Natural_Artificial := ifelse(DepthZone == "AR","Artificial reef","Natural reef")]

#Add in MPA
MPA_site_key <- readRDS(file.path("keys","MPA_site_key.rds"))
PCoA_sqrt_grouping <- MPA_site_key[PCoA_sqrt_grouping, on = "Site"]

#Visualize by Region
PCoA_allspp_allsite_byreg <- ggplot(PCoA_sqrt_grouping) +
  geom_point(aes(Dim1, Dim2, color = Region, shape = Natural_Artificial), size = 4) +
  scale_color_manual(values = c("maroon2"   ,   "deepskyblue2", "deepskyblue4","rosybrown2" , "green2"   ,   "darkorchid2" , "darkseagreen4"    ,    "brown1" )) +
  scale_shape_manual(values = c(7,19)) +
  labs(x = "PCoA dim 1",y = "PCoA dim 2", shape = "Reef type") +
  guides(color = guide_legend(ncol = 2)) +
  theme_classic() +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.position = "top")

ggsave(PCoA_allspp_allsite_byreg, path = "figures", filename = "PCoA_allspp_allsite_byreg.jpg", height = 8, width = 12, unit = "in", dpi = 300)


#Visualize by MPA status
PCoA_allspp_allsite_byMPA <- ggplot(PCoA_sqrt_grouping) +
  geom_point(aes(Dim1, Dim2, fill = MPA_overlap, shape = MPA_overlap), size = 4, alpha = 0.8) +
  labs(x = "PCoA dim 1",y = "PCoA dim 2", fill = "MPA", shape = "MPA") +
  scale_fill_manual(values = c("black","white")) +
  scale_shape_manual(values = c(23,21)) +
  theme_classic() +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.position = "top")

ggsave(PCoA_allspp_allsite_byMPA, path = "figures", filename = "PCoA_allspp_allsite_byMPA.jpg", height = 8, width = 12, unit = "in", dpi = 300)

#Merge side by side
PCoA_allspp_allsite_MPA_REG_merge <- cowplot::plot_grid(PCoA_allspp_allsite_byMPA, PCoA_allspp_allsite_byreg, ncol = 1, labels = c("a.","b."))

#Save
ggsave(PCoA_allspp_allsite_MPA_REG_merge, path = "figures", filename = "PCoA_allspp_allsite_MPA_REG_merge.jpg", height = 17, width = 11, unit = "in", dpi = 300)

#####################
#PERMANOVA, Permutational Multivariate Analysis of Variance (perMANOVA)
#####################

#Spp only
dat_averages_bysite.wide.spp.l

#Spp and attributes
dat_averages_bysite.wide.spp.attributes

#which are artificial reef sites?
AR_ref <- which(dat_averages_bysite.wide.spp.attributes$DepthZone == "ARM")

#which are outer or deep sites (for comparison with ARs)
outer_deep_ref <- which(dat_averages_bysite.wide.spp.attributes$DepthZone %in% c("Outer","Deep"))

#excluding AR
dat_averages_bysite.wide.spp.l.NATURAL <- dat_averages_bysite.wide.spp.l[-AR_ref,]
dat_averages_bysite.wide.spp.attributes_NATURAL <- dat_averages_bysite.wide.spp.attributes[DepthZone!="ARM"]

#PERMANOVA FOR ALL SPECIES
#permanova not including ARM as depth zone
permanova_allspp_natural_depthzone <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL ~ dat_averages_bysite.wide.spp.attributes_NATURAL$DepthZone,
  method = "bray",
  permutations = 9999
  ,
  by = "margin"
)

permanova_allspp_natural_depthzone

#Depth accounts for 19% of variation in natural reef sites#PERMANOVA FOR ALL SPECIES
#permanova not including ARM as depth zone
permanova_allspp_natural_depthzone <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL ~ dat_averages_bysite.wide.spp.attributes_NATURAL$DepthZone,
  method = "bray",
  permutations = 9999
  ,
  by = "margin"
)

permanova_allspp_natural_depthzone

#Depth accounts for 19% of variation in natural reef sites#PERMANOVA FOR ALL SPECIES

#PERMANOVA comparing MPA to non-MPA sites
permanova_allspp_natural_MPAoverlap <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL ~ dat_averages_bysite.wide.spp.attributes_NATURAL$MPA_overlap,
  method = "bray",
  permutations = 9999
  ,
  by = "margin"
)

permanova_allspp_natural_MPAoverlap

#Depth accounts for 19% of variation in natural reef sites

#island mainland and depth zone
permanova_allspp_natural_isl_main <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL ~ dat_averages_bysite.wide.spp.attributes_NATURAL$type, #20% of variation
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_allspp_natural_isl_main

#island/mainland accounts for 14% of variation

#Interaction
#island mainland and depth zone
permanova_allspp_natural_isl_main_depth <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL ~ dat_averages_bysite.wide.spp.attributes_NATURAL$type * dat_averages_bysite.wide.spp.attributes_NATURAL$DepthZone, #20% of variation
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_allspp_natural_isl_main_depth

#island/mainland * depth interaction = 2% of variation

#PAIRED TESTS (just repeat PERMANOVA but for only two groups at a time)
#Deep and outer zones
rows_keep <- which(dat_averages_bysite.wide.spp.attributes_NATURAL[,DepthZone] %in% c("Deep","Outer"))
dat_averages_bysite.wide.spp.attributes_NATURAL_deep_outer <- dat_averages_bysite.wide.spp.attributes_NATURAL[rows_keep,]
dat_averages_bysite.wide.spp.l.NATURAL.trim.deep_outer <- dat_averages_bysite.wide.spp.l.NATURAL[rows_keep,]

permanova_allspp_naturalreef_deepouter <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL.trim.deep_outer ~ dat_averages_bysite.wide.spp.attributes_NATURAL_deep_outer$DepthZone,
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_allspp_naturalreef_deepouter #R2 = 16, p-value = 0.0004

#Deep and middle zones
rows_keep <- which(dat_averages_bysite.wide.spp.attributes_NATURAL[,DepthZone] %in% c("Deep","Middle"))
dat_averages_bysite.wide.spp.attributes_NATURAL_deep_middle <- dat_averages_bysite.wide.spp.attributes_NATURAL[rows_keep,]
dat_averages_bysite.wide.spp.l.NATURAL.trim.deep_middle <- dat_averages_bysite.wide.spp.l.NATURAL[rows_keep,]

permanova_allspp_naturalreef_deepmiddle <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL.trim.deep_middle ~ dat_averages_bysite.wide.spp.attributes_NATURAL_deep_middle$DepthZone,
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_allspp_naturalreef_deepmiddle #R2 = 17, p-value = 0.0004

#Deep and inner
rows_keep <- which(dat_averages_bysite.wide.spp.attributes_NATURAL[,DepthZone] %in% c("Deep","Inner"))
dat_averages_bysite.wide.spp.attributes_NATURAL_deep_inner <- dat_averages_bysite.wide.spp.attributes_NATURAL[rows_keep,]
dat_averages_bysite.wide.spp.l.NATURAL.trim.deep_inner <- dat_averages_bysite.wide.spp.l.NATURAL[rows_keep,]

permanova_allspp_naturalreef_deepinner <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL.trim.deep_inner ~ dat_averages_bysite.wide.spp.attributes_NATURAL_deep_inner$DepthZone,
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_allspp_naturalreef_deepinner #R2 = 18, p-value = 0.0004

#Outer and inner
rows_keep <- which(dat_averages_bysite.wide.spp.attributes_NATURAL[,DepthZone] %in% c("Outer","Inner"))
dat_averages_bysite.wide.spp.attributes_NATURAL_outer_inner <- dat_averages_bysite.wide.spp.attributes_NATURAL[rows_keep,]
dat_averages_bysite.wide.spp.l.NATURAL.trim.outer_inner <- dat_averages_bysite.wide.spp.l.NATURAL[rows_keep,]

permanova_allspp_naturalreef_outerinner <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL.trim.outer_inner ~ dat_averages_bysite.wide.spp.attributes_NATURAL_outer_inner$DepthZone,
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_allspp_naturalreef_outerinner #R2 = 17, p-value = 0.0004

#Outer and middle
rows_keep <- which(dat_averages_bysite.wide.spp.attributes_NATURAL[,DepthZone] %in% c("Outer","Middle"))
dat_averages_bysite.wide.spp.attributes_NATURAL_outer_middle <- dat_averages_bysite.wide.spp.attributes_NATURAL[rows_keep,]
dat_averages_bysite.wide.spp.l.NATURAL.trim.outer_middle <- dat_averages_bysite.wide.spp.l.NATURAL[rows_keep,]

permanova_allspp_naturalreef_outermiddle <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL.trim.outer_middle ~ dat_averages_bysite.wide.spp.attributes_NATURAL_outer_middle$DepthZone,
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_allspp_naturalreef_outermiddle #R2 = 17, p-value = 0.0004

#Middle and inner
rows_keep <- which(dat_averages_bysite.wide.spp.attributes_NATURAL[,DepthZone] %in% c("Middle","Inner"))
dat_averages_bysite.wide.spp.attributes_NATURAL_middle_inner <- dat_averages_bysite.wide.spp.attributes_NATURAL[rows_keep,]
dat_averages_bysite.wide.spp.l.NATURAL.trim.middle_inner <- dat_averages_bysite.wide.spp.l.NATURAL[rows_keep,]

permanova_allspp_naturalreef_middleinner <- adonis2(
  dat_averages_bysite.wide.spp.l.NATURAL.trim.middle_inner ~ dat_averages_bysite.wide.spp.attributes_NATURAL_middle_inner$DepthZone,
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_allspp_naturalreef_middleinner #R2 = 18, p-value = 0.0004

#PERMDIS FOR DISPERSION
#distance matrix
AllSpp_NATURALREEF.dist <- vegdist(dat_averages_bysite.wide.spp.l.NATURAL, method = "bray")

#Island/mainland beta dispersion
PERMDISP2_AllSpp_NATURALREEF_isl_main <- betadisper(d = AllSpp_NATURALREEF.dist, group = dat_averages_bysite.wide.spp.attributes_NATURAL$type, type = "centroid")

#MPA non-MPA beta dispersion
PERMDISP2_AllSpp_NATURALREEF_MPAoverlap <- betadisper(d = AllSpp_NATURALREEF.dist, group = dat_averages_bysite.wide.spp.attributes_NATURAL$MPA_overlap, type = "centroid")

#Depth zone beta dispersion
PERMDISP2_AllSpp_NATURALREEF_depthzone <- betadisper(d = AllSpp_NATURALREEF.dist, group = dat_averages_bysite.wide.spp.attributes_NATURAL$DepthZone, type = "centroid")

#anova of depth zones
anova(PERMDISP2_AllSpp_NATURALREEF_depthzone)

#No significant difference between depth zones

#ANOVA for MPA and non-MPA
anova(PERMDISP2_AllSpp_NATURALREEF_MPAoverlap)

#ANOVA (compare dispersion) of island versus mainland reefs
anova(PERMDISP2_AllSpp_NATURALREEF_isl_main)

  #No significant difference between island and mainland reefs


################TAXA SPECIFIC PERMANOVA ####

#Fish PERMANOVA ####

community_matrix_fish_sqrt_attributes <- cbind(dat_fish_averages_bysite.wide[,c(1:3)], dat_fish_averages_bysite.wide.spp.l)

#Add island.mainland
community_matrix_fish_sqrt_attributes[,type := ifelse(DepthZone %in% c("ARM","Module"),"ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]

#MPA designation
community_matrix_fish_sqrt_attributes <- MPA_site_key[community_matrix_fish_sqrt_attributes, on = "Site"]

#excluding ARs
fish_natural_reef_rows <- which(community_matrix_fish_sqrt_attributes[,DepthZone] != "ARM")

community_matrix_fish_sqrt_NATURAL_attributes <- community_matrix_fish_sqrt_attributes[fish_natural_reef_rows,]
dat_fish_averages_bysite.wide.NATURAL.spp.l <- dat_fish_averages_bysite.wide.spp.l[fish_natural_reef_rows,]

permanova_fish_NATURAL_depthzone <- adonis2(
  dat_fish_averages_bysite.wide.NATURAL.spp.l ~ community_matrix_fish_sqrt_NATURAL_attributes$DepthZone,
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_fish_NATURAL_depthzone

#Depth accounts for 18% of variation in fish composition at natural reef sites

#island mainland
permanova_fish_NATURAL_isl_main <- adonis2(
  dat_fish_averages_bysite.wide.NATURAL.spp.l ~ community_matrix_fish_sqrt_NATURAL_attributes$type, 
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_fish_NATURAL_isl_main


#island/mainland accounts for 11% of variation
permanova_fish_NATURAL_isl_main_depthzone <- adonis2(
  dat_fish_averages_bysite.wide.NATURAL.spp.l ~ community_matrix_fish_sqrt_NATURAL_attributes$type * community_matrix_fish_sqrt_NATURAL_attributes$DepthZone, 
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_fish_NATURAL_isl_main_depthzone

#interaction is not significant

#PERMDIS FOR DISPERSION
#distance matrix
Fish_NATURALREEF.dist <- vegdist(dat_fish_averages_bysite.wide.NATURAL.spp.l, method = "bray")

#Island/mainland beta dispersion
PERMDISP2_Fish_NATURALREEF_isl_main <- betadisper(d = Fish_NATURALREEF.dist, group = community_matrix_fish_sqrt_NATURAL_attributes$type, type = "centroid")

#Anova of island vs. mainland
anova(PERMDISP2_Fish_NATURALREEF_isl_main)

#depth betadisp
PERMDISP2_Fish_NATURALREEF_depthzone <- betadisper(d = Fish_NATURALREEF.dist, group = community_matrix_fish_sqrt_NATURAL_attributes$DepthZone, type = "centroid")

#anova of island vs. mainland
anova(PERMDISP2_Fish_NATURALREEF_depthzone)


#Macroinvert PERMANOVA ####

community_matrix_macroinvert_sqrt_attributes <- cbind(dat_macroinvert_averages_bysite.wide[,c(1:3)], dat_macroinvert_averages_bysite.wide.spp.l)

#Add island.mainland
community_matrix_macroinvert_sqrt_attributes[,type := ifelse(DepthZone %in% c("ARM","Module"),"ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]

#MPA designation
community_matrix_macroinvert_sqrt_attributes <- MPA_site_key[community_matrix_macroinvert_sqrt_attributes, on = "Site"]

#excluding ARs
macroinvert_natural_reef_rows <- which(community_matrix_macroinvert_sqrt_attributes[,DepthZone] != "ARM")

community_matrix_macroinvert_sqrt_NATURAL_attributes <- community_matrix_macroinvert_sqrt_attributes[macroinvert_natural_reef_rows,]
dat_macroinvert_averages_bysite.wide.NATURAL.spp.l <- dat_macroinvert_averages_bysite.wide.spp.l[macroinvert_natural_reef_rows,]

permanova_macroinvert_NATURAL_depthzone <- adonis2(
  dat_macroinvert_averages_bysite.wide.NATURAL.spp.l ~ community_matrix_macroinvert_sqrt_NATURAL_attributes$DepthZone,
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_macroinvert_NATURAL_depthzone

#Depth accounts for 16% of variation in macroinvert composition at natural reef sites

#Island versus mainland
permanova_macroinvert_NATURAL_isl_main <- adonis2(
  dat_macroinvert_averages_bysite.wide.NATURAL.spp.l ~ community_matrix_macroinvert_sqrt_NATURAL_attributes$type, 
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_macroinvert_NATURAL_isl_main

#island/mainland accounts for 14% of variation

permanova_macroinvert_NATURAL_isl_main_depthzone <- adonis2(
  dat_macroinvert_averages_bysite.wide.NATURAL.spp.l ~ community_matrix_macroinvert_sqrt_NATURAL_attributes$type * community_matrix_macroinvert_sqrt_NATURAL_attributes$DepthZone, 
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_macroinvert_NATURAL_isl_main_depthzone

#interaction is significant

#PERMDIS FOR DISPERSION
#distance matrix
Macroinvert_NATURALREEF.dist <- vegdist(dat_macroinvert_averages_bysite.wide.NATURAL.spp.l, method = "bray")

#Island/mainland beta dispersion
PERMDISP2_Macroinvert_NATURALREEF_isl_main <- betadisper(d = Macroinvert_NATURALREEF.dist, group = community_matrix_macroinvert_sqrt_NATURAL_attributes$type, type = "centroid")

#Anova of island vs. mainland
anova(PERMDISP2_Macroinvert_NATURALREEF_isl_main)

#depth betadisp
PERMDISP2_Macroinvert_NATURALREEF_depthzone <- betadisper(d = Macroinvert_NATURALREEF.dist, group = community_matrix_macroinvert_sqrt_NATURAL_attributes$DepthZone, type = "centroid")

#anova of island vs. mainland
anova(PERMDISP2_Macroinvert_NATURALREEF_depthzone)

############ARM vs NAT

#Reduce to mainland outer and deep reefs only
dat_averages_bysite.wide.spp.l.OUTERDEEPARM <- dat_averages_bysite.wide.spp.l[c(outer_deep_ref,AR_ref),]
dat_averages_bysite.wide.spp.attributes_OUTERDEEPARM <- dat_averages_bysite.wide.spp.attributes[c(outer_deep_ref,AR_ref),]

#Designate arm vs natural
dat_averages_bysite.wide.spp.attributes_OUTERDEEPARM[,ARM_NATURAL := ifelse(DepthZone == "ARM","ARM","NATURAL")]

#Perform PERMANOVA
permanova_allspp_ARvsNAT <- adonis2(
  #Reduce to mainland outer and deep reefs only
  dat_averages_bysite.wide.spp.l.OUTERDEEPARM ~ dat_averages_bysite.wide.spp.attributes_OUTERDEEPARM$ARM_NATURAL, 
  method = "bray",
  permutations = 9999,
  by = "margin"
)
permanova_allspp_ARvsNAT

#distance matrix
  dat_averages_bysite.wide.spp.l.OUTERDEEPARM.dist <- vegdist(  dat_averages_bysite.wide.spp.l.OUTERDEEPARM, method = "bray")

#natural vs. AR betadisp
permdis_allspp_ar_vs_nat <- betadisper(d =   dat_averages_bysite.wide.spp.l.OUTERDEEPARM.dist, group =  dat_averages_bysite.wide.spp.attributes_OUTERDEEPARM$ARM_NATURAL, type = "centroid")

#anova of arm vs. not arm
anova(permdis_allspp_ar_vs_nat)

