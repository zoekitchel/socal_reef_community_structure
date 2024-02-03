# CREATION DATE 28 Jan 2023

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Create depth zone summaries

#############################
##Setup
#############################
library(ggplot2)
library(data.table)
library(RColorBrewer)

source(file.path("functions","return_spptaxonomy_function.R"))

########################
##Load data
########################
dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))
dat_fish_site_averages <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages.rds"))
dat_macroinvert_site_averages <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages.rds"))
dat_kelp_site_averages <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages.rds"))

########################
##Averaged across all sites, top species per depth zone
########################

dat_fish_averages <- dat_fish_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2),
                                               mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                                       .(Species, DepthZone)] 
dat_macroinvert_averages <- dat_macroinvert_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                        .(BenthicReefSpecies, DepthZone)]  
dat_kelp_averages <- dat_kelp_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                        .(BenthicReefSpecies, DepthZone)]  

#Add column for whether species is in top 5 species (by density) for that zone
setkey(dat_fish_averages, mean_depthzone_density_m2)
dat_fish_averages_top5density <- dat_fish_averages[, tail(.SD, 5), .(DepthZone)]
dat_fish_averages_top5density <- dat_fish_averages_top5density[,1:3]

setkey(dat_macroinvert_averages, mean_depthzone_density_m2)
dat_macroinvert_averages_top5density <- dat_macroinvert_averages[, tail(.SD, 5), .(DepthZone)]
dat_macroinvert_averages_top5density <- dat_macroinvert_averages_top5density[,1:3]

setkey(dat_kelp_averages, mean_depthzone_density_m2)
dat_kelp_averages_top5density <- dat_kelp_averages[, tail(.SD, 5), .(DepthZone)]
dat_kelp_averages_top5density <- dat_kelp_averages_top5density[,1:3]

#Add column for whether species is in top 5 species (by biomass) for that zone
setkey(dat_fish_averages, mean_depthzone_wt_density_g_m2)
dat_fish_averages_top5biomass <- dat_fish_averages[, tail(.SD, 5), .(DepthZone)]
dat_fish_averages_top5biomass <- dat_fish_averages_top5biomass[,c(1:2,4)]

############################################################
#Stacked barplots for all sites together, only top 5 species
############################################################

#fish density
fish_density_top5_stacked <- ggplot(dat_fish_averages_top5density, aes(fill=Species, y=mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Average density (per 100m^2)") +
  scale_fill_manual(values = brewer.pal(n = length(unique(dat_fish_averages_top5density$Species)),"Set3")) +
  theme_classic()

ggsave(fish_density_top5_stacked, path = file.path("figures"), filename = "fish_density_top5_stacked.jpg", height = 3, width = 5, units = "in")

#fish biomass
fish_biomass_top5_stacked <- ggplot(dat_fish_averages_top5biomass, aes(fill=Species, y=mean_depthzone_wt_density_g_m2*100/1000, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Average biomass (kg per 100m^2)") +
  scale_fill_manual(values = brewer.pal(n = length(unique(dat_fish_averages_top5biomass$Species)),"Set3")) +
  theme_classic()

ggsave(fish_biomass_top5_stacked, path = file.path("figures"), filename = "fish_biomass_top5_stacked.jpg", height = 3, width = 5, units = "in")

#macroinvert density
macroinvert_density_top5_stacked <- ggplot(dat_macroinvert_averages_top5density, aes(fill=BenthicReefSpecies, y=mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Average density (per 100m^2)", fill = "Macroinvert species") +
  scale_fill_manual(values = brewer.pal(n = length(unique(dat_macroinvert_averages_top5density$BenthicReefSpecies)),"Set3")) +
  theme_classic()

ggsave(macroinvert_density_top5_stacked, path = file.path("figures"), filename = "macroinvert_density_top5_stacked.jpg", height = 3, width = 5, units = "in")

#kelp density
kelp_density_top5_stacked <- ggplot(dat_kelp_averages_top5density, aes(fill=BenthicReefSpecies, y=mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Average density (per 100m^2)", fill = "Kelp species") +
  scale_fill_manual(values = brewer.pal(n = length(unique(dat_kelp_averages_top5density$BenthicReefSpecies)),"Set3")) +
  theme_classic()

ggsave(kelp_density_top5_stacked, path = file.path("figures"), filename = "kelp_density_top5_stacked.jpg", height = 3, width = 5, units = "in")


#fish density
fish_density_top5_stacked <- ggplot(dat_fish_averages_top5density, aes(fill=Species, y=mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Average density (per 100m^2)") +
  scale_fill_manual(values = brewer.pal(n = length(unique(dat_fish_averages_top5density$Species)),"Set3")) +
  theme_classic()

ggsave(fish_density_top5_stacked, path = file.path("figures"), filename = "fish_density_top5_stacked.jpg", height = 3, width = 5, units = "in")

############################################################
#Stacked barplots for all sites together, all species grouped by taxonomy
############################################################

#palette with 15 colors
pal_15 <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
            "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
            "#CCEBC5", "#FFED6F","darksalmon","lightskyblue","turquoise")

#link spp names to spp taxonomy
fish_list <- get_taxa(taxon_list = unique(dat_fish_averages$Species)) #CHECK WEIRD ORDER
fish_list[,Species := query]

dat_fish_averages <- dat_fish_averages[fish_list, on = "Species"]

#Sum density across orders
dat_fish_summed_order <- dat_fish_averages[,.(mean_depthzone_density_m2_summed_by_order = sum(mean_depthzone_density_m2),
                                               mean_depthzone_wt_density_g_m2_summed_by_order=sum(mean_depthzone_wt_density_g_m2)),
                                            .(order,DepthZone)] 

#fish density
fish_density_order_stacked <- ggplot(dat_fish_summed_order,
                                     aes(fill=order, y=mean_depthzone_density_m2_summed_by_order*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Summed average density (per 100m^2)", fill = "Order") +
  scale_fill_manual(values = pal_15) +
  theme_classic()

ggsave(fish_density_order_stacked, path = file.path("figures"), filename = "fish_density_order_stacked.jpg",
       height = 5, width = 7, units = "in")

#fish biomass
fish_biomass_order_stacked <- ggplot(dat_fish_summed_order, aes(fill=order, y=mean_depthzone_wt_density_g_m2_summed_by_order*100/1000, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Summed average biomass (kg per 100m^2)") +
  scale_fill_manual(values = pal_15) +
  theme_classic()

ggsave(fish_biomass_order_stacked, path = file.path("figures"), filename = "fish_biomass_order_stacked.jpg",
       height = 5, width = 7, units = "in")


############################################################
#Now, average biomass and density by Island vs ARM vs Coast
############################################################

#New Column Identifying ARM vs Island vs Natural Coast
dat_fish_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island"),"Island","Natural mainland"))]
dat_macroinvert_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island"),"Island","Natural mainland"))]
dat_kelp_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island"),"Island","Natural mainland"))]

#average by type and depth zone
dat_fish_averages_sitetype <- dat_fish_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2),
                                               mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                            .(Species, DepthZone, type)] 
dat_macroinvert_averages_sitetype <- dat_macroinvert_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                          .(BenthicReefSpecies, DepthZone, type)]  
dat_kelp_averages_sitetype <- dat_kelp_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                            .(BenthicReefSpecies, DepthZone, type)]  

#Add column for whether species is in top 5 species (by density) for that zone
setkey(dat_fish_averages_sitetype, mean_depthzone_density_m2)
dat_fish_averages_sitetype_top5density <- dat_fish_averages_sitetype[, tail(.SD, 5), .(DepthZone, type)]
dat_fish_averages_sitetype_top5density <- dat_fish_averages_sitetype_top5density[,1:4]

setkey(dat_macroinvert_averages_sitetype, mean_depthzone_density_m2)
dat_macroinvert_averages_sitetype_top5density <- dat_macroinvert_averages_sitetype[, tail(.SD, 5), .(DepthZone, type)]

setkey(dat_kelp_averages_sitetype, mean_depthzone_density_m2)
dat_kelp_averages_sitetype_top5density <- dat_kelp_averages_sitetype[, tail(.SD, 5), .(DepthZone, type)]

#Add column for whether species is in top 5 species (by biomass) for that zone and sitetype
setkey(dat_fish_averages_sitetype, mean_depthzone_wt_density_g_m2)
dat_fish_averages_sitetype_top5biomass <- dat_fish_averages_sitetype[, tail(.SD, 5), .(DepthZone, type)]
dat_fish_averages_sitetype_top5biomass <- dat_fish_averages_sitetype_top5biomass[,c(1:3,5)]

############################################################
#Stacked barplots for average biomass and density by site type and depth, only top 5 species
############################################################

#palette with 13 colors
pal_13 <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
            "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
            "#CCEBC5", "#FFED6F","darksalmon")

#fish density
fish_density_top5_sitetype_stacked <- ggplot(dat_fish_averages_sitetype_top5density, aes(fill=Species, y=mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Average density (per 100m^2)") +
  scale_fill_manual(values = brewer.pal(n = length(unique(dat_fish_averages_sitetype_top5density$Species)),"Set3")) +
  facet_grid(~type, scales = "free_x", space = "free") +
  theme_classic()

ggsave(fish_density_top5_sitetype_stacked, path = file.path("figures"), filename = "fish_density_top5_sitetype_stacked.jpg", height = 5, width = 8, units = "in")

#fish biomass
fish_biomass_top5_sitetype_stacked <- ggplot(dat_fish_averages_sitetype_top5biomass, aes(fill=Species, y=mean_depthzone_wt_density_g_m2*100/1000, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Average biomass (kg per 100m^2)") +
  scale_fill_manual(values = brewer.pal(n = length(unique(dat_fish_averages_sitetype_top5biomass$Species)),"Set3")) +
  facet_grid(~type, scales = "free_x", space = "free") +
  theme_classic()

ggsave(fish_biomass_top5_sitetype_stacked, path = file.path("figures"), filename = "fish_biomass_top5_sitetype_stacked.jpg", height = 5, width = 8, units = "in")

#macroinvert density
macroinvert_density_top5_sitetype_stacked <- ggplot(dat_macroinvert_averages_sitetype_top5density, aes(fill=BenthicReefSpecies, y=mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Average density (per 100m^2)", fill = "Macroinvert species") +
  scale_fill_manual(values = pal_13) +
  facet_grid(~type, scales = "free_x", space = "free") +
  theme_classic()

ggsave(macroinvert_density_top5_sitetype_stacked, path = file.path("figures"), filename = "macroinvert_top5_sitetype_stacked.jpg", height = 5, width = 8, units = "in")

#kelp density
kelp_density_top5_sitetype_stacked <- ggplot(dat_kelp_averages_sitetype_top5density, aes(fill=BenthicReefSpecies, y=mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Average density (per 100m^2)", fill = "Kelp species") +
  scale_fill_manual(values = brewer.pal(n = length(unique(dat_kelp_averages_sitetype_top5density$BenthicReefSpecies)),"Set3")) +
  facet_grid(~type, scales = "free_x", space = "free") +
  theme_classic()

ggsave(kelp_density_top5_sitetype_stacked, path = file.path("figures"), filename = "kelp_density_top5_sitetype_stacked.jpg", height = 5, width = 8, units = "in")

##################################################################################
#Visual summaries for OSM poster
##################################################################################
##################################################################################
#TO DO: average biomass and density for natural reefs only, MPA vs non-MPA (not sure how to determine this right now)
##################################################################################
