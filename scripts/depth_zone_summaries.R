# CREATION DATE 28 Jan 2023

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Create depth zone summaries

#############################
##Setup
#############################
library(ggplot2)
library(data.table)
library(cowplot)
library(RColorBrewer)

source(file.path("functions","return_spptaxonomy_function.R"))

########################
##Load data
########################
dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))
dat_fish_site_averages <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages.rds"))
dat_macroinvert_site_averages <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages.rds"))
dat_kelp_site_averages <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages.rds"))

#Pull in lat lon from priority sites (I use 2022 list here)
VRG_priority_lat_lon <- fread(file.path("VRG_sites","2022_DiveSitePriority_Coor_List.csv"))

#take average lat and lon across matching sites
VRG_priority_lat_lon[,Latitude:= mean(Lat),.(Site)][,Longitude := mean(Lon),.(Site)]

VRG_priority_lat_lon.r <- unique(VRG_priority_lat_lon[,.(Site,Latitude, Longitude, ARM)])

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
#fish density
dat_fish_averages[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,Species,"Other"), .(DepthZone)]
dat_fish_averages[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
dat_fishdensity_averages_deep_ar.u <- unique(dat_fish_averages[,.(Species_top5, DepthZone, summed_mean_depthzone_density_m2)])

#manually set factor order for plotting
dat_fishdensity_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Brachyistius frenatus", "Chromis punctipinnis", "Girella nigricans", "Halichoeres semicinctus", "Hypsypops rubicundus",
                                                                                   "Lythrypnus dalli", "Oxyjulis californica", "Paralabrax clathratus", "Paralabrax nebulifer", "Semicossyphus pulcher",
                                                                                 "Other"),
                                                                        labels = c("Brachyistius frenatus\n(kelp perch)", "Chromis punctipinnis\n(blacksmith damselfish)", "Girella nigricans\n(opaleye)", "Halichoeres semicinctus\n(rock wrasse)", "Hypsypops rubicundus\n(garibaldi)",
                                                                                   "Lythrypnus dalli\n(blue-banded goby)", "Oxyjulis californica\n(señorita wrasse)", "Paralabrax clathratus\n(kelp bass)", "Paralabrax nebulifer\n(barred sand bass)", "Semicossyphus pulcher\n(sheephead)",
                                                                                 "Other"))]

#macroinvert density
dat_macroinvert_averages[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,BenthicReefSpecies,"Other"), .(DepthZone)]
dat_macroinvert_averages[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
dat_macroinvert_averages_deep_ar.u <- unique(dat_macroinvert_averages[,.(Species_top5, DepthZone, summed_mean_depthzone_density_m2)])

dat_macroinvert_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Anthopleura sola", "Apostichopus parvimensis", "Centrostephanus coronatus", "Leptogorgia chilensis",
                                                                                "Megastraea undosa", "Mesocentrotus franciscanus", "Muricea californica", "Pachycerianthus fimbriatus",
                                                                                "Patiria miniata","Strongylocentrotus purpuratus", "Other"),
                                                       labels = c("Anthopleura sola\n(starburst anenome)", "Apostichopus parvimensis\n(warty sea cucumber)", "Centrostephanus coronatus\n(crowned urchin)", "Leptogorgia chilensis\n(red gorgonian)",
                                                                  "Megastraea undosa\n(wavy turban snail)", "Mesocentrotus franciscanus\n(red urchin)", "Muricea californica\n(golden gregorian)", "Pachycerianthus fimbriatus\n(tube dwelling anenome)",
                                                                  "Patiria miniata\n(bat star)","Strongylocentrotus purpuratus\n(purple urchin)", "Other"))]


#Kelp density
dat_kelp_averages[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,BenthicReefSpecies,"Other"), .(DepthZone)]
dat_kelp_averages[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
dat_kelp_averages_deep_ar.u <- unique(dat_kelp_averages[,.(Species_top5, DepthZone, summed_mean_depthzone_density_m2)])

dat_kelp_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Agarum fimbriatum", "Eisenia arborea", "Laminaria farlowii", "Macrocystis pyrifera",
                                                                         "Pterygophora californica", "Sargassum horneri","Sargassum palmeri", "Stephanocystis spp.",  "Other"))]


#Add column for whether species is in top 5 species (by biomass) for that zone, and then sum biomass for all others
#fish biomass
dat_fish_averages[, Species_top5 := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,Species,"Other"), .(DepthZone)]
dat_fish_averages[, summed_mean_depthzone_wt_density_g_m2 := sum(mean_depthzone_wt_density_g_m2), .(DepthZone, Species_top5)]
dat_fishbiomass_averages_deep_ar.u <- unique(dat_fish_averages[,.(Species_top5, DepthZone, summed_mean_depthzone_wt_density_g_m2)])

dat_fishbiomass_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Anisotremus davidsonii","Chromis punctipinnis", "Girella nigricans", "Hypsypops rubicundus",
                                                                               "Paralabrax clathratus", "Paralabrax nebulifer","Sebastes serranoides", "Semicossyphus pulcher","Stereolepis gigas",
                                                                                "Other"),
                                                       labels = c("Anisotremus davidsonii\n(sargo)", "Chromis punctipinnis\n(blacksmith damselfish)", "Girella nigricans\n(opaleye)", "Hypsypops rubicundus\n(garibaldi)",
                                                                 "Paralabrax clathratus\n(kelp bass)", "Paralabrax nebulifer\n(barred sand bass)","Sebastes serranoides\n(olive rockfish)", "Semicossyphus pulcher\n(sheephead)","Stereolepis gigas\n(giant sea bass)",
                                                                  "Other"))]



############################################################
#Stacked barplots for all sites together, if not in top 5 species, summed into OTHER category
############################################################

#fish density palette with 10 colors
pal_fishdens10 <- c("#CCEBC5",  #"Brachyistius frenatus"
                    "#80B1D3",  #"Chromis punctipinnis"
                    "#8ea489",  #"Girella nigricans"
                    "#FCCDE5",  #"Halichoeres semicinctus"
                    "#FDB462",  #"Hypsypops rubicundus"
                    "lightskyblue",  #"Lythrypnus dalli"
                   "#FFFFB3",  #"Oxyjulis californica"
                   "#8DD3C7",  #"Paralabrax clathratus"
                   "#FFED6F",  #"Paralabrax nebulifer"
                   "darksalmon",  #"Semicossyphus pulcher"
                   "#D9D9D9")  #"Other"

#fish density palette with 10 colors
pal_fishbio10 <- c("#CEAD64",#"Anisotremus davidsonii" 
                    "#80B1D3",#"Chromis punctipinnis"
                    "#8ea489",#"Girella nigricans" 
                    "#FDB462",#"Hypsypops rubicundus"
                    "#8DD3C7",#"Paralabrax clathratus"
                    "#FFED6F", #"Paralabrax nebulifer"
                   "#9D9E39", #"Sebastes serranoides"
                    "darksalmon",#"Semicossyphus pulcher"
                   "#A38389",#"Stereolepsis gigas"
                    "#D9D9D9")

#macro density palette with 11 colors
pal_macro10 <- c("#97D4BA", "#B67436","#3F4965", "#DA7E80", "#9C8074",
                   "#BA4C61", "#BF9D5D", "#D7CDAA", "#E08454", "#BC80BD","#D9D9D9")

#kelp density palette with 9 colors
pal_kelp9 <- c("#8DD3C7", "#E9D677", "#BEBADA", "#FB8072", "#3E8AAF",
               "turquoise", "#B3DE69", "#CDA3AB", "#D9D9D9")

#fish density
fish_density_top5_stacked <- ggplot(dat_fishdensity_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Average density per 100m" ^2)), fill = "Fish species") +
  scale_fill_manual(values = pal_fishdens10) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))

ggsave(fish_density_top5_stacked, path = file.path("figures"), filename = "fish_density_top5_stacked.jpg", height = 3.5, width = 5, units = "in")

#fish biomass
fish_biomass_top5_stacked <- ggplot(dat_fishbiomass_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_wt_density_g_m2*100/1000, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Average biomass in kg per 100m" ^2)), fill = "Fish species") +
  scale_fill_manual(values = pal_fishbio10) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))

ggsave(fish_biomass_top5_stacked, path = file.path("figures"), filename = "fish_biomass_top5_stacked.jpg", height = 3, width = 5, units = "in")

#macroinvert density
macroinvert_density_top5_stacked <- ggplot(dat_macroinvert_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Average density per 100m" ^2)), fill = "Macroinvertebrate species") +
  scale_fill_manual(values = pal_macro10) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))

ggsave(macroinvert_density_top5_stacked, path = file.path("figures"), filename = "macroinvert_density_top5_stacked.jpg", height = 3.5, width = 5, units = "in")

#kelp density
kelp_density_top5_stacked <- ggplot(dat_kelp_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Average density per 100m" ^2)), fill = "Kelp species") +
  scale_fill_manual(values = pal_kelp9) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))

ggsave(kelp_density_top5_stacked, path = file.path("figures"), filename = "kelp_density_top5_stacked.jpg", height = 3.5, width = 5, units = "in")


############################################################
#Stacked barplots for all sites together, all species grouped by taxonomy
############################################################

#palette with 16 colors
pal_16 <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
            "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
            "#CCEBC5", "#FFED6F","darksalmon","lightskyblue","turquoise",
            "#C5C8FD")

#link spp names to spp taxonomy
fish_list <- get_taxa(taxon_list = unique(dat_fish_averages$Species)) #CHECK WEIRD ORDER
fish_list[,Species := query]

#some strange order identification, will specify for following species
fish_list[taxa == "Anisotremus davidsonii", order := "Perciformes"]
fish_list[taxa == "Caulolatilus princeps", order := "Perciformes"]
fish_list[taxa == "Cheilotrema saturnum", order := "Acanthuriformes"]
fish_list[taxa == "Halichoeres semicinctus", order := "Labriformes"]
fish_list[taxa == "Oxyjulis californica", order := "Labriformes"]
fish_list[taxa == "Pristigenys serrula", order := "Perciformes"]
fish_list[taxa == "Bodianus pulcher", order := "Labriformes"]
fish_list[taxa == "Brachyistius frenatus", order := "Perciformes"]
fish_list[taxa == "Chromis punctipinnis", order := "Perciformes"]
fish_list[taxa == "Cymatogaster aggregata", order := "Perciformes"]
fish_list[taxa == "Embiotoca jacksoni", order := "Perciformes"]
fish_list[taxa == "Embiotoca lateralis", order := "Perciformes"]
fish_list[taxa == "Hyperprosopon argenteum", order := "Perciformes"]
fish_list[taxa == "Hypsurus caryi", order := "Perciformes"]
fish_list[taxa == "Hypsypops rubicundus", order := "Perciformes"]
fish_list[taxa == "Micrometrus minimus", order := "Perciformes"]
fish_list[taxa == "Phanerodon atripes", order := "Perciformes"]
fish_list[taxa == "Phanerodon furcatus", order := "Perciformes"]
fish_list[taxa == "Rhacochilus toxotes", order := "Perciformes"]
fish_list[taxa == "Zalembius rosaceus", order := "Perciformes"]
fish_list[taxa == "Phanerodon vacca", order := "Perciformes"]

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
  scale_fill_manual(values = pal_16) +
  theme_classic()

ggsave(fish_density_order_stacked, path = file.path("figures"), filename = "fish_density_order_stacked.jpg",
       height = 5, width = 7, units = "in")

#fish biomass
fish_biomass_order_stacked <- ggplot(dat_fish_summed_order, aes(fill=order, y=mean_depthzone_wt_density_g_m2_summed_by_order*100/1000, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = "Summed average biomass (kg per 100m^2)", fill = "Order") +
  scale_fill_manual(values = pal_16) +
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
          
          ########################
          ##Load data
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
          ##Averaged across all sites, top species per depth zone
          ########################
          
          dat_fish_averages_deep_ar <- dat_fish_site_averages_OSM.deep_ar[,.(mean_depthzone_density_m2 = mean(mean_density_m2),
                                                         mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                                      .(Species, DepthZone)] 
          dat_macroinvert_averages_deep_ar <- dat_macroinvert_site_averages_OSM.deep_ar[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                                    .(BenthicReefSpecies, DepthZone)]  
          dat_kelp_averages_deep_ar <- dat_kelp_site_averages_OSM.deep_ar[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                      .(BenthicReefSpecies, DepthZone)]  
          
            ###################
            #Same, but for SoS only
            ###################
          dat_fish_averages_SoS <- dat_fish_site_averages_OSM.deep_ar[Site == "Star of Scotland",.(mean_depthzone_density_m2 = mean(mean_density_m2),
                                                                             mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                                                          .(Species, DepthZone)][,DepthZone := "Star of Scotland"] 
          dat_macroinvert_averages_SoS <- dat_macroinvert_site_averages_OSM.deep_ar[Site == "Star of Scotland",.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                                                        .(BenthicReefSpecies, DepthZone)][,DepthZone := "Star of Scotland"]  
          dat_kelp_averages_SoS <- dat_kelp_site_averages_OSM.deep_ar[Site == "Star of Scotland",.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                                          .(BenthicReefSpecies, DepthZone)][,DepthZone := "Star of Scotland"]  
          
          ######################
          #Row bind to add Star of Scotland in
          ######################
          dat_fish_averages_deep_ar <- rbind(dat_fish_averages_deep_ar, dat_fish_averages_SoS)
          dat_macroinvert_averages_deep_ar <- rbind(dat_macroinvert_averages_deep_ar, dat_macroinvert_averages_SoS)
          dat_kelp_averages_deep_ar <- rbind(dat_kelp_averages_deep_ar, dat_kelp_averages_SoS)
          
          #Add column for whether species is in top 5 species (by density) for that zone
          #fish density
          dat_fish_averages_deep_ar[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,Species,"Other"), .(DepthZone)]
          dat_fish_averages_deep_ar[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
          dat_fishdensity_averages_deep_ar.u <- unique(dat_fish_averages_deep_ar[,.(Species_top5, DepthZone, summed_mean_depthzone_density_m2)])
          
          #manually set factor order for plotting
          dat_fishdensity_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Brachyistius frenatus", "Chromis punctipinnis", "Girella nigricans", "Halichoeres semicinctus", "Hypsypops rubicundus",
                                                                                          "Lythrypnus dalli", "Oxyjulis californica", "Paralabrax clathratus", "Paralabrax nebulifer", "Semicossyphus pulcher",
                                                                                          "Other"),
                                                                 labels = c("Brachyistius frenatus\n(kelp perch)", "Chromis punctipinnis\n(blacksmith damselfish)", "Girella nigricans\n(opaleye)", "Halichoeres semicinctus\n(rock wrasse)", "Hypsypops rubicundus\n(garibaldi)",
                                                                            "Lythrypnus dalli\n(blue-banded goby)", "Oxyjulis californica\n(señorita wrasse)", "Paralabrax clathratus\n(kelp bass)", "Paralabrax nebulifer\n(barred sand bass)", "Semicossyphus pulcher\n(sheephead)",
                                                                            "Other"))]
          
          #site type column
          dat_fishdensity_averages_deep_ar.u[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]
          
          #reorder these factors
          dat_fishdensity_averages_deep_ar.u[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
                                                                    labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]
          
          #macroinvert density
          dat_macroinvert_averages_deep_ar[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,BenthicReefSpecies,"Other"), .(DepthZone)]
          dat_macroinvert_averages_deep_ar[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
          dat_macroinvert_averages_deep_ar.u <- unique(dat_macroinvert_averages_deep_ar[,.(Species_top5, DepthZone, summed_mean_depthzone_density_m2)])
          
          dat_macroinvert_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Anthopleura sola", "Apostichopus parvimensis", "Centrostephanus coronatus","Crassadoma gigantea","Kelletia kelletii", "Leptogorgia chilensis",
                                                                                          "Megastraea undosa", "Mesocentrotus franciscanus", "Muricea californica","Muricea fruticosa","Neobernaya spadicea", "Pachycerianthus fimbriatus",
                                                                                          "Patiria miniata","Strongylocentrotus purpuratus", "Other"),
                                                                 labels = c("Anthopleura sola\n(starburst anenome)", "Apostichopus parvimensis\n(warty sea cucumber)", "Centrostephanus coronatus\n(crowned urchin)","Crassadoma gigantea\n(giant rock scallop)","Kelletia kelletii\n(Kellet's whelk)", "Leptogorgia chilensis\n(red gorgonian)",
                                                                            "Megastraea undosa\n(wavy turban snail)", "Mesocentrotus franciscanus\n(red urchin)", "Muricea californica\n(golden gregorian)","Muricea fruticosa\n(brown gregorian)","Neobernaya spadicea\n(chesnut cowrie)", "Pachycerianthus fimbriatus\n(tube dwelling anenome)",
                                                                            "Patiria miniata\n(bat star)","Strongylocentrotus purpuratus\n(purple urchin)", "Other"))]
          #site type column
          dat_macroinvert_averages_deep_ar.u[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]
          
          #reorder these factors
          dat_macroinvert_averages_deep_ar.u[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
                                                                    labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]
          
          
          #Kelp density
          dat_kelp_averages_deep_ar[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,BenthicReefSpecies,"Other"), .(DepthZone)]
          dat_kelp_averages_deep_ar[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
          dat_kelp_averages_deep_ar.u <- unique(dat_kelp_averages_deep_ar[,.(Species_top5, DepthZone, summed_mean_depthzone_density_m2)])
          
          dat_kelp_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Agarum fimbriatum", "Eisenia arborea","Egregia menziesii" , "Laminaria farlowii", "Macrocystis pyrifera","Pelagophycus porra",
                                                                                   "Pterygophora californica", "Sargassum horneri","Sargassum palmeri", "Stephanocystis spp.",  "Other"))]
          
          #site type column
          dat_kelp_averages_deep_ar.u[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]
          
          #reorder these factors
          dat_kelp_averages_deep_ar.u[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
                                                                    labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]
          
          #Add column for whether species is in top 5 species (by biomass) for that zone, and then sum biomass for all others
          #fish biomass
          dat_fish_averages_deep_ar[, Species_top5 := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,Species,"Other"), .(DepthZone)]
          dat_fish_averages_deep_ar[, summed_mean_depthzone_wt_density_g_m2 := sum(mean_depthzone_wt_density_g_m2), .(DepthZone, Species_top5)]
          dat_fishbiomass_averages_deep_ar.u <- unique(dat_fish_averages_deep_ar[,.(Species_top5, DepthZone, summed_mean_depthzone_wt_density_g_m2)])
          
          dat_fishbiomass_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Anisotremus davidsonii","Chromis punctipinnis", "Girella nigricans", "Hypsypops rubicundus",
                                                                                          "Paralabrax clathratus", "Paralabrax nebulifer","Sebastes serranoides", "Semicossyphus pulcher","Stereolepis gigas",
                                                                                          "Other"),
                                                                 labels = c("Anisotremus davidsonii\n(sargo)", "Chromis punctipinnis\n(blacksmith damselfish)", "Girella nigricans\n(opaleye)", "Hypsypops rubicundus\n(garibaldi)",
                                                                            "Paralabrax clathratus\n(kelp bass)", "Paralabrax nebulifer\n(barred sand bass)","Sebastes serranoides\n(olive rockfish)", "Semicossyphus pulcher\n(sheephead)","Stereolepis gigas\n(giant sea bass)",
                                                                            "Other"))]
          
          #site type column
          dat_fishbiomass_averages_deep_ar.u[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]
          
          #reorder these factors
          dat_fishbiomass_averages_deep_ar.u[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
                                                                    labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]
          
          
          ############################################################
          #Stacked barplots for all sites together, if not in top 5 species, summed into OTHER category
          ############################################################
          
          #fish density palette with 9 colors
          pal_fishdens9 <- c("#CCEBC5",  #"Brachyistius frenatus"
                              "#80B1D3",  #"Chromis punctipinnis"
                              "#8ea489",  #"Girella nigricans"
                           #   "#FCCDE5",  #"Halichoeres semicinctus"
                           #   "#FDB462",  #"Hypsypops rubicundus"
                              "lightskyblue",  #"Lythrypnus dalli"
                              "#FFFFB3",  #"Oxyjulis californica"
                              "#8DD3C7",  #"Paralabrax clathratus"
                              "#FFED6F",  #"Paralabrax nebulifer"
                              "darksalmon",  #"Semicossyphus pulcher"
                              "#D9D9D9")  #"Other"
          
          #fish biomass palette with 7 colors
          pal_fishbio7 <- c(
                             #"#CEAD64",#"Anisotremus davidsonii" 
                             "#80B1D3",#"Chromis punctipinnis"
                             "#8ea489",#"Girella nigricans" 
                            # "#FDB462",#"Hypsypops rubicundus"
                             "#8DD3C7",#"Paralabrax clathratus"
                             "#FFED6F", #"Paralabrax nebulifer"
                             #"#9D9E39", #"Sebastes serranoides"
                             "darksalmon",#"Semicossyphus pulcher"
                             "#A38389",#"Stereolepsis gigas"
                             "#D9D9D9")
          
          #macro density palette with 11 colors
          pal_macro11 <- c(
                            #"#97D4BA",#"Anthopleura sola"
                           #"#B67436",#"Apostichopus parvimensis"
                           "#3F4965",#"Centrostephanus coronatus"
                            "#EDC75C", #Crassadoma gigantea ADD 
                           "#814A23", # "Kelletia kelletii"
                           "#DA7E80",#"Leptogorgia chilensis"
                           "#9C8074", #"Megastraea undosa"
                           "#BA4C61", #"Mesocentrotus franciscanus"
                           "#BF9D5D",#"Muricea californica"
                           "#C9664B",#"Muricea fruticosa" Brown grogorian
                           "#411009", #"Neobernaya spadicea" Chesnut cowrie
                          # "#D7CDAA",# "Pachycerianthus fimbriatus"
                          # "#E08454",# "Patiria miniata"
                           "#BC80BD",#"Strongylocentrotus purpuratus"
                           "#D9D9D9")#"Other"

          
          #kelp density palette with 9 colors
          pal_kelp8 <- c("#8DD3C7",#"Agarum fimbriatum"
                        # "#E9D677",#"Eisenia arborea"
                         "#FCEFDB", # "Egregia menziesii" 
                         "#BEBADA",#"Laminaria farlowii"
                         "#FB8072",#"Macrocystis pyrifera"
                        "#38597A", # "Pelagophycus porra"  
                         "#3E8AAF",#"Pterygophora californica"
                        # "turquoise",# "Sargassum horneri"
                        # "#B3DE69",#"Sargassum palmeri"
                         "#CDA3AB",#"Stephanocystis spp."
                         "#D9D9D9") #other
  
          
          #fish density
          fish_density_top5deep_ar_stacked <- ggplot(dat_fishdensity_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_density_m2*100, x=`Site type`)) + 
            geom_bar(position="stack", stat="identity") +
            labs(x = "Site type", y = expression(paste("Average density per 100m" ^2)), fill = "Fish species") +
            scale_fill_manual(values = pal_fishdens9) +
            scale_y_continuous(expand = c(0,0)) +
            theme_classic() +
            theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))
          
          ggsave(fish_density_top5deep_ar_stacked, path = file.path("figures"), filename = "fish_density_top5deep_ar_stacked.jpg", height = 3.5, width = 5, units = "in")
          
          #fish biomass
          fish_biomass_top5deep_ar_stacked <- ggplot(dat_fishbiomass_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_wt_density_g_m2*100/1000, x=`Site type`)) + 
            geom_bar(position="stack", stat="identity") +
            labs(x = "Site type", y = expression(paste("Average biomass in kg per 100m" ^2)), fill = "Fish species") +
            scale_fill_manual(values = pal_fishbio7) +
            scale_y_continuous(expand = c(0,0)) +
            theme_classic() +
            theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))
          
          ggsave(fish_biomass_top5deep_ar_stacked, path = file.path("figures"), filename = "fish_biomass_top5deep_ar_stacked.jpg", height = 3, width = 5, units = "in")
          
          #macroinvert density
          macroinvert_density_top5deep_ar_stacked <- ggplot(dat_macroinvert_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_density_m2*100, x=`Site type`)) + 
            geom_bar(position="stack", stat="identity") +
            labs(x = "Site type", y = expression(paste("Average density per 100m" ^2)), fill = "Macroinvertebrate species") +
            scale_fill_manual(values = pal_macro10) +
            scale_y_continuous(expand = c(0,0)) +
            theme_classic() +
            theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))
          
          ggsave(macroinvert_density_top5deep_ar_stacked, path = file.path("figures"), filename = "macroinvert_density_top5deep_ar_stacked.jpg", height = 3.5, width = 5, units = "in")
          
          #kelp density
          kelp_density_top5deep_ar_stacked <- ggplot(dat_kelp_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_density_m2*100, x=`Site type`)) + 
            geom_bar(position="stack", stat="identity") +
            labs(x = "Site type", y = expression(paste("Average density per 100m" ^2)), fill = "Kelp species") +
            scale_fill_manual(values = pal_kelp8) +
            scale_y_continuous(expand = c(0,0)) +
            theme_classic() +
            theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))
          
          ggsave(kelp_density_top5deep_ar_stacked, path = file.path("figures"), filename = "kelp_density_top5deep_ar_stacked.jpg", height = 3.5, width = 5, units = "in")
          
          #Merge figures
      CRANE_deep_artificial_SoS_2022_23 <- plot_grid(fish_density_top5deep_ar_stacked + theme(axis.title.x = element_blank()),
                    fish_biomass_top5deep_ar_stacked + theme(axis.title.x = element_blank()),
                    macroinvert_density_top5deep_ar_stacked,
                    kelp_density_top5deep_ar_stacked + labs(y=" "))
      
      ggsave(CRANE_deep_artificial_SoS_2022_23, path = file.path("figures"), filename = "CRANE_deep_artificial_SoS_2022_23.jpg",
             height = 8, width = 11, units = "in")
      ggsave(CRANE_deep_artificial_SoS_2022_23, path = file.path("figures"), filename = "CRANE_deep_artificial_SoS_2022_23.pdf",
             height = 8, width = 11, units = "in")
      


##################################################################################
#TO DO: average biomass and density for natural reefs only, MPA vs non-MPA (not sure how to determine this right now)
##################################################################################
