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

source(file.path("functions","add_kruskal_wallis_dunn_test_letter_boxplot.R"))

########################
##Load data
########################
dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))
dat_fish_site_averages <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages.rds"))
dat_macroinvert_site_averages <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages.rds"))
dat_kelp_site_averages <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages.rds"))

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

#New Column Identifying Island vs Mainland
dat_fish_site_averages[,type := factor(ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland"))]
dat_macroinvert_site_averages[,type := factor(ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland"))]
dat_kelp_site_averages[,type := factor(ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland"))]

#New column identifying ar mainland, natural island, and natural mainland
#New Column Identifying ARM vs Island vs Natural Coast
dat_fish_site_averages[,type_wAR := factor(ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island",ifelse(DepthZone %in% c("ARM","Module"),"Artificial_reef","Natural_mainland")))]
dat_macroinvert_site_averages[,type_wAR := factor(ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island",ifelse(DepthZone %in% c("ARM","Module"),"Artificial_reef","Natural_mainland")))]
dat_kelp_site_averages[,type_wAR := factor(ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island",ifelse(DepthZone %in% c("ARM","Module"),"Artificial_reef","Natural_mainland")))]

#Split depth zone for artificial reefs in PV and SM
dat_fish_site_averages[,DepthZone := ifelse(DepthZone %in% c( "ARM", "Module") & AR_Complex == "Palos Verdes Reef","AR_PVR",ifelse(DepthZone %in% c( "ARM", "Module") & AR_Complex != "Palos Verdes", "AR_SM",as.character(DepthZone)))]
dat_kelp_site_averages[,DepthZone := ifelse(DepthZone %in% c( "ARM", "Module") & AR_Complex == "Palos Verdes Reef","AR_PVR",ifelse(DepthZone %in% c( "ARM", "Module") & AR_Complex != "Palos Verdes", "AR_SM",as.character(DepthZone)))]
dat_macroinvert_site_averages[,DepthZone := ifelse(DepthZone %in% c( "ARM", "Module") & AR_Complex == "Palos Verdes Reef","AR_PVR",ifelse(DepthZone %in% c( "ARM", "Module") & AR_Complex != "Palos Verdes", "AR_SM",as.character(DepthZone)))]

#Add in MPAs
MPA_site_key <- readRDS(file.path("keys","MPA_site_key.rds"))
dat_fish_site_averages <- MPA_site_key[dat_fish_site_averages, on = "Site"]
dat_kelp_site_averages <- MPA_site_key[dat_kelp_site_averages, on = "Site"]
dat_macroinvert_site_averages <- MPA_site_key[dat_macroinvert_site_averages, on = "Site"]

########################
##Total abundance (summed density of all taxa, summed biomass of all taxa) by depth zone
########################

#abundance for fish
dat_fish_site_averages[,total_abundance_depthzone_site := sum(mean_density_m2)*100,.(Site, DepthZone)] #multiply by 100 to get # per 100m^2
dat_fish_site_averages[,total_biomass_depthzone_site := sum(mean_wt_density_g_m2)*100/1000,.(Site, DepthZone)] #total fish biomass in kg per 100m^2 (multiply by 100, divide by 1000)

dat_fish_total_abundances <- unique(dat_fish_site_averages[,.(Site, DepthZone, type, type_wAR, total_abundance_depthzone_site, total_biomass_depthzone_site, MPA_overlap)])

#Set factor order for depth zone
dat_fish_total_abundances[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]

#Set factor order for type
dat_fish_total_abundances[,type_wAR := factor(type_wAR, levels = c("Island","Natural_mainland","Artificial_reef"))]

#mean island mainland depth zone
dat_fish_abundance_summary_stats <- dat_fish_total_abundances[,.(mean_dens = mean(total_abundance_depthzone_site),SD_dens = sd(total_abundance_depthzone_site),
                          mean_biomass_dens = mean(total_biomass_depthzone_site),SD_biomass_dens = sd(total_biomass_depthzone_site)),.(DepthZone, type_wAR, MPA_overlap)]

#abundance for macroinverts
dat_macroinvert_site_averages[,total_abundance_depthzone_site := sum(mean_density_m2)*100,.(Site, DepthZone)] #multiply by 100 to get # per 100m^2
dat_macroinvert_total_abundances <- unique(dat_macroinvert_site_averages[,.(Site, DepthZone, type, type_wAR, total_abundance_depthzone_site, MPA_overlap)])

#Set factor order for depth zone
dat_macroinvert_total_abundances[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]

#Set factor order for type
dat_macroinvert_total_abundances[,type_wAR := factor(type_wAR, levels = c("Island","Natural_mainland","Artificial_reef"))]

#abundance for kelp
dat_kelp_site_averages[,total_abundance_depthzone_site := sum(mean_density_m2)*100,.(Site, DepthZone)] #multiply by 100 to get # per 100m^2
dat_kelp_total_abundances <- unique(dat_kelp_site_averages[,.(Site, DepthZone, type, type_wAR, total_abundance_depthzone_site, MPA_overlap)])
#Set factor order for depth zone
dat_kelp_total_abundances[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]

#Set factor order for type
dat_kelp_total_abundances[,type_wAR := factor(type_wAR, levels = c("Island","Natural_mainland","Artificial_reef"))]


#Visualize all fish density ####
##Island vs. mainland ####
fish_abundance_depthzone <- ggplot(dat_fish_total_abundances) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_abundance_depthzone_site, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2,
  ) +
  scale_fill_manual(values = c("#7FB1D3","#BDBAD9","#FB8071"), labels = c("Natural island","Natural mainland", "Artificial mainland")) +
  scale_pattern_manual(values = c("stripe","none","none"), labels = c("Natural island","Natural mainland", "Artificial mainland")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","AR\nPV","AR\nSMB")) +
  labs(x = "", y = bquote("   Density\n(count per 100 m"^2*")") , fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
 ggtitle("Fish")+
   theme(
    legend.position = c(0.97, 0.95),
    legend.justification = c("right", "top"),
    legend.key.size = unit(2,"line"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave(fish_abundance_depthzone, path = file.path("figures"), filename ="fish_abundance_depthzone.jpg", height = 4.5, width = 6, units = "in")

(fish_abundance_depthzone_dunn <- generate_boxplot_nonpar_dunn_letters(dt = dat_fish_total_abundances, var = "total_abundance_depthzone_site", plot = fish_abundance_depthzone))

#Compare fish densities between islands and mainland with Wilcox rank-sum test (Mann-Whitney U test)
wilcox.test(dat_fish_total_abundances[type_wAR == "Natural_mainland"]$total_abundance_depthzone_site, dat_fish_total_abundances[type_wAR == "Island"]$total_abundance_depthzone_site) #Fish density significantly higher at island sites
#The two groups are significantly different (W = 1619, p < 0.0001)

#Means and SD
mean(dat_fish_total_abundances[type_wAR == "Island"]$total_abundance_depthzone_site); sd(dat_fish_total_abundances[type_wAR == "Island"]$total_abundance_depthzone_site)
#Mean = 144, SD = 136
mean(dat_fish_total_abundances[type_wAR == "Natural_mainland"]$total_abundance_depthzone_site); sd(dat_fish_total_abundances[type_wAR == "Natural_mainland"]$total_abundance_depthzone_site)
#Mean = 39, SD = 37

##Island vs. mainland vs. MPA status ####
fish_abundance_depthzone_MPA <- ggplot(dat_fish_total_abundances) +
  geom_boxplot(aes(x = DepthZone, y = total_abundance_depthzone_site, color = MPA_overlap, fill = type_wAR), position = position_dodge2(preserve = "single"),
                       outlier.size = 1.2,
  ) +
  scale_fill_manual(values = c("#7FB1D3","#BDBAD9","#FB8071"), labels = c("Natural island","Natural mainland", "Artificial mainland")) +
  scale_color_manual(values = c("black","darkgrey")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","AR\nPV","AR\nSMB")) +
  labs(x = "", y = bquote("   Density\n(count per 100 m"^2*")") , fill = "Reef type", pattern = "Reef type", color = "Within MPA") +
  theme_classic() +
  ggtitle("Fish")+
  theme(
    legend.position = c(0.97, 0.95),
    legend.justification = c("right", "top"),
    legend.key.size = unit(2,"line"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave(fish_abundance_depthzone_MPA, path = file.path("figures"), filename ="fish_abundance_depthzone_MPA.jpg", height = 4.5, width = 6, units = "in")


#Visualize all fish biomass ####

fish_biomass_depthzone <- ggplot(dat_fish_total_abundances) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_biomass_depthzone_site, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#7FB1D3","#BDBAD9","#FB8071"), labels = c("Natural island","Natural mainland","Artificial reef")) +
  scale_pattern_manual(values = c("stripe","none","none"), labels = c("Natural island","Natural mainland","Artificial reef")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos\nVerdes","Santa\nMonica\nBay")) +
  labs(x = "                        Depth Zone                AR Complex", y = bquote("   Density\n(kg per 100 m"^2*")"), fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.97, 0.95),
    legend.justification = c("right", "top")
  )

ggsave(fish_biomass_depthzone, path = file.path("figures"), filename ="fish_biomass_depthzone.jpg", height = 4.5, width = 6, units = "in")

(fish_biomass_depthzone_dunn <- generate_boxplot_nonpar_dunn_letters(dat_fish_total_abundances, "total_biomass_depthzone_site", fish_biomass_depthzone))

#Just Island vs. Mainland comparison
wilcox.test(dat_fish_total_abundances[type_wAR == "Natural_mainland"]$total_biomass_depthzone_site, dat_fish_total_abundances[type_wAR == "Island"]$total_biomass_depthzone_site) #Fish density significantly higher at island sites
#W = 2628, p < 0.0001 (significantly different)


#Means and SD
mean(dat_fish_total_abundances[type == "Island"]$total_biomass_depthzone_site); sd(dat_fish_total_abundances[type == "Island"]$total_biomass_depthzone_site)
mean(dat_fish_total_abundances[type == "Natural mainland"]$total_biomass_depthzone_site); sd(dat_fish_total_abundances[type == "Natural mainland"]$total_biomass_depthzone_site)

##Island vs. mainland vs. MPA status ####
fish_biomass_depthzone_MPA <- ggplot(dat_fish_total_abundances) +
  geom_boxplot(aes(x = DepthZone, y = total_biomass_depthzone_site, color = MPA_overlap, fill = type_wAR), position = position_dodge2(preserve = "single"),
               outlier.size = 1.2,
  ) +
  scale_fill_manual(values = c("#7FB1D3","#BDBAD9","#FB8071"), labels = c("Natural island","Natural mainland", "Artificial mainland")) +
  scale_color_manual(values = c("black","darkgrey")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","AR\nPV","AR\nSMB")) +
  labs(x = "", y = bquote("   Density\n(kg per 100 m"^2*")"), fill = "Reef type", pattern = "Reef type", color = "Within MPA") +
  theme_classic() +
  ggtitle("Fish")+
  theme(
    legend.position = c(0.97, 0.95),
    legend.justification = c("right", "top"),
    legend.key.size = unit(2,"line"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave(fish_biomass_depthzone_MPA, path = file.path("figures"), filename ="fish_biomass_depthzone_MPA.jpg", height = 4.5, width = 6, units = "in")

#Visualize all macroinvertebrate density####
##Island vs. mainland ####
macroinvert_abundance_depthzone <- ggplot(dat_macroinvert_total_abundances) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_abundance_depthzone_site, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos\nVerdes","Santa\nMonica\nBay")) +
  labs(x = "                        Depth Zone                AR Complex", y = bquote("Density (count per 100 m"^2*")"), fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  ggtitle("Macroinvertebrate") +
  theme(
    legend.position = c(0.85, 0.99),
    legend.justification = c("right", "top"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave(macroinvert_abundance_depthzone, path = file.path("figures"), filename ="macroinvert_abundance_depthzone.jpg", height = 4.5, width = 6, units = "in")

macroinvert_abundance_depthzone_dunn <- generate_boxplot_nonpar_dunn_letters(dat_macroinvert_total_abundances, "total_abundance_depthzone_site", macroinvert_abundance_depthzone)

#Just Island vs. Mainland comparison
wilcox.test(dat_macroinvert_total_abundances[type_wAR == "Natural_mainland"]$total_abundance_depthzone_site, dat_macroinvert_total_abundances[type_wAR == "Island"]$total_abundance_depthzone_site) #macroinvert density significantly higher at island sites

#Means and SD
mean(dat_macroinvert_total_abundances[type_wAR == "Island"]$total_abundance_depthzone_site); sd(dat_macroinvert_total_abundances[type_wAR == "Island"]$total_abundance_depthzone_site)
mean(dat_macroinvert_total_abundances[type_wAR == "Natural_mainland"]$total_abundance_depthzone_site); sd(dat_macroinvert_total_abundances[type_wAR == "Natural_mainland"]$total_abundance_depthzone_site)

##Island vs. mainland vs. MPA status ####
macroinvert_abundance_depthzone_MPA <- ggplot(dat_macroinvert_total_abundances) +
  geom_boxplot(aes(x = DepthZone, y = total_abundance_depthzone_site, color = MPA_overlap, fill = type_wAR), position = position_dodge2(preserve = "single"),
               outlier.size = 1.2,
  ) +
  scale_fill_manual(values = c("#7FB1D3","#BDBAD9","#FB8071"), labels = c("Natural island","Natural mainland", "Artificial mainland")) +
  scale_color_manual(values = c("black","darkgrey")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","AR\nPV","AR\nSMB")) +
  labs(x = "", y = bquote("   Density\n(count per 100 m"^2*")"), fill = "Reef type", pattern = "Reef type", color = "Within MPA") +
  theme_classic() +
  ggtitle("Macroinvertebrate")+
  theme(
    legend.position = c(0.97, 0.95),
    legend.justification = c("right", "top"),
    legend.key.size = unit(2,"line"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave(macroinvert_abundance_depthzone_MPA, path = file.path("figures"), filename ="macroinvert_abundance_depthzone_MPA.jpg", height = 4.5, width = 6, units = "in")

#Visualize all kelp density ####
#Island vs. mainland ####
kelp_abundance_depthzone <- ggplot(dat_kelp_total_abundances) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_abundance_depthzone_site, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos\nVerdes","Santa\nMonica\nBay")) +
  labs(x = "                        Depth Zone                AR Complex", y = bquote("Density (count per 100 m"^2*")"), fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  ggtitle("Macroalgae") +
  theme(
    legend.position = c(0.97, 0.95),
    legend.justification = c("right", "top"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave(kelp_abundance_depthzone, path = file.path("figures"), filename ="kelp_abundance_depthzone.jpg", height = 4.5, width = 6, units = "in")

kelp_abundance_depthzone_dunn <- generate_boxplot_nonpar_dunn_letters(dat_kelp_total_abundances, "total_abundance_depthzone_site", kelp_abundance_depthzone)

#Just Island vs. Mainland comparison
wilcox.test(dat_kelp_total_abundances[type_wAR == "Natural_mainland"]$total_abundance_depthzone_site, dat_kelp_total_abundances[type_wAR == "Island"]$total_abundance_depthzone_site) #kelp density significantly higher at island sites

#Means and SD
mean(dat_kelp_total_abundances[type_wAR == "Island"]$total_abundance_depthzone_site); sd(dat_kelp_total_abundances[type_wAR == "Island"]$total_abundance_depthzone_site)
mean(dat_kelp_total_abundances[type_wAR == "Natural_mainland"]$total_abundance_depthzone_site); sd(dat_kelp_total_abundances[type_wAR == "Natural_mainland"]$total_abundance_depthzone_site)

##Island vs. mainland vs. MPA status ####
kelp_abundance_depthzone_MPA <- ggplot(dat_kelp_total_abundances) +
  geom_boxplot(aes(x = DepthZone, y = total_abundance_depthzone_site, color = MPA_overlap, fill = type_wAR), position = position_dodge2(preserve = "single"),
               outlier.size = 1.2,
  ) +
  scale_fill_manual(values = c("#7FB1D3","#BDBAD9","#FB8071"), labels = c("Natural island","Natural mainland", "Artificial mainland")) +
  scale_color_manual(values = c("black","darkgrey")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","AR\nPV","AR\nSMB")) +
  labs(x = "", y = bquote("   Density\n(count per 100 m"^2*")"), fill = "Reef type", pattern = "Reef type", color = "Within MPA") +
  theme_classic() +
  ggtitle("Macroalgae")+
  theme(
    legend.position = c(0.97, 0.95),
    legend.justification = c("right", "top"),
    legend.key.size = unit(2,"line"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
  )

ggsave(kelp_abundance_depthzone_MPA, path = file.path("figures"), filename ="kelp_abundance_depthzone_MPA.jpg", height = 4.5, width = 6, units = "in")


########################
##Total richness (summed # of taxa, summed biomass of all taxa) by depth zone ####
########################
#Calculate total # of species across all sites
#Fish

dat_fish_site_averages.r <- dat_fish_site_averages[mean_density_m2 >0,] #only rows with values
dat_fish_site_averages.r[,total_richness_depthzone_site := uniqueN(taxa), by= .(Site, DepthZone, type_wAR)]
dat_fish_richness_depthzone_site <- unique(dat_fish_site_averages.r[,.(Site, DepthZone, type_wAR, total_richness_depthzone_site)])
#Set factor order for depth zone
dat_fish_richness_depthzone_site[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]

#Macroinverts
dat_macroinvert_site_averages.r <- dat_macroinvert_site_averages[mean_density_m2 >0,] #only rows with values
dat_macroinvert_site_averages.r[,total_richness_depthzone_site := uniqueN(taxa), by= .(Site, DepthZone, type_wAR)]
dat_macroinvert_richness_depthzone_site <- unique(dat_macroinvert_site_averages.r[,.(Site, DepthZone, type_wAR, total_richness_depthzone_site)])
#Set factor order for depth zone
dat_macroinvert_richness_depthzone_site[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]

#kelp
dat_kelp_site_averages.r <- dat_kelp_site_averages[mean_density_m2 >0,] #only rows with values
dat_kelp_site_averages.r[,total_richness_depthzone_site := uniqueN(taxa), by= .(Site, DepthZone, type_wAR)]
dat_kelp_richness_depthzone_site <- unique(dat_kelp_site_averages.r[,.(Site, DepthZone, type_wAR, total_richness_depthzone_site)])
#Set factor order for depth zone
dat_kelp_richness_depthzone_site[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]



#Visualize all fish richness ####
##Island vs. mainland
fish_richness_depthzone <- ggplot(dat_fish_richness_depthzone_site) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_richness_depthzone_site, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos\nVerdes","Santa\nMonica\nBay")) +
  labs(x = "                        Depth Zone                AR Complex", y = "Richness\n", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.4, 1),
    legend.justification = c("right", "top")
  )

ggsave(fish_richness_depthzone, path = file.path("figures"), filename ="fish_richness_depthzone.jpg", height = 4.5, width = 6, units = "in")

fish_richness_depthzone_dunn <- generate_boxplot_nonpar_dunn_letters(dat_fish_richness_depthzone_site, "total_richness_depthzone_site", fish_richness_depthzone, richness = T)

#Just Island vs. Mainland comparison
wilcox.test(dat_fish_richness_depthzone_site[type_wAR == "Natural_mainland"]$total_richness_depthzone_site, dat_fish_richness_depthzone_site[type_wAR == "Island"]$total_richness_depthzone_site) #Fish richness significantly higher at mainland sites

#Means and SD
mean(dat_fish_richness_depthzone_site[type_wAR == "Island"]$total_richness_depthzone_site); sd(dat_fish_richness_depthzone_site[type_wAR == "Island"]$total_richness_depthzone_site)
mean(dat_fish_richness_depthzone_site[type_wAR == "Natural_mainland"]$total_richness_depthzone_site); sd(dat_fish_richness_depthzone_site[type_wAR == "Natural_mainland"]$total_richness_depthzone_site)

#Visualize all macroinvertebrate richness ####
##Island vs. mainland ####
macroinvert_richness_depthzone <- ggplot(dat_macroinvert_richness_depthzone_site) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_richness_depthzone_site, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos\nVerdes","Santa\nMonica\nBay")) +
  labs(x = "                        Depth Zone                AR Complex", y = "Richness", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.4, 1),
    legend.justification = c("right", "top")
  )

ggsave(macroinvert_richness_depthzone, path = file.path("figures"), filename ="macroinvert_richness_depthzone.jpg", height = 4.5, width = 6, units = "in")

macroinvert_richness_depthzone_dunn <- generate_boxplot_nonpar_dunn_letters(dat_macroinvert_richness_depthzone_site, "total_richness_depthzone_site", macroinvert_richness_depthzone, richness = T)

#Just Island vs. Mainland comparison
wilcox.test(dat_macroinvert_richness_depthzone_site[type_wAR == "Natural_mainland"]$total_richness_depthzone_site, dat_macroinvert_richness_depthzone_site[type_wAR == "Island"]$total_richness_depthzone_site) #macroinvert density significantly higher at island sites

#Means and SD
mean(dat_macroinvert_richness_depthzone_site[type_wAR == "Island"]$total_richness_depthzone_site); sd(dat_macroinvert_richness_depthzone_site[type_wAR == "Island"]$total_richness_depthzone_site)
mean(dat_macroinvert_richness_depthzone_site[type_wAR == "Natural_mainland"]$total_richness_depthzone_site); sd(dat_macroinvert_richness_depthzone_site[type_wAR == "Natural_mainland"]$total_richness_depthzone_site)


#Visualize all kelp richness ####
##Island vs. mainland ####
kelp_richness_depthzone <- ggplot(dat_kelp_richness_depthzone_site) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_richness_depthzone_site, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos\nVerdes","Santa\nMonica\nBay")) +
  labs(x = "                        Depth Zone                AR Complex", y = "Richness", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.95, 1),
    legend.justification = c("right", "top")
  )

ggsave(kelp_richness_depthzone, path = file.path("figures"), filename ="kelp_richness_depthzone.jpg", height = 4.5, width = 6, units = "in")


kelp_richness_depthzone_dunn <- generate_boxplot_nonpar_dunn_letters(dat_kelp_richness_depthzone_site, "total_richness_depthzone_site", kelp_richness_depthzone, richness = T)

#Just Island vs. Mainland comparison
wilcox.test(dat_kelp_richness_depthzone_site[type_wAR == "Natural_mainland"]$total_richness_depthzone_site, dat_kelp_richness_depthzone_site[type_wAR == "Island"]$total_richness_depthzone_site) 

#Means and SD
mean(dat_kelp_richness_depthzone_site[type_wAR == "Island"]$total_richness_depthzone_site); sd(dat_kelp_richness_depthzone_site[type_wAR == "Island"]$total_richness_depthzone_site)
mean(dat_kelp_richness_depthzone_site[type_wAR == "Natural_mainland"]$total_richness_depthzone_site); sd(dat_kelp_richness_depthzone_site[type_wAR == "Natural_mainland"]$total_richness_depthzone_site)


########################
#Simpson diversity index ####
########################
##Manually calculate index ####
#first, take proportion of each species
#then, square the proportions and sum these for each site
#fish
dat_fish_site_averages.r[,total_site_depthzone_density := sum(mean_density_m2),.(Site, DepthZone)]
dat_fish_site_averages.r[,total_site_depthzone_biomass := sum(mean_wt_density_g_m2),.(Site, DepthZone)]
dat_fish_site_averages.r[,taxa_proportion_density := mean_density_m2/total_site_depthzone_density]
dat_fish_site_averages.r[,taxa_proportion_biomass := mean_wt_density_g_m2/total_site_depthzone_biomass]
dat_fish_site_averages.r[,Simpson_index_density := 1-sum(taxa_proportion_density^2),.(Site,DepthZone)]
dat_fish_site_averages.r[,Simpson_index_biomass := 1-sum(taxa_proportion_biomass^2),.(Site,DepthZone)]

dat_fish_simpson <- unique(dat_fish_site_averages.r[,.(Site, DepthZone, type_wAR, Simpson_index_density, Simpson_index_biomass)])
#Set factor order for depth zones and AR complexes
dat_fish_simpson[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]

#macroinvert
dat_macroinvert_site_averages.r[,total_site_depthzone_density := sum(mean_density_m2),.(Site, DepthZone)]
dat_macroinvert_site_averages.r[,taxa_proportion_density := mean_density_m2/total_site_depthzone_density]
dat_macroinvert_site_averages.r[,Simpson_index_density := 1-sum(taxa_proportion_density^2),.(Site,DepthZone)]

dat_macroinvert_simpson <- unique(dat_macroinvert_site_averages.r[,.(Site, DepthZone, type_wAR, Simpson_index_density)])
#Set factor order for depth zones and AR complexes
dat_macroinvert_simpson[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]

#macroalgae
dat_kelp_site_averages.r[,total_site_depthzone_density := sum(mean_density_m2),.(Site, DepthZone)]
dat_kelp_site_averages.r[,taxa_proportion_density := mean_density_m2/total_site_depthzone_density]
dat_kelp_site_averages.r[,Simpson_index_density := 1-sum(taxa_proportion_density^2),.(Site,DepthZone)]

dat_kelp_simpson <- unique(dat_kelp_site_averages.r[,.(Site, DepthZone, type_wAR, Simpson_index_density)])
#Set factor order for depth zones and AR complexes
dat_kelp_simpson[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]

#Visualize fish diversity (calculated with biomass) ####
##Island vs. mainland
fish_simpson_depthzone_density <- ggplot(dat_fish_simpson) +
  geom_boxplot_pattern(aes(x = DepthZone, y = Simpson_index_density, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos\nVerdes","Santa\nMonica\nBay")) +
  labs(x = "                        Depth Zone                AR Complex",  y = "Simpson index\n(count-based)", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.87, 0.15)
  )

ggsave(fish_simpson_depthzone_density, path = file.path("figures"), filename ="fish_simpson_depthzone_density.jpg", height = 4.5, width = 6, units = "in")

fish_simpson_depthzone_density_dunn <- generate_boxplot_nonpar_dunn_letters(dat_fish_simpson, "Simpson_index_density", fish_simpson_depthzone_density)

#Just Island vs. Mainland comparison
wilcox.test(dat_fish_simpson[type_wAR == "Natural_mainland"]$Simpson_index_density, dat_fish_simpson[type_wAR == "Island"]$Simpson_index_density) #Fish density significantly higher at island sites

#Means and SD
mean(dat_fish_simpson[type_wAR == "Island"]$Simpson_index_density); sd(dat_fish_simpson[type_wAR == "Island"]$Simpson_index_density)
mean(dat_fish_simpson[type_wAR == "Natural_mainland"]$Simpson_index_density); sd(dat_fish_simpson[type_wAR == "Natural_mainland"]$Simpson_index_density)

#same but for fish biomass

fish_simpson_depthzone_biomass <- ggplot(dat_fish_simpson) +
  geom_boxplot_pattern(aes(x = DepthZone, y = Simpson_index_biomass, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos\nVerdes","Santa\nMonica\nBay")) +
  labs(x = "                        Depth Zone                AR Complex", y = "Simpson index\n(biomass-based)", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.87, 0.15)
  )

ggsave(fish_simpson_depthzone_biomass, path = file.path("figures"), filename ="fish_simpson_depthzone_biomass.jpg", height = 4.5, width = 6, units = "in")

summary(lm(data = dat_fish_simpson[type_wAR != "Artificial_reef"], Simpson_index_biomass~DepthZone))

#When using biomass instead of abundance to calculate diversity, diversity is much more consistent across Depth Zones
#and exhibit little difference between island and mainland reefs, largely due to the small size of high abundance species such as blacksmith and gobies

summary(lm(data = dat_fish_simpson[type_wAR != "Artificial_reef"], Simpson_index_biomass~DepthZone*type_wAR))

fish_simpson_depthzone_biomass_dunn <- generate_boxplot_nonpar_dunn_letters(dat_fish_simpson, "Simpson_index_biomass", fish_simpson_depthzone_biomass)

#Just Island vs. Mainland comparison
wilcox.test(dat_fish_simpson[type_wAR == "Natural_mainland"]$Simpson_index_biomass, dat_fish_simpson[type_wAR == "Island"]$Simpson_index_biomass) #Fish biomass significantly higher at island sites

#Means and SD
mean(dat_fish_simpson[type_wAR == "Island"]$Simpson_index_biomass); sd(dat_fish_simpson[type_wAR == "Island"]$Simpson_index_biomass)
mean(dat_fish_simpson[type_wAR == "Natural_mainland"]$Simpson_index_biomass); sd(dat_fish_simpson[type_wAR == "Natural_mainland"]$Simpson_index_biomass)


#visualize macroinvert diversity (calculated with density)
macroinvert_simpson_depthzone_density <- ggplot(dat_macroinvert_simpson) +
  geom_boxplot_pattern(aes(x = DepthZone, y = Simpson_index_density, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial mainland","Natural island","Natural mainland")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos\nVerdes","Santa\nMonica\nBay")) +
  labs(x = "                        Depth Zone                AR Complex", y = "Simpson index", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.87, 0.15)
  )

macroinvert_simpson_depthzone_dunn <- generate_boxplot_nonpar_dunn_letters(dat_macroinvert_simpson, "Simpson_index_density", macroinvert_simpson_depthzone_density)

#Just Island vs. Mainland comparison
wilcox.test(dat_macroinvert_simpson[type_wAR == "Natural_mainland"]$Simpson_index_density, dat_macroinvert_simpson[type_wAR == "Island"]$Simpson_index_density) #macroinvert density significantly higher at island sites

#Means and SD
mean(dat_macroinvert_simpson[type_wAR == "Island"]$Simpson_index_density); sd(dat_macroinvert_simpson[type_wAR == "Island"]$Simpson_index_density)
mean(dat_macroinvert_simpson[type_wAR == "Natural_mainland"]$Simpson_index_density); sd(dat_macroinvert_simpson[type_wAR == "Natural_mainland"]$Simpson_index_density)


#visualize kelp diversity (calculated with density)
kelp_simpson_depthzone_density <- ggplot(dat_kelp_simpson) +
  geom_boxplot_pattern(aes(x = DepthZone, y = Simpson_index_density, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial mainland","Natural island","Natural_mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial mainland","Natural island","Natural_mainland")) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos\nVerdes","Santa\nMonica\nBay")) +
  labs(x = "                        Depth Zone                AR Complex", y = "Simpson index", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.87, 0.85)
  )

ggsave(kelp_simpson_depthzone_density, path = file.path("figures"), filename ="kelp_simpson_depthzone_density.jpg", height = 4.5, width = 6, units = "in")

kelp_simpson_depthzone_density_dunn <- generate_boxplot_nonpar_dunn_letters(dat_kelp_simpson, "Simpson_index_density", kelp_simpson_depthzone_density)

#Just Island vs. Mainland comparison
wilcox.test(dat_kelp_simpson[type_wAR == "Natural_mainland"]$Simpson_index_density, dat_kelp_simpson[type_wAR == "Island"]$Simpson_index_density) #kelp density significantly higher at island sites

#Means and SD
mean(dat_kelp_simpson[type_wAR == "Island"]$Simpson_index_density); sd(dat_kelp_simpson[type_wAR == "Island"]$Simpson_index_density)
mean(dat_kelp_simpson[type_wAR == "Natural_mainland"]$Simpson_index_density); sd(dat_kelp_simpson[type_wAR == "Natural_mainland"]$Simpson_index_density)


#add legend on top
megamerge_dunn_leg <- get_legend(fish_simpson_depthzone_density_dunn + theme(legend.position = "top", legend.direction = "horizontal"))

#Plot all using patchwork package
# Arrange the plots
depthzone_abun_rich_div_merge_dunn <- (fish_abundance_depthzone_dunn + theme(axis.title.x = element_blank(), legend.position = "none") + 
                   kelp_abundance_depthzone_dunn + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
                   macroinvert_abundance_depthzone_dunn + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")) /      # Row 1: 3 plots across
  (fish_biomass_depthzone_dunn + theme(axis.title.x = element_blank(), legend.position = "none") + plot_spacer() +plot_spacer())/ # Row 2: 1 plot taking 1/3 of the width
  (fish_richness_depthzone_dunn + theme(axis.title.x = element_blank(), legend.position = "none") +
     kelp_richness_depthzone_dunn + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none") +
     macroinvert_richness_depthzone_dunn + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")) /      # Row 3: 3 plots across
  (fish_simpson_depthzone_density_dunn + theme(axis.title.x = element_blank(), legend.position = "none") +
     kelp_simpson_depthzone_density_dunn + theme(axis.title.y = element_blank(), legend.position = "none") +
     macroinvert_simpson_depthzone_dunn + theme(axis.title.y = element_blank(), legend.position = "none"))/      # Row 4: 3 plots across
  (fish_simpson_depthzone_biomass_dunn + theme(legend.position = "none") + plot_spacer() +plot_spacer()) # Row 5: 1 plot taking 1/3 of the width

depthzone_abun_rich_div_merge_dunn <- depthzone_abun_rich_div_merge_dunn + plot_annotation(tag_levels = 'a') + plot_layout(guides = "collect")

# Combine the legend grob and the arranged plots
depthzone_abun_rich_div_merge_dunn.l <- plot_grid(megamerge_dunn_leg, depthzone_abun_rich_div_merge_dunn, ncol = 1, rel_heights = c(1,10))

#save
ggsave(depthzone_abun_rich_div_merge_dunn.l, path = file.path("figures"), filename = "depthzone_abun_rich_div_merge_dunn.l.jpg", height = 14, width = 12, unit = "in")


#########
#Instead, plot by transforming metrics from wide to long in a single DT ####
#as an alternate route for plotting (facet), try to merge all data tables and then flip from wide to long

##Link together metrics across taxa ####
#abundance/density
dat_fish_total_abundances[, taxa_type := "Fish"]
dat_fish_total_abundances.l <- melt(dat_fish_total_abundances, id.vars = c("taxa_type","Site","DepthZone","type_wAR"),measure.vars = c("total_abundance_depthzone_site","total_biomass_depthzone_site"), variable.name = "metric",value.name = "value")
dat_macroinvert_total_abundances[, taxa_type := "Macroinvertebrate"]
dat_macroinvert_total_abundances.l <- melt(dat_macroinvert_total_abundances, id.vars = c("taxa_type","Site","DepthZone","type_wAR"),measure.vars = c("total_abundance_depthzone_site"), variable.name = "metric",value.name = "value")
dat_kelp_total_abundances[, taxa_type := "Macroalgae"]
dat_kelp_total_abundances.l <- melt(dat_kelp_total_abundances, id.vars = c("taxa_type","Site","DepthZone","type_wAR"), variable.name = "metric",measure.vars = c("total_abundance_depthzone_site"), value.name = "value")


#richness
dat_fish_richness_depthzone_site[,taxa_type := "Fish"]
dat_fish_richness_depthzone_site.l <- melt(dat_fish_richness_depthzone_site, id.vars = c("taxa_type","Site","DepthZone","type_wAR"),measure.vars = c("total_richness_depthzone_site"), variable.name = "metric",value.name = "value")
dat_macroinvert_richness_depthzone_site[,taxa_type := "Macroinvertebrate"]
dat_macroinvert_richness_depthzone_site.l <- melt(dat_macroinvert_richness_depthzone_site, id.vars = c("taxa_type","Site","DepthZone","type_wAR"),measure.vars = c("total_richness_depthzone_site"), variable.name = "metric",value.name = "value")
dat_kelp_richness_depthzone_site[,taxa_type := "Macroalgae"]
dat_kelp_richness_depthzone_site.l <- melt(dat_kelp_richness_depthzone_site, id.vars = c("taxa_type","Site","DepthZone","type_wAR"),measure.vars = c("total_richness_depthzone_site"), variable.name = "metric",value.name = "value")

#simpson diversity index
dat_fish_simpson[,taxa_type := "Fish"]
dat_fish_simpson.l <- melt(dat_fish_simpson, id.vars = c("taxa_type","Site","DepthZone","type_wAR"),measure.vars = c("Simpson_index_density","Simpson_index_biomass"), variable.name = "metric",value.name = "value")
dat_macroinvert_simpson[,taxa_type := "Macroinvertebrate"]
dat_macroinvert_simpson.l <- melt(dat_macroinvert_simpson, id.vars = c("taxa_type","Site","DepthZone","type_wAR"),measure.vars = c("Simpson_index_density"), variable.name = "metric",value.name = "value")
dat_kelp_simpson[,taxa_type := "Macroalgae"]
dat_kelp_simpson.l <- melt(dat_kelp_simpson, id.vars = c("taxa_type","Site","DepthZone","type_wAR"),measure.vars = c("Simpson_index_density"), variable.name = "metric",value.name = "value")

#rbind all above

merged_abundance_richness_diversity_datatables <- rbind(dat_fish_total_abundances.l, dat_macroinvert_total_abundances.l, dat_kelp_total_abundances.l,
                                                        dat_fish_richness_depthzone_site.l, dat_macroinvert_richness_depthzone_site.l, dat_kelp_richness_depthzone_site.l,
                                                        dat_fish_simpson.l, dat_macroinvert_simpson.l, dat_kelp_simpson.l)

# Custom labels for the 'metric' variable
metric_labels <- c(
  total_abundance_depthzone_site = "Density\n(count-based)",
  total_biomass_depthzone_site = "Density\n(biomass-based)",
  total_richness_depthzone_site = "Richness\n ",
  Simpson_index_density = "Simpson diversity\n(count-based)",
  Simpson_index_biomass = "Simpson diversity\n(biomass-based)"
)

#change ARM to AR only
merged_abundance_richness_diversity_datatables[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"), labels = c("Inner","Middle","Outer","Deep","AR\nPV", "AR\nSMB"))]

# Convert 'metric' column to factor
merged_abundance_richness_diversity_datatables[, metric := factor(metric, levels = names(metric_labels))]

#Edit type_wAR, and set order
merged_abundance_richness_diversity_datatables[, type_wAR := ifelse(type_wAR == "Island","Island",ifelse(DepthZone %in% c("AR\nPV", "AR\nSMB"),"Artificial mainland","Natural_mainland"))][,type_wAR := factor(type_wAR, c("Island","Natural_mainland","Artificial mainland"))]


#merged plot
depthzone_density_abundance_richness <- ggplot(merged_abundance_richness_diversity_datatables) +
  geom_boxplot_pattern(aes(x = DepthZone, y = value, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#7FB1D3","#BDBAD9","#FB8071"), labels = c("Natural island","Natural mainland","Artificial mainland")) +
  scale_pattern_manual(values = c("stripe","none","none"), labels = c("Natural island","Natural mainland","Artificial mainland")) +
  labs(x = "Depth Zone", y = "Value", fill = "Reef type", pattern = "Reef type") +
  facet_grid(rows = vars(metric),cols = vars(taxa_type), scales = "free", labeller = labeller(metric = as_labeller(metric_labels)), switch = "y") +
  theme_classic() +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  #coord_flip() 
  theme(
    #legend.position = c(0.67, 0.7),
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    axis.title.x = element_blank(),
    legend.key.size = unit(2.5,"lines"),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 13),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_blank()
  )
  

depthzone_density_abundance_richness

ggsave(depthzone_density_abundance_richness, path = "figures",filename = "depthzone_density_abundance_richness.jpg", height = 12, width = 10, units = "in")

#Merged plot, but don't split by depth zone

mainlandvsisland_density_abundance_richness <- ggplot(merged_abundance_richness_diversity_datatables) +
  geom_boxplot_pattern(aes(x = taxa_type, y = value, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#7FB1D3","#BDBAD9","#FB8071"), labels = c("Natural island","Natural mainland","Artificial mainland")) +
  scale_pattern_manual(values = c("stripe","none","none"), labels = c("Natural island","Natural mainland","Artificial mainland")) +
  labs(x = "Taxa", y = "Value", fill = "Reef type", pattern = "Reef type") +
  facet_wrap(~metric, scales = "free", labeller = labeller(metric = as_labeller(metric_labels)), ncol = 5) +
  theme_classic() +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  #coord_flip() 
  theme(
    #legend.position = c(0.67, 0.7),
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    axis.title.x = element_blank(),
    legend.key.size = unit(2.5,"lines"),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 13),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.title.y = element_blank()
  )

ggsave(mainlandvsisland_density_abundance_richness, path = "figures",filename = "mainlandvsisland_density_abundance_richness.jpg", height = 6, width = 12, units = "in")


#################################################################################
#How do the predictors of abundance, richness, and diversity vary by depth zone ####
#################################################################################

#pull environmental data
all_env_lat_lon <- fread(file.path("data","enviro_predictors","all_env_lat_lon.csv"))

#Adjust depthzones to have AR names
all_env_lat_lon[,DepthZone := ifelse(DepthZone %in% c( "ARM", "Module") & grepl("PVR",Site) == T,"AR_PVR",ifelse(DepthZone %in% c( "ARM", "Module") & grepl("PVR",Site) != T, "AR_SM",as.character(DepthZone)))]


#make list of predictors
predictors <- c(
  #"Latitude",
                "dist_200m_bath",
                "min_sst_C",
                "mean_sst_C",
                "max_sst_C",
                "min_chl_mg_m3",
                "mean_chl_mg_m3",
                "max_chl_mg_m3",
                "giantkelp_stipe_density_m2",
                "Relief_index",
                "Relief_SD",
                "Substrate_index",
                "Substrate_SD")


#link abundance/richness/diversity calculations with environmental data for each site
#abundance/density
dat_fish_total_abundances.env <- all_env_lat_lon[dat_fish_total_abundances, on = c("Site","DepthZone")]
dat_fish_total_abundances.env <- dat_fish_total_abundances.env[complete.cases(dat_fish_total_abundances.env)]
dat_fish_total_abundances.env[, (paste0(predictors, ".s")) := lapply(.SD, scale), .SDcols = predictors]
dat_macroinvert_total_abundances.env <- all_env_lat_lon[dat_macroinvert_total_abundances, on = c("Site","DepthZone")]
dat_macroinvert_total_abundances.env <- dat_macroinvert_total_abundances.env[complete.cases(dat_macroinvert_total_abundances.env)]
dat_macroinvert_total_abundances.env[, (paste0(predictors, ".s")) := lapply(.SD, scale), .SDcols = predictors]
dat_kelp_total_abundances.env <- all_env_lat_lon[dat_kelp_total_abundances, on = c("Site","DepthZone")]
dat_kelp_total_abundances.env <- dat_kelp_total_abundances.env[complete.cases(dat_kelp_total_abundances.env)]
dat_kelp_total_abundances.env[, (paste0(predictors, ".s")) := lapply(.SD, scale), .SDcols = predictors]

#richness
dat_fish_richness_depthzone_site.env <- all_env_lat_lon[dat_fish_richness_depthzone_site, on = c("Site","DepthZone")]
dat_fish_richness_depthzone_site.env <- dat_fish_richness_depthzone_site.env[complete.cases(dat_fish_richness_depthzone_site.env)]
dat_fish_richness_depthzone_site.env[, (paste0(predictors, ".s")) := lapply(.SD, scale), .SDcols = predictors]
dat_macroinvert_richness_depthzone_site.env <- all_env_lat_lon[dat_macroinvert_richness_depthzone_site, on = c("Site","DepthZone")]
dat_macroinvert_richness_depthzone_site.env <- dat_macroinvert_richness_depthzone_site.env[complete.cases(dat_macroinvert_richness_depthzone_site.env)]
dat_macroinvert_richness_depthzone_site.env[, (paste0(predictors, ".s")) := lapply(.SD, scale), .SDcols = predictors]
dat_kelp_richness_depthzone_site.env <- all_env_lat_lon[dat_kelp_richness_depthzone_site, on = c("Site","DepthZone")]
dat_kelp_richness_depthzone_site.env <- dat_kelp_richness_depthzone_site.env[complete.cases(dat_kelp_richness_depthzone_site.env)]
dat_kelp_richness_depthzone_site.env[, (paste0(predictors, ".s")) := lapply(.SD, scale), .SDcols = predictors]

#simpson diversity index
dat_fish_simpson.env <- all_env_lat_lon[dat_fish_simpson, on = c("Site","DepthZone")]
dat_fish_simpson.env[, (paste0(predictors, ".s")) := lapply(.SD, scale), .SDcols = predictors]
dat_fish_simpson.env <- dat_fish_simpson.env[complete.cases(dat_fish_simpson.env)]
dat_macroinvert_simpson.env <- all_env_lat_lon[dat_macroinvert_simpson, on = c("Site","DepthZone")]
dat_macroinvert_simpson.env <- dat_macroinvert_simpson.env[complete.cases(dat_macroinvert_simpson.env)]
dat_macroinvert_simpson.env[, (paste0(predictors, ".s")) := lapply(.SD, scale), .SDcols = predictors]
dat_kelp_simpson.env <- all_env_lat_lon[dat_kelp_simpson, on = c("Site","DepthZone")]
dat_kelp_simpson.env <- dat_kelp_simpson.env[complete.cases(dat_kelp_simpson.env)]
dat_kelp_simpson.env[, (paste0(predictors, ".s")) := lapply(.SD, scale), .SDcols = predictors]

#How does substrate and relief vary with depth?
Env_Site <- unique(dat_fish_total_abundances.env[,.(Site, Region, DepthZone, type, type_wAR, Relief_index, Relief_SD, Substrate_index, Substrate_SD)])

#Wide to long
Env_Site.l <- melt(Env_Site, id.vars = c("Site", "Region", "DepthZone", "type", "type_wAR"), variable.name = "Reef metric",value.name = "value")

#Change factor order
Env_Site.l[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"),
                                           labels = c("Inner","Middle","Outer","Deep","AR\nPV","AR\nSMB"))]

#Facet labels
facet_labels <- c("Relief_index" = "a.                                Relief index",
                  "Relief_SD" = "b.                                Relief SD",
                  "Substrate_index" = "c.                                Substrate index",
                  "Substrate_SD" = "d.                                Substrate SD")

#Make plot of relief and substrate for island and mainland reefs ####
#merged plot
depthzone_relief_substrate <- ggplot(Env_Site.l) +
  geom_boxplot_pattern(aes(x = DepthZone, y = value, fill = type_wAR, pattern = type_wAR), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#7FB1D3","#BDBAD9","#FB8071"), labels = c("Natural island","Natural mainland","Artificial mainland")) +
  scale_pattern_manual(values = c("stripe","none","none"), labels = c("Natural island","Natural mainland","Artificial mainland")) +
  labs(x = "Depth Zone/Reef complex", y = "", fill = "Reef type", pattern = "Reef type") +
  facet_wrap(~`Reef metric`, scales = "free", labeller = as_labeller(facet_labels)) +
  theme_classic() +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  #coord_flip() 
  theme(
    #legend.position = c(0.67, 0.7),
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 13, hjust = 0),
    legend.key.size = unit(2.5,"lines"),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 13),
    axis.text.x = element_text(size = 10)
  )


depthzone_relief_substrate 

ggsave(depthzone_relief_substrate, path = "figures",filename = "depthzone_relief_substrate.jpg", height = 8, width = 10, units = "in")

#####################
#Function to make Kendall correlation plots and data tables of correlations for each taxa and metric
make_kendall_cor_plot <- function(taxa, metric_input, metric, metric_name, data.table){
  
  data.table <- data.table(data.table)
  
  metric_data <- data.table[[metric_name]]
  
  data.table[,metric_named := metric_data]
  
  
  #deep
  data.table.deep <- data.table[DepthZone == "Deep",.(metric_named,Latitude, dist_200m_bath, min_sst_C, mean_sst_C, max_sst_C, min_chl_mg_m3, mean_chl_mg_m3, max_chl_mg_m3, giantkelp_stipe_density_m2, Relief_index, Relief_SD, Substrate_index, Substrate_SD)]
  
  cor_d <- data.table(DepthZone = "Deep", predictor = predictors, estimate = numeric(), p_value = numeric())
  for(i in 1:length(predictors)){
    test_correlation <- cor.test(y = data.table.deep$metric_named, x = data.table.deep[,get(predictors[i])], method = "kendall")
    rho <- test_correlation$estimate
    pvalue <- test_correlation$p.value
    cor_d[i,"estimate"] <- rho
    cor_d[i,"p_value"] <- pvalue
    
  }
  
  
  #outer
  data.table.outer <- data.table[DepthZone == "Outer",.(metric_named,Latitude, dist_200m_bath, min_sst_C, mean_sst_C, max_sst_C, min_chl_mg_m3, mean_chl_mg_m3, max_chl_mg_m3, giantkelp_stipe_density_m2, Relief_index, Relief_SD, Substrate_index, Substrate_SD)]
  
  
  cor_o <- data.table(DepthZone = "Outer", predictor = predictors, estimate = numeric(), p_value = numeric())
  for(i in 1:length(predictors)){
    test_correlation <- cor.test(y = data.table.outer$metric_named, x = data.table.outer[,get(predictors[i])], method = "kendall")
    rho <- test_correlation$estimate
    pvalue <- test_correlation$p.value
    cor_o[i,"estimate"] <- rho
    cor_o[i,"p_value"] <- pvalue
    
  }
  
  
  #middle
  data.table.middle <- data.table[DepthZone == "Middle",.(metric_named,Latitude, dist_200m_bath, min_sst_C, mean_sst_C, max_sst_C, min_chl_mg_m3, mean_chl_mg_m3, max_chl_mg_m3, giantkelp_stipe_density_m2, Relief_index, Relief_SD, Substrate_index, Substrate_SD)]
  
  
  cor_m <- data.table(DepthZone = "Middle", predictor = predictors, estimate = numeric(), p_value = numeric())
  for(i in 1:length(predictors)){
    test_correlation <- cor.test(y = data.table.middle$metric_named, x = data.table.middle[,get(predictors[i])], method = "kendall")
    rho <- test_correlation$estimate
    pvalue <- test_correlation$p.value
    cor_m[i,"estimate"] <- rho
    cor_m[i,"p_value"] <- pvalue
    
  }
  
  
  #inner
  data.table.inner <- data.table[DepthZone == "Inner",.(metric_named,Latitude, dist_200m_bath, min_sst_C, mean_sst_C, max_sst_C, min_chl_mg_m3, mean_chl_mg_m3, max_chl_mg_m3, giantkelp_stipe_density_m2, Relief_index, Relief_SD, Substrate_index, Substrate_SD)]
  
  cor_i <- data.table(DepthZone = "Inner", predictor = predictors, estimate = numeric(), p_value = numeric())
  for(i in 1:length(predictors)){
    test_correlation <- cor.test(y = data.table.inner$metric_named, x = data.table.inner[,get(predictors[i])], method = "kendall")
    rho <- test_correlation$estimate
    pvalue <- test_correlation$p.value
    cor_i[i,"estimate"] <- rho
    cor_i[i,"p_value"] <- pvalue
    
  }
  
  #ARM
  data.table.ARM <- data.table[DepthZone == "ARM",.(metric_named,Latitude, dist_200m_bath, min_sst_C, mean_sst_C, max_sst_C, min_chl_mg_m3, mean_chl_mg_m3, max_chl_mg_m3, giantkelp_stipe_density_m2, Relief_index, Relief_SD, Substrate_index, Substrate_SD)]
  
  
  cor_ARM <- data.table(DepthZone = "ARM", predictor = predictors, estimate = numeric(), p_value = numeric())
  for(i in 1:length(predictors)){
    test_correlation <- cor.test(y = data.table.ARM$metric_named, x = data.table.ARM[,get(predictors[i])], method = "kendall")
    rho <- test_correlation$estimate
    pvalue <- test_correlation$p.value
    cor_ARM[i,"estimate"] <- rho
    cor_ARM[i,"p_value"] <- pvalue
    
  }
  
  #merge
  cor_merge <- rbind(cor_d, cor_o,cor_m, cor_i, cor_ARM)
  
  #Plot
  
  #flip wide to long
  #add identifying columns
  cor_merge[,taxa  := taxa][,metric_input := metric_input][,metric := metric][,metric_name := metric_name]
  cor_merge.l <- melt(cor_merge, id.vars = c("taxa","metric_input","metric", "predictor","metric_name","DepthZone"))
  cor_merge.l <- dcast(cor_merge.l, formula = taxa + metric_input + metric + predictor + metric_name + DepthZone ~ variable, value.var = "value")
  
  cor_merge.l[,DepthZone := factor(DepthZone, levels = c("ARM","Deep","Outer","Middle","Inner"), labels = c("AR","Deep","Outer","Middle","Inner"))]
  
  #significance
  cor_merge.l[,significant := ifelse(p_value<0.05,TRUE,FALSE)]
  
  #alter estimates to match significance
  cor_merge.l[,estimate_sig := ifelse(significant==T,estimate,0)]
  
  
  #add better labels and order for environmental variables
cor_merge.l[,predictor := factor(predictor, levels = c(
                                                 "Latitude", "dist_200m_bath", "min_sst_C","mean_sst_C", "max_sst_C","min_chl_mg_m3" ,
                                                 "mean_chl_mg_m3" ,"max_chl_mg_m3","giantkelp_stipe_density_m2", "Relief_index", "Relief_SD","Substrate_index" ,
                                                 "Substrate_SD"),
                                    labels = c("Latitude", "Distance to 200 m isobath", "Minimum SST","Mean SST", "Maximum SST","Minimum chlorophyll" ,
                                               "Mean chlorophyll" ,"Maximum chlorophyll","Giant kelp stipe density", "Relief index", "Relief SD","Substrate index" ,
                                               "Substrate SD"))]

plot_title <- paste(taxa, metric_input, metric, sep = " ")  

  single_plot <- ggplot(cor_merge.l) +
    geom_tile(aes(y = DepthZone, x = predictor, fill = estimate_sig)) +
    geom_text(aes(y = DepthZone, x = predictor, label = round(estimate_sig,2)), size = 2) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1,1)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Habitat characteristic", y = "Depth zone", fill = "Kendall's\ncorrelation\ncoefficient", title = plot_title) +
    geom_segment(aes(x=0.5,xend = 13.5, y = 1.5, yend = 1.5)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = -0.1), axis.title.x = element_blank(), plot.title = element_text(face = "bold"))

plot_name <- paste0(taxa,"_",metric_input,"_", metric,"_Kendall_corr_heatmap")

ggsave(path = file.path("figures","kendall_corr"), filename = paste0(plot_name,".jpg"), height = 5, width = 5, unit = "in")

return(cor_merge.l)
}

#fish
   #fish abundance density
fish_abundance_abundensity_corr.l <- make_kendall_cor_plot(taxa = "Fish",metric_input = "Count density",metric = "Abundance",metric_name = "total_abundance_depthzone_site",data.table = dat_fish_total_abundances.env)
   #fish abundance biomass
fish_abundance_biodensity_corr.l <- make_kendall_cor_plot(taxa = "Fish",metric_input = "Biomass density",metric = "Abundance",metric_name = "total_biomass_depthzone_site",data.table = dat_fish_total_abundances.env)
   #fish richness
fish_richness_corr.l <- make_kendall_cor_plot(taxa = "Fish",metric_input = "",metric = "Species richness",metric_name = "total_richness_depthzone_site",data.table = dat_fish_richness_depthzone_site.env)
   #fish diversity density
fish_simpson_abundensity_corr.l <- make_kendall_cor_plot(taxa = "Fish",metric_input = "Count density",metric = "Simpson diversity",metric_name = "Simpson_index_density",data.table = dat_fish_simpson.env)
   #fish diversity biomass
fish_simpson_biodensity_corr.l <- make_kendall_cor_plot(taxa = "Fish",metric_input = "Biomass density",metric = "Simpson diversity",metric_name = "Simpson_index_biomass",data.table = dat_fish_simpson.env)

#macroinvertebrates
   #macroinvert abundance density
macroinvert_abundance_abundensity_corr.l <- make_kendall_cor_plot(taxa = "Macroinvertebrate",metric_input = "Count density",metric = "Abundance",metric_name = "total_abundance_depthzone_site",data.table = dat_macroinvert_total_abundances.env)

   #macroinvert richness
macroinvert_richness_corr.l <- make_kendall_cor_plot(taxa = "Macroinvertebrate",metric_input = "",metric = "Species richness",metric_name = "total_richness_depthzone_site",data.table = dat_macroinvert_richness_depthzone_site.env)

   #macroinvert diversity density
macroinvert_simpson_abundensity_corr.l <- make_kendall_cor_plot(taxa = "Macroinvertebrate",metric_input = "Count density",metric = "Simpson diversity",metric_name = "Simpson_index_density",data.table = dat_fish_simpson.env)


#kelp
#kelp abundance density
kelp_abundance_abundensity_corr.l <- make_kendall_cor_plot(taxa = "Macroalgae",metric_input = "Count density",metric = "Abundance",metric_name = "total_abundance_depthzone_site",data.table = dat_kelp_total_abundances.env)

#kelp richness
kelp_richness_corr.l <- make_kendall_cor_plot(taxa = "Macroalgae",metric_input = "",metric = "Species richness",metric_name = "total_richness_depthzone_site",data.table = dat_kelp_richness_depthzone_site.env)

#kelp diversity density
kelp_simpson_abundensity_corr.l <- make_kendall_cor_plot(taxa = "Macroalgae",metric_input = "Count density",metric = "Simpson diversity",metric_name = "Simpson_index_density",data.table = dat_kelp_simpson.env)

#Merge long data tables
all_kendall_correlations_list <- mget(ls(pattern = "*_corr.l"))

#check to make sure this pulls in 11 individual data tables
stopifnot(length(all_kendall_correlations_list) == 11)

#merge this list of data tables into one data table
all_kendall_correlations <- rbindlist(all_kendall_correlations_list)

#set levels to improve order on full facetted plot
all_kendall_correlations[,metric_name := factor(metric_name,levels = c("total_richness_depthzone_site","total_abundance_depthzone_site", "total_biomass_depthzone_site","Simpson_index_density","Simpson_index_biomass"))]
all_kendall_correlations[,metric := factor(metric,levels = c("Species richness","Abundance","Simpson diversity"))]

#Make facet plot of correlations
(all_kendall_correlations_plot <- ggplot(all_kendall_correlations) +
  geom_tile(aes(y = DepthZone, x = predictor, fill = estimate_sig)) +
#  geom_text(aes(y = DepthZone, x = predictor, label = round(estimate_sig,2)), size = 2) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1,1)) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Habitat characteristic", y = "Depth zone", fill = "Kendall's tau-b") +
  geom_segment(aes(x=0.5,xend = 13.5, y = 1.5, yend = 1.5)) +
  facet_grid(rows = vars(taxa), cols = vars(metric_name), labeller = labeller(metric_name = as_labeller(metric_labels)), switch = "y") +
  theme_classic() +
  guides(fill = guide_colorbar(ticks.colour = "black", frame.colour = "black", frame.linewidth = 0.3)) +
  theme(axis.text.x = element_text(angle = 90, hjust = -0.005),
        axis.title.x = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        strip.placement = "outside",
  strip.background = element_blank(),
  strip.text = element_text(size = 13),
  axis.title.y = element_blank()
  ))

#save plot
ggsave(all_kendall_correlations_plot, path = file.path("figures","kendall_corr"), filename = "all_kendall_correlations_plot.jpg",
       width = 13, height = 8, units = "in")

#Influence of habitat on abundance and species richness NATURAL REEFS ONLY, and should exclude giant kelp cover for kelp analyses

#GLMMs are used for identifying species-habitat associates
#used to assess the influence of reef habitat factors on the abundance, richness, and diversity of assemblages
#include site as random effect which accounts for population sample bias within a site
#this assumes that most correlation occurs at the within-site scale, and not at the between-site scale (not necessarily true for us)
#Don't include co-linear variables in the same model
#evaluate for overdispersion using ratio of Pearson statistic and residual degrees of freedom to make sure parameter is <1.5, if not, use MCMCglmm
#full GLMMs run, but vary the inclusion of single variables from sets of collinear variables
#models ranked according to AIC

#I calculate relative variable importance for each taxa ~ parameter combo (table or plot)

#dredge
options(na.action = "na.fail")

#function that will output row of data table that will make matrix of coefficients and akaike weights
#output row of data table that will make correlation matrix
depthzones <- c("Inner","Middle","Outer","Deep")

rank_mod_output_coef_akaike <- function(dt, response_var, metric_category, taxa, metric){
  
  #empty data table for output
  output <- data.table()

  for (i in 1:length(depthzones)) {
    
    dt.bydepth <- dt[DepthZone == depthzones[i]]

#extract coefficient value from all model
for(j in 1:length(predictors)) {
  if(metric_category != "Simpson_diversity"){
  mod <- glm(data = dt.bydepth, get(response_var) ~ get(paste0(predictors[j],".s")), family = quasipoisson) #quasipoisson okay with overdispersed count data
  }else{
    mod <- glm(data = dt.bydepth, get(response_var) ~ get(paste0(predictors[j],".s")), family = binomial) #binomial 
  }
  #Extract r_squared
  r_squared <- r.squaredGLMM(mod)[[2,1]] #Delta method
  
  
  #Extract coefficient
  coefficient <- mod$coefficients[[2]]
  
  #Pvalue
  pvalue <- coef(summary(mod))[2,4]
  
  #Make row in data table
  row.dt <- data.table(metric_category = metric_category,
                       taxa = taxa,
                       metric = metric,
                       DepthZone = depthzones[i],
                       env_var = predictors[j],
                       coefficient = coefficient,
                       pvalue = pvalue,
                       r_squared = r_squared)
  
  output <- rbind(output,row.dt)
}
  }
  return(output)
}

#execute function on fish abundance
fish_abundance_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_fish_total_abundances.env,
                                                     response_var = "total_abundance_depthzone_site",
                                                     metric_category = "abundance",
                                                     metric = "density",
                                                     taxa = "fish")

#execute function on fish abundance (biomass)
fish_abundance_biomass_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_fish_total_abundances.env,
                                                     response_var = "total_biomass_depthzone_site",
                                                     metric_category = "abundance",
                                                     metric = "biomass",
                                                     taxa = "fish")

#execute function on kelp abundance
kelp_abundance_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_kelp_total_abundances.env,
                                                     response_var = "total_abundance_depthzone_site",
                                                     metric_category = "abundance",
                                                     metric = "density",
                                                     taxa = "kelp")
#execute function on macroinvert abundance
kelp_macroinvert_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_macroinvert_total_abundances.env,
                                                     response_var = "total_abundance_depthzone_site",
                                                     metric_category = "abundance",
                                                     metric = "density",
                                                     taxa = "macroinvert")

#execute function on fish richness
fish_richness_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_fish_richness_depthzone_site.env,
                                                                    response_var = "total_richness_depthzone_site",
                                                                    metric_category = "richness",
                                                                    metric = "density",
                                                                    taxa = "fish")
#execute function on kelp richness
kelp_richness_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_kelp_richness_depthzone_site.env,
                                                                   response_var = "total_richness_depthzone_site",
                                                                   metric_category = "richness",
                                                                   metric = "density",
                                                                   taxa = "kelp")
#execute function on macroinvert richness
macroinvert_richness_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_macroinvert_richness_depthzone_site.env,
                                                                   response_var = "total_richness_depthzone_site",
                                                                   metric_category = "richness",
                                                                   metric = "density",
                                                                   taxa = "macroinvert")

#execute function on fish Simpson diversity
fish_simpson_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_fish_simpson.env,
                                                                         response_var = "Simpson_index_density",
                                                                         metric_category = "Simpson_diversity",
                                                                         metric = "density",
                                                                         taxa = "fish")
#execute function on fish Simpson diversity (abundance)
fish_simpson_biomass_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_fish_simpson.env,
                                                                         response_var = "Simpson_index_biomass",
                                                                         metric_category = "Simpson_diversity",
                                                                         metric = "biomass",
                                                                         taxa = "fish")
#execute function on kelp Simpson diversity
kelp_simpson_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_kelp_simpson.env,
                                                                         response_var = "Simpson_index_density",
                                                                         metric_category = "Simpson_diversity",
                                                                         metric = "density",
                                                                         taxa = "kelp")
#execute function on macroinvert Simpson diversity
macroinvert_simpson_env_model_rank_output <- rank_mod_output_coef_akaike(dt = dat_macroinvert_simpson.env,
                                                                          response_var = "Simpson_index_density",
                                                                          metric_category = "Simpson_diversity",
                                                                          metric = "density",
                                                                          taxa = "macroinvert")

#Merge all
#First, list all
if(exists("all_taxa_env_model_rank_output.list")){
  rm(all_taxa_env_model_rank_output.list) # do avoid including mega list in the list of dt names
}

all_taxa_env_model_rank_output.list <- ls(pattern = "*env_model_rank_output")

#Then, merge data tables
if(exists("all_taxa_env_model_rank_output")){
  rm(all_taxa_env_model_rank_output) #if it exists, delete it
}
all_taxa_env_model_rank_output <- do.call(rbind, mget(all_taxa_env_model_rank_output.list))

#Make factor edits to each 

#Reorder depth zones
all_taxa_env_model_rank_output[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep"))]

#Rename taxa
all_taxa_env_model_rank_output[,taxa := factor(taxa, levels = c("fish","kelp","macroinvert"), labels = c("Fish","Macroalgae","Macroinvertebrates"))]

#Rename environmental variables and group by local vs. regional
all_taxa_env_model_rank_output[,env_var := factor(env_var, levels = c(
  "dist_200m_bath", "min_sst_C","mean_sst_C", "max_sst_C","min_chl_mg_m3" ,
  "mean_chl_mg_m3" ,"max_chl_mg_m3","giantkelp_stipe_density_m2", "Relief_SD","Relief_index", 
  "Substrate_SD","Substrate_index"),
  labels = c("Distance to 200 m isobath", "Minimum SST","Mean SST", "Maximum SST","Minimum chlorophyll" ,
             "Mean chlorophyll" ,"Maximum chlorophyll","Giant kelp stipe density", "Relief SD","Relief index","Substrate SD","Substrate index" 
             ))]

#Rename metrics
all_taxa_env_model_rank_output[,metric := factor(metric, levels = c(
  "density","biomass"),
  labels = c("Density (count)", "Density (biomass)"))]

#Rename metric categories
all_taxa_env_model_rank_output[,metric_category := factor(metric_category, levels = c(
  "abundance","richness","Simpson_diversity"),
  labels = c("Abundance","Richness","Simpson diversity index"))]


#visualize
all_taxa_env_model_rank_plot <- ggplot() +
  annotate(geom = 'rect', xmin=-Inf, xmax=Inf, ymin=8-0.5, ymax=Inf, alpha=0.4, fill = 'darkgrey')+ 
  geom_point(data = all_taxa_env_model_rank_output , aes(x = DepthZone, y = env_var, fill = r_squared, size = r_squared), shape = 21) +
  geom_text(data = all_taxa_env_model_rank_output[r_squared > 0.4 & pvalue <= 0.05], aes(x = DepthZone, y = env_var, label = round(coefficient,2)), size = 2) +
  scale_size_continuous(range = c(0,10)) +
  scale_fill_gradient(low = "white",high = "turquoise", guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  labs(x = "Depth zone", y = "Environmental variable", fill = "Delta R-squared", size = "") +
  facet_grid(metric_category + metric ~ taxa) +
    #guides(size = FALSE) +
  theme_classic() +
  scale_y_discrete() + #this allows the grey box to shade where I want it to!!!
  theme(legend.position = "top",
        legend.justification = "left",
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

#make rectangle
rect1 <- rectGrob(
  x = 0.62,
  y = 0.99,
  width = 0.06,
  height = 0.018,
  hjust = 0, vjust = 1,
  gp = gpar(fill = "darkgrey", alpha = 0.4)
)

rect2 <- rectGrob(
  x = 0.62,
  y = 0.973,
  width = 0.06,
  height = 0.018,
  hjust = 0, vjust = 1,
  gp = gpar(fill = "white", alpha = 0.4)
)

all_taxa_env_model_rank_plot <- ggdraw(all_taxa_env_model_rank_plot) +
  draw_grob(rect1) +
  draw_grob(rect2) +
  draw_text("Local", 
            x = 0.65,
            y = 0.983,
            size = 12) +
  draw_text("Regional", 
            x = 0.65,
              y = 0.963,
            size = 12) 
  

ggsave(all_taxa_env_model_rank_plot, path = "figures", filename = "all_taxa_env_model_rank_plot.jpg", height = 14, width = 14, units = "in", dpi = 500)


