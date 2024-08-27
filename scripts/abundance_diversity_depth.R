# CREATION DATE 7 July 2024
# MODIFIED DATE 26 August 2024

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

source(file.path("functions","add_tukey_letter_boxplot.R"))

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

#New Column Identifying ARM vs Island vs Natural Coast
dat_fish_site_averages[,type := factor(ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland")))]
dat_macroinvert_site_averages[,type := factor(ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland")))]
dat_kelp_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]


########################
##Total abundance (summed density of all taxa, summed biomass of all taxa) by depth zone
########################

#abundance for fish
dat_fish_site_averages[,total_abundance_depthzone_site := sum(mean_density_m2)*100,.(Site, DepthZone)] #multiply by 100 to get # per 100m^2
dat_fish_site_averages[,total_biomass_depthzone_site := sum(mean_wt_density_g_m2)*100/1000,.(Site, DepthZone)] #total fish biomass in kg per 100m^2 (multiply by 100, divide by 1000)

dat_fish_total_abundances <- unique(dat_fish_site_averages[,.(Site, DepthZone, type, total_abundance_depthzone_site, total_biomass_depthzone_site)])

#abundance for macroinverts
dat_macroinvert_site_averages[,total_abundance_depthzone_site := sum(mean_density_m2)*100,.(Site, DepthZone)] #multiply by 100 to get # per 100m^2
dat_macroinvert_total_abundances <- unique(dat_macroinvert_site_averages[,.(Site, DepthZone, type, total_abundance_depthzone_site)])


#abundance for kelp
dat_kelp_site_averages[,total_abundance_depthzone_site := sum(mean_density_m2)*100,.(Site, DepthZone)] #multiply by 100 to get # per 100m^2
dat_kelp_total_abundances <- unique(dat_kelp_site_averages[,.(Site, DepthZone, type, total_abundance_depthzone_site)])


#visualize all fish density
fish_abundance_depthzone <- ggplot(dat_fish_total_abundances) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_abundance_depthzone_site, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2,
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone",y = bquote("   Fish density\n(count per 100 m"^2*")") , fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.key.size = unit(2,"line")
  )

ggsave(fish_abundance_depthzone, path = file.path("figures"), filename ="fish_abundance_depthzone.jpg", height = 4.5, width = 5, units = "in")

(fish_abundance_depthzone_tukey <- generate_boxplot_tukey_labels(dat_fish_total_abundances, "total_abundance_depthzone_site", fish_abundance_depthzone))

#visualize all fish biomass

fish_biomass_depthzone <- ggplot(dat_fish_total_abundances) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_biomass_depthzone_site, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = bquote("   Fish density\n(kg per 100 m"^2*")"), fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )

ggsave(fish_biomass_depthzone, path = file.path("figures"), filename ="fish_biomass_depthzone.jpg", height = 4.5, width = 5, units = "in")

(fish_biomass_depthzone_tukey <- generate_boxplot_tukey_labels(dat_fish_total_abundances, "total_biomass_depthzone_site", fish_biomass_depthzone))

#visualize all macroinvertebrate density
macroinvert_abundance_depthzone <- ggplot(dat_macroinvert_total_abundances) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_abundance_depthzone_site, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = bquote("Macroinvertebrate density (count per 100 m"^2*")"), fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.85, 1.1),
    legend.justification = c("right", "top")
  )

ggsave(macroinvert_abundance_depthzone, path = file.path("figures"), filename ="macroinvert_abundance_depthzone.jpg", height = 4.5, width = 5, units = "in")

macroinvert_abundance_depthzone_tukey <- generate_boxplot_tukey_labels(dat_macroinvert_total_abundances, "total_abundance_depthzone_site", macroinvert_abundance_depthzone)

#visualize all kelp density
kelp_abundance_depthzone <- ggplot(dat_kelp_total_abundances) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_abundance_depthzone_site, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = bquote("Macroalgae density (count per 100 m"^2*")"), fill = "Reef type", pattern = "Reef type") +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","AR")) +
  theme_classic() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )

ggsave(kelp_abundance_depthzone, path = file.path("figures"), filename ="kelp_abundance_depthzone.jpg", height = 4.5, width = 5, units = "in")

kelp_abundance_depthzone_tukey <- generate_boxplot_tukey_labels(dat_kelp_total_abundances, "total_abundance_depthzone_site", kelp_abundance_depthzone)

#Merge plots
#first, stack fish biomass and density
fish_abundance_merge <- cowplot::plot_grid(fish_abundance_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank(), plot.margin = margin(0,0,0,1,unit = "cm")), 
                                           fish_biomass_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank(), plot.margin = margin(0,0,0,1,unit = "cm")),
                                                                          ncol = 1, align = "v", labels = c("a.","b."), label_y = c(1,1.1))

abundance_depthzone_merge <- cowplot::plot_grid(fish_abundance_merge,
                                                kelp_abundance_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
                                                macroinvert_abundance_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
                                                
                                                                                                nrow = 1, labels = c("","c.", "d."), align = "h", rel_widths = c(1.1,1,1))

########################
##Total richness (summed # of taxa, summed biomass of all taxa) by depth zone
########################
#total # across all sites
#fish
#only rows with values
dat_fish_site_averages.r <- dat_fish_site_averages[mean_density_m2 >0,]
dat_fish_site_averages.r[,total_richness_depthzone_site := uniqueN(taxa), by= .(Site, DepthZone, type)]
dat_fish_richness_depthzone_site <- unique(dat_fish_site_averages.r[,.(Site, DepthZone, type, total_richness_depthzone_site)])
#macroinvert
dat_macroinvert_site_averages.r <- dat_macroinvert_site_averages[mean_density_m2 >0,]
dat_macroinvert_site_averages.r[,total_richness_depthzone_site := uniqueN(taxa), by= .(Site, DepthZone, type)]
dat_macroinvert_richness_depthzone_site <- unique(dat_macroinvert_site_averages.r[,.(Site, DepthZone, type, total_richness_depthzone_site)])
#kelp
dat_kelp_site_averages.r <- dat_kelp_site_averages[mean_density_m2 >0,]
dat_kelp_site_averages.r[,total_richness_depthzone_site := uniqueN(taxa), by= .(Site, DepthZone, type)]
dat_kelp_richness_depthzone_site <- unique(dat_kelp_site_averages.r[,.(Site, DepthZone, type, total_richness_depthzone_site)])


#visualize all fish richness
fish_richness_depthzone <- ggplot(dat_fish_richness_depthzone_site) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_richness_depthzone_site, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = "Fish richness\n", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.4, 1),
    legend.justification = c("right", "top")
  )

ggsave(fish_richness_depthzone, path = file.path("figures"), filename ="fish_richness_depthzone.jpg", height = 4.5, width = 5, units = "in")

fish_richness_depthzone_tukey <- generate_boxplot_tukey_labels(dat_fish_richness_depthzone_site, "total_richness_depthzone_site", fish_richness_depthzone, richness = T)

#visualize all macroinvertebrate richness
macroinvert_richness_depthzone <- ggplot(dat_macroinvert_richness_depthzone_site) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_richness_depthzone_site, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = "Macroinvertebrate richness", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.4, 1),
    legend.justification = c("right", "top")
  )

ggsave(macroinvert_richness_depthzone, path = file.path("figures"), filename ="macroinvert_richness_depthzone.jpg", height = 4.5, width = 5, units = "in")

macroinvert_richness_depthzone_tukey <- generate_boxplot_tukey_labels(dat_macroinvert_richness_depthzone_site, "total_richness_depthzone_site", macroinvert_richness_depthzone, richness = T)

#visualize all kelp richness
kelp_richness_depthzone <- ggplot(dat_kelp_richness_depthzone_site) +
  geom_boxplot_pattern(aes(x = DepthZone, y = total_richness_depthzone_site, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = "Macroalgae richness", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.95, 1),
    legend.justification = c("right", "top")
  )

ggsave(kelp_richness_depthzone, path = file.path("figures"), filename ="kelp_richness_depthzone.jpg", height = 4.5, width = 5, units = "in")


kelp_richness_depthzone_tukey <- generate_boxplot_tukey_labels(dat_kelp_richness_depthzone_site, "total_richness_depthzone_site", kelp_richness_depthzone, richness = T)


#Merge plots
richness_depthzone_merge <- cowplot::plot_grid(fish_richness_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank(), plot.margin = margin(0,0,0,1,unit = "cm")),
                                                kelp_richness_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank(), plot.margin = margin(0,0,0,0.8,unit = "cm")),
                                               macroinvert_richness_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank(), plot.margin = margin(0,0,0,0.8,unit = "cm")),
                                               
                                                                                               nrow = 1, labels = c("e.","f.", "g."), align = "h")

########################
##Simpson diversity index
########################

#first, take proportion of each species
#then, square the proportions and sum these for each site
#fish
dat_fish_site_averages.r[,total_site_depthzone_density := sum(mean_density_m2),.(Site, DepthZone)]
dat_fish_site_averages.r[,total_site_depthzone_biomass := sum(mean_wt_density_g_m2),.(Site, DepthZone)]
dat_fish_site_averages.r[,taxa_proportion_density := mean_density_m2/total_site_depthzone_density]
dat_fish_site_averages.r[,taxa_proportion_biomass := mean_wt_density_g_m2/total_site_depthzone_biomass]
dat_fish_site_averages.r[,Simpson_index_density := 1-sum(taxa_proportion_density^2),.(Site,DepthZone)]
dat_fish_site_averages.r[,Simpson_index_biomass := 1-sum(taxa_proportion_biomass^2),.(Site,DepthZone)]

dat_fish_simpson <- unique(dat_fish_site_averages.r[,.(Site, DepthZone, type, Simpson_index_density, Simpson_index_biomass)])

#macroinvert
dat_macroinvert_site_averages.r[,total_site_depthzone_density := sum(mean_density_m2),.(Site, DepthZone)]
dat_macroinvert_site_averages.r[,taxa_proportion_density := mean_density_m2/total_site_depthzone_density]
dat_macroinvert_site_averages.r[,Simpson_index_density := 1-sum(taxa_proportion_density^2),.(Site,DepthZone)]

dat_macroinvert_simpson <- unique(dat_macroinvert_site_averages.r[,.(Site, DepthZone, type, Simpson_index_density)])

#macroalgae
dat_kelp_site_averages.r[,total_site_depthzone_density := sum(mean_density_m2),.(Site, DepthZone)]
dat_kelp_site_averages.r[,taxa_proportion_density := mean_density_m2/total_site_depthzone_density]
dat_kelp_site_averages.r[,Simpson_index_density := 1-sum(taxa_proportion_density^2),.(Site,DepthZone)]

dat_kelp_simpson <- unique(dat_kelp_site_averages.r[,.(Site, DepthZone, type, Simpson_index_density)])


###VISUALIZE

#visualize fish diversity (calculated with biomass)
fish_simpson_depthzone_density <- ggplot(dat_fish_simpson) +
  geom_boxplot_pattern(aes(x = DepthZone, y = Simpson_index_density, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = "Fish Simpson index\n(count-based)", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.87, 0.15)
  )

ggsave(fish_simpson_depthzone_density, path = file.path("figures"), filename ="fish_simpson_depthzone_density.jpg", height = 4.5, width = 6, units = "in")

fish_simpson_depthzone_density_tukey <- generate_boxplot_tukey_labels(dat_fish_simpson, "Simpson_index_density", fish_simpson_depthzone_density)

#same but for fish biomass

fish_simpson_depthzone_biomass <- ggplot(dat_fish_simpson) +
  geom_boxplot_pattern(aes(x = DepthZone, y = Simpson_index_biomass, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = "Fish Simpson index\n(biomass-based)", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.87, 0.15)
  )

ggsave(fish_simpson_depthzone_biomass, path = file.path("figures"), filename ="fish_simpson_depthzone_biomass.jpg", height = 4.5, width = 6, units = "in")

summary(lm(data = dat_fish_simpson[type != "ARM"], Simpson_index_biomass~DepthZone))

#When using biomass instead of abundance to calculate diversity, diversity is much more consistent across Depth Zones
#and exhibit little difference between island and mainland reefs, largely due to the small size of high abundance species such as blacksmith and gobies

summary(lm(data = dat_fish_simpson[type != "ARM"], Simpson_index_biomass~DepthZone*type))

fish_simpson_depthzone_biomass_tukey <- generate_boxplot_tukey_labels(dat_fish_simpson, "Simpson_index_biomass", fish_simpson_depthzone_biomass)


#visualize macroinvert diversity (calculated with density)
macroinvert_simpson_depthzone_density <- ggplot(dat_macroinvert_simpson) +
  geom_boxplot_pattern(aes(x = DepthZone, y = Simpson_index_density, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = "Macroinvertebrate Simpson index", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.87, 0.15)
  )

macroinvert_simpson_depthzone_tukey <- generate_boxplot_tukey_labels(dat_macroinvert_simpson, "Simpson_index_density", macroinvert_simpson_depthzone_density)


#visualize kelp diversity (calculated with density)
kelp_simpson_depthzone_density <- ggplot(dat_kelp_simpson) +
  geom_boxplot_pattern(aes(x = DepthZone, y = Simpson_index_density, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = "Macroalgae Simpson index", fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.87, 0.85)
  )

ggsave(kelp_simpson_depthzone_density, path = file.path("figures"), filename ="kelp_simpson_depthzone_density.jpg", height = 4.5, width = 6, units = "in")

kelp_simpson_depthzone_density_tukey <- generate_boxplot_tukey_labels(dat_kelp_simpson, "Simpson_index_density", kelp_simpson_depthzone_density)

#Merge plots
#first, stack fish biomass and density
fish_diversity_merge <- cowplot::plot_grid(fish_simpson_depthzone_density_tukey  + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank(),  plot.margin = margin(0,0,0,0.6,unit = "cm"))  , 
                                           fish_simpson_depthzone_biomass_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank(),  plot.margin = margin(0,0,0,0.6,unit = "cm")),
                                           ncol = 1, align = "v", labels = c("h.","i."), label_y = c(1,1.1))

diversity_depthzone_merge <- cowplot::plot_grid(fish_diversity_merge,
                                                kelp_simpson_depthzone_density_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
                                                macroinvert_simpson_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
                                                
                                                                                                nrow = 1, labels = c("","j.", "k."), align = "h", rel_widths = c(1.1,1,1))

#mega-merge
depthzone_abun_rich_div_merge_tukey <- cowplot::plot_grid(abundance_depthzone_merge, richness_depthzone_merge, diversity_depthzone_merge, ncol = 1, align = "hv")

#add legend on top
megamerge_tukey_leg <- get_legend(fish_simpson_depthzone_density_tukey + theme(legend.position = "top", legend.direction = "horizontal"))

depthzone_abun_rich_div_merge_tukey.l <- plot_grid(megamerge_tukey_leg, depthzone_abun_rich_div_merge_tukey, rel_heights = c(1,10), ncol = 1)

#save
ggsave(depthzone_abun_rich_div_merge_tukey.l, path = file.path("figures"), filename = "depthzone_abun_rich_div_merge_tukey.l.jpg", height = 12, width = 10, unit = "in")





#as an alternate route for plotting (facet), try to merge all data tables and then flip from wide to long

#link together
#abundance/density
dat_fish_total_abundances[, taxa_type := "Fish"]
dat_fish_total_abundances.l <- melt(dat_fish_total_abundances, id.vars = c("taxa_type","Site","DepthZone","type"),measure.vars = c("total_abundance_depthzone_site","total_biomass_depthzone_site"), variable.name = "metric",value.name = "value")
dat_macroinvert_total_abundances[, taxa_type := "Macroinvertebrate"]
dat_macroinvert_total_abundances.l <- melt(dat_macroinvert_total_abundances, id.vars = c("taxa_type","Site","DepthZone","type"),measure.vars = c("total_abundance_depthzone_site"), variable.name = "metric",value.name = "value")
dat_kelp_total_abundances[, taxa_type := "Macroalgae"]
dat_kelp_total_abundances.l <- melt(dat_kelp_total_abundances, id.vars = c("taxa_type","Site","DepthZone","type"), variable.name = "metric",measure.vars = c("total_abundance_depthzone_site"), value.name = "value")


#richness
dat_fish_richness_depthzone_site[,taxa_type := "Fish"]
dat_fish_richness_depthzone_site.l <- melt(dat_fish_richness_depthzone_site, id.vars = c("taxa_type","Site","DepthZone","type"),measure.vars = c("total_richness_depthzone_site"), variable.name = "metric",value.name = "value")
dat_macroinvert_richness_depthzone_site[,taxa_type := "Macroinvertebrate"]
dat_macroinvert_richness_depthzone_site.l <- melt(dat_macroinvert_richness_depthzone_site, id.vars = c("taxa_type","Site","DepthZone","type"),measure.vars = c("total_richness_depthzone_site"), variable.name = "metric",value.name = "value")
dat_kelp_richness_depthzone_site[,taxa_type := "Macroalgae"]
dat_kelp_richness_depthzone_site.l <- melt(dat_kelp_richness_depthzone_site, id.vars = c("taxa_type","Site","DepthZone","type"),measure.vars = c("total_richness_depthzone_site"), variable.name = "metric",value.name = "value")

#simpson diversity index
dat_fish_simpson[,taxa_type := "Fish"]
dat_fish_simpson.l <- melt(dat_fish_simpson, id.vars = c("taxa_type","Site","DepthZone","type"),measure.vars = c("Simpson_index_density","Simpson_index_biomass"), variable.name = "metric",value.name = "value")
dat_macroinvert_simpson[,taxa_type := "Macroinvertebrate"]
dat_macroinvert_simpson.l <- melt(dat_macroinvert_simpson, id.vars = c("taxa_type","Site","DepthZone","type"),measure.vars = c("Simpson_index_density"), variable.name = "metric",value.name = "value")
dat_kelp_simpson[,taxa_type := "Macroalgae"]
dat_kelp_simpson.l <- melt(dat_kelp_simpson, id.vars = c("taxa_type","Site","DepthZone","type"),measure.vars = c("Simpson_index_density"), variable.name = "metric",value.name = "value")

#rbind all above

merged_abundance_richness_diversity_datatables <- rbind(dat_fish_total_abundances.l, dat_macroinvert_total_abundances.l, dat_kelp_total_abundances.l,
                                                        dat_fish_richness_depthzone_site.l, dat_macroinvert_richness_depthzone_site.l, dat_kelp_richness_depthzone_site.l,
                                                        dat_fish_simpson.l, dat_macroinvert_simpson.l, dat_kelp_simpson.l)

# Custom labels for the 'metric' variable
metric_labels <- c(
  total_abundance_depthzone_site = "Density\n(count-based)",
  total_biomass_depthzone_site = "Density\n(biomass-based)",
  total_richness_depthzone_site = "Richness\n ",
  Simpson_index_density = "Simpson index\n(count-based)",
  Simpson_index_biomass = "Simpson index\n(biomass-based)"
)

#change ARM to AR only
merged_abundance_richness_diversity_datatables[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","ARM"), labels = c("Inner","Middle","Outer","Deep","AR"))]

# Convert 'metric' column to factor
merged_abundance_richness_diversity_datatables[, metric := factor(metric, levels = names(metric_labels))]


#merged plot
#visualize kelp diversity (calculated with density)
depthzone_density_abundance_richness <- ggplot(merged_abundance_richness_diversity_datatables) +
  geom_boxplot_pattern(aes(x = DepthZone, y = value, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = "Value", fill = "Reef type", pattern = "Reef type") +
  facet_grid(rows = vars(metric),cols = vars(taxa_type), scales = "free", labeller = labeller(metric = as_labeller(metric_labels)), switch = "y") +
  theme_classic() +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  geom_vline(aes(xintercept = 4.65), color = "black", linewidth = 0.8) +
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

#################################################################################
#How do the predictors of abundance, richness, and diversity vary by depth zone
#################################################################################

#pull environmental data
all_env_lat_lon <- fread(file.path("data","enviro_predictors","all_env_lat_lon.csv"))

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


