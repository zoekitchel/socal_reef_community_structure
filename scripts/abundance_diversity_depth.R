# CREATION DATE 7 July 2024
# MODIFIED DATE 11 July 2024

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
dat_fish_site_averages[,total_abundance_depthzone_site := sum(mean_density_m2),.(Site, DepthZone)]
dat_fish_site_averages[,total_biomass_depthzone_site := sum(mean_wt_density_g_m2),.(Site, DepthZone)] #total fish biomass in grams

dat_fish_total_abundances <- unique(dat_fish_site_averages[,.(Site, DepthZone, type, total_abundance_depthzone_site, total_biomass_depthzone_site)])

#abundance for macroinverts
dat_macroinvert_site_averages[,total_abundance_depthzone_site := sum(mean_density_m2),.(Site, DepthZone)]
dat_macroinvert_total_abundances <- unique(dat_macroinvert_site_averages[,.(Site, DepthZone, type, total_abundance_depthzone_site)])


#abundance for kelp
dat_kelp_site_averages[,total_abundance_depthzone_site := sum(mean_density_m2),.(Site, DepthZone)]
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
  labs(x = "Depth Zone",y = bquote("Fish density (#/m"^2*")") , fill = "Reef type", pattern = "Reef type") +
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
  geom_boxplot_pattern(aes(x = DepthZone, y = total_biomass_depthzone_site/1000, fill = type, pattern = type), position = position_dodge2(preserve = "single"),
                       colour          = 'black', 
                       pattern_density = 0.01, 
                       pattern_fill    = 'black',
                       pattern_colour  = 'black',
                       pattern_spacing = 0.05, outlier.alpha = 0.3, outlier.size = 1.2
  ) +
  scale_fill_manual(values = c("#FB8071","#7FB1D3","#BDBAD9"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  scale_pattern_manual(values = c("none","stripe","none"), labels = c("Artificial reef","Natural island","Natural mainland")) +
  labs(x = "Depth Zone", y = bquote("Fish biomass (kg/m"^2*")"), fill = "Reef type", pattern = "Reef type") +
  theme_classic() +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top")
  )

ggsave(fish_biomass_depthzone, path = file.path("figures"), filename ="fish_biomass_depthzone.jpg", height = 4.5, width = 5, units = "in")

(fish_biomass_depthzone_tukey <- generate_boxplot_tukey_labels(dat_fish_total_abundances, "total_biomass_depthzone_site", fish_biomass_depthzone, divide_var_1000 = T))

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
  labs(x = "Depth Zone", y = bquote("Macroinvertebrate density (#/m"^2*")"), fill = "Reef type", pattern = "Reef type") +
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
  labs(x = "Depth Zone", y = bquote("Macroalgae density (#/m"^2*")"), fill = "Reef type", pattern = "Reef type") +
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
fish_abundance_merge <- cowplot::plot_grid(fish_abundance_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()), 
                                           fish_biomass_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
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
  labs(x = "Depth Zone", y = "Fish richness", fill = "Reef type", pattern = "Reef type") +
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
richness_depthzone_merge <- cowplot::plot_grid(fish_richness_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
                                                kelp_richness_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
                                               macroinvert_richness_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
                                               
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
  labs(x = "Depth Zone", y = "Density", fill = "Reef type", pattern = "Reef type") +
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
  labs(x = "Depth Zone", y = "Biomass", fill = "Reef type", pattern = "Reef type") +
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
fish_diversity_merge <- cowplot::plot_grid(fish_simpson_depthzone_density_tukey  + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()), 
                                           fish_simpson_depthzone_biomass_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
                                           ncol = 1, align = "v", labels = c("h.","i."), label_y = c(1,1.1))

diversity_depthzone_merge <- cowplot::plot_grid(fish_diversity_merge,
                                                kelp_simpson_depthzone_density_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
                                                macroinvert_simpson_depthzone_tukey + theme(legend.position = "null", axis.title.x = element_blank(), axis.text = element_text(size = 12), axis.text.x = element_blank()),
                                                
                                                                                                nrow = 1, labels = c("","j.", "k."), align = "h", rel_widths = c(1.1,1,1))

#mega-merge
depthzone_abun_rich_div_merge_tukey <- cowplot::plot_grid(abundance_depthzone_merge, richness_depthzone_merge, diversity_depthzone_merge, ncol = 1, align = "hv")

#save
ggsave(depthzone_abun_rich_div_merge_tukey, path = file.path("figures"), filename = "depthzone_abun_rich_div_merge_tukey.jpg", height = 9, width = 10, unit = "in")





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
  total_abundance_depthzone_site = "Density (#/m^2)",
  total_biomass_depthzone_site = "Density (kg/m^2)",
  total_richness_depthzone_site = "Richness",
  Simpson_index_density = "Simpson index (#)",
  Simpson_index_biomass = "Simpson index (kg)"
)

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
  theme(
    legend.position = c(0.67, 0.7),
    panel.border = element_rect(color = "black", fill = NA),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    axis.title.y = element_blank(),
    legend.key.size = unit(2.5,"lines"),
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 13),
    axis.text.y = element_text(size = 10)
  )

ggsave(depthzone_density_abundance_richness, path = "figures",filename = "depthzone_density_abundance_richness.jpg", height = 10, width = 8.5, units = "in")

#################################################################################
#How do the predictors of abundance, richness, and diversity vary by depth zone
#################################################################################

#pull environmental data
all_env_lat_lon <- fread(file.path("data","enviro_predictors","all_env_lat_lon.csv"))

#NEED TO REMOVE ARM SITES FOR THESE

#link
#abundance/density
dat_fish_total_abundances.env <- all_env_lat_lon[dat_fish_total_abundances, on = c("Site","DepthZone")]
dat_fish_total_abundances.env <- dat_fish_total_abundances.env[complete.cases(dat_fish_total_abundances.env)]
dat_macroinvert_total_abundances.env <- all_env_lat_lon[dat_macroinvert_total_abundances, on = c("Site","DepthZone")]
dat_macroinvert_total_abundances.env <- dat_macroinvert_total_abundances.env[complete.cases(dat_macroinvert_total_abundances.env)]
dat_kelp_total_abundances.env <- all_env_lat_lon[dat_kelp_total_abundances, on = c("Site","DepthZone")]
dat_kelp_total_abundances.env <- dat_kelp_total_abundances.env[complete.cases(dat_kelp_total_abundances.env)]
#remove ARM from these analyses
#dat_fish_total_abundances.env <- dat_fish_total_abundances.env[type != "ARM"]
#dat_macroinvert_total_abundances.env <- dat_macroinvert_total_abundances.env[type != "ARM"]
#dat_kelp_total_abundances.env <- dat_kelp_total_abundances.env[type != "ARM"]

#richness
dat_fish_richness_depthzone_site.env <- all_env_lat_lon[dat_fish_richness_depthzone_site, on = c("Site","DepthZone")]
dat_fish_richness_depthzone_site.env <- dat_fish_richness_depthzone_site.env[complete.cases(dat_fish_richness_depthzone_site.env)]
dat_macroinvert_richness_depthzone_site.env <- all_env_lat_lon[dat_macroinvert_richness_depthzone_site, on = c("Site","DepthZone")]
dat_macroinvert_richness_depthzone_site.env <- dat_macroinvert_richness_depthzone_site.env[complete.cases(dat_macroinvert_richness_depthzone_site.env)]
dat_kelp_richness_depthzone_site.env <- all_env_lat_lon[dat_kelp_richness_depthzone_site, on = c("Site","DepthZone")]
dat_kelp_richness_depthzone_site.env <- dat_kelp_richness_depthzone_site.env[complete.cases(dat_kelp_richness_depthzone_site.env)]
#remove ARM from these analyses
#dat_fish_richness_depthzone_site.env <- dat_fish_richness_depthzone_site.env[type != (#)
#dat_macroinvert_richness_depthzone_site.env <- dat_macroinvert_richness_depthzone_site.env[type != (#)
#dat_kelp_richness_depthzone_site.env <- dat_kelp_richness_depthzone_site.env[type != (#)

#simpson diversity index
dat_fish_simpson.env <- all_env_lat_lon[dat_fish_simpson, on = c("Site","DepthZone")]
dat_fish_simpson.env <- dat_fish_simpson.env[complete.cases(dat_fish_simpson.env)]
dat_macroinvert_simpson.env <- all_env_lat_lon[dat_macroinvert_simpson, on = c("Site","DepthZone")]
dat_macroinvert_simpson.env <- dat_macroinvert_simpson.env[complete.cases(dat_macroinvert_simpson.env)]
dat_kelp_simpson.env <- all_env_lat_lon[dat_kelp_simpson, on = c("Site","DepthZone")]
dat_kelp_simpson.env <- dat_kelp_simpson.env[complete.cases(dat_kelp_simpson.env)]
#remove ARM from these analyses
#dat_fish_simpson.env <- dat_fish_simpson.env[type != "ARM"]
#dat_macroinvert_simpson.env <- dat_macroinvert_simpson.env[type != "ARM"]
#dat_kelp_simpson.env <- dat_kelp_simpson.env[type != "ARM"]

library(ggcorrplot)

data.table <- dat_fish_total_abundances.env
taxa <- "Fish"
metric_input <- "density"
metric <- "abundance"
metric_name <- "total_abundance_depthzone_site"

#Spearman correlation plots by variable
make_spearman_cor_plot <- function(taxa, metric_input, metric, metric_name, data.table){
  
  data.table <- data.table(data.table)
  
  metric_data <- data.table[[metric_name]]
  
  data.table[,metric_named := metric_data]
  
  #deep
  data.table.deep <- data.table[DepthZone == "Deep",.(metric_named,Latitude, dist_200m_bath, min_sst_C, mean_sst_C, max_sst_C, min_chl_mg_m3, mean_chl_mg_m3, max_chl_mg_m3, giantkelp_stipe_density_m2, Relief_index, Relief_SD, Substrate_index, Substrate_SD)]
  
  cor_d <- cor(y = data.table.deep$metric_named, x = data.table.deep[,2:14], use = "everything", method = "spearman")
  colnames(cor_d) <- "Deep"
  
  
  #outer
  data.table.outer <- data.table[DepthZone == "Outer",.(metric_named,Latitude, dist_200m_bath, min_sst_C, mean_sst_C, max_sst_C, min_chl_mg_m3, mean_chl_mg_m3, max_chl_mg_m3, giantkelp_stipe_density_m2, Relief_index, Relief_SD, Substrate_index, Substrate_SD)]
  
  
  cor_o <- cor(y = data.table.outer$metric_named, x = data.table.outer[,2:14], use = "everything", method = "spearman")
  colnames(cor_o) <- "Outer"
  
  
  #middle
  data.table.middle <- data.table[DepthZone == "Middle",.(metric_named,Latitude, dist_200m_bath, min_sst_C, mean_sst_C, max_sst_C, min_chl_mg_m3, mean_chl_mg_m3, max_chl_mg_m3, giantkelp_stipe_density_m2, Relief_index, Relief_SD, Substrate_index, Substrate_SD)]
  
  
  cor_m <- cor(y = data.table.middle$metric_named, x = data.table.middle[,2:14], use = "everything", method = "spearman")
  colnames(cor_m) <- "Middle"
  
  
  #inner
  data.table.inner <- data.table[DepthZone == "Inner",.(metric_named,Latitude, dist_200m_bath, min_sst_C, mean_sst_C, max_sst_C, min_chl_mg_m3, mean_chl_mg_m3, max_chl_mg_m3, giantkelp_stipe_density_m2, Relief_index, Relief_SD, Substrate_index, Substrate_SD)]
  
  
  cor_i <- cor(y = data.table.inner$metric_named, x = data.table.inner[,2:14], use = "everything", method = "spearman")
  colnames(cor_i) <- "Inner"
  
  #merge
  cor_merge <- cbind(cor_d, cor_o,cor_m, cor_i)
  cor_merge.dt <- data.table(cor_merge)
  cor_merge.dt[,env_variable := rownames(cor_merge)]
  
  #Plot
  
  #flip wide to long
  #add identifying columns
  cor_merge.dt[,taxa  := taxa][,metric_input := metric_input][,metric := metric]
  cor_merge.l <- melt(cor_merge.dt, id.vars = c("taxa","metric_input","metric", "env_variable"), variable.name = "DepthZone",value.name = "coefficient")
  
  cor_merge.l[,DepthZone := factor(DepthZone, levels = c("Deep","Outer","Middle","Inner"))]
  
  #add better labels and order for environmental variables
cor_merge.l[,env_variable := factor(env_variable, levels = c(
                                                 "Latitude", "dist_200m_bath", "min_sst_C","mean_sst_C", "max_sst_C","min_chl_mg_m3" ,
                                                 "mean_chl_mg_m3" ,"max_chl_mg_m3","giantkelp_stipe_density_m2", "Relief_index", "Relief_SD","Substrate_index" ,
                                                 "Substrate_SD"),
                                    labels = c("Latitude", "Distance to 200 m isobath", "Minimum SST","Mean SST", "Maximum SST","Minimum chlorophyll" ,
                                               "Mean chlorophyll" ,"Maximum chlorophyll","Giant kelp stipe density", "Relief index", "Relief SD","Substrate index" ,
                                               "Substrate SD"))]

plot_title <- paste(taxa, metric_input, metric, sep = " ")  

  single_plot <- ggplot(cor_merge.l) +
    geom_tile(aes(y = DepthZone, x = env_variable, fill = coefficient)) +
    geom_text(aes(y = DepthZone, x = env_variable, label = round(coefficient,2)), size = 2) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1,1)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Habitat characteristic", y = "Depth zone", fill = "Spearman's\ncorrelation\ncoefficient", title = plot_title) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = -0.1), axis.title.x = element_blank(), plot.title = element_text(face = "bold"))

plot_name <- paste0(taxa,"_",metric_input,"_", metric,"_spearmans_corr_heatmap")

ggsave(path = file.path("figures","spearman_corr"), filename = paste0(plot_name,".jpg"), height = 5, width = 5, unit = "in")

return(cor_merge.l)
}

#fish
   #fish abundance density
fish_abundance_abundensity_corr.l <- make_spearman_cor_plot(taxa = "Fish",metric_input = "count density",metric = "abundance",metric_name = "total_abundance_depthzone_site",data.table = dat_fish_total_abundances.env)
   #fish abundance biomass
fish_abundance_biodensity_corr.l <- make_spearman_cor_plot(taxa = "Fish",metric_input = "biomass density",metric = "abundance",metric_name = "total_biomass_depthzone_site",data.table = dat_fish_total_abundances.env)
   #fish richness
fish_richness_corr.l <- make_spearman_cor_plot(taxa = "Fish",metric_input = "",metric = "richness",metric_name = "total_richness_depthzone_site",data.table = dat_fish_richness_depthzone_site.env)
   #fish diversity density
fish_simpson_abundensity_corr.l <- make_spearman_cor_plot(taxa = "Fish",metric_input = "count density",metric = "Simpson diversity",metric_name = "Simpson_index_density",data.table = dat_fish_simpson.env)
   #fish diversity biomass
fish_simpson_biodensity_corr.l <- make_spearman_cor_plot(taxa = "Fish",metric_input = "biomass density",metric = "Simpson diversity",metric_name = "Simpson_index_biomass",data.table = dat_fish_simpson.env)

#macroinvertebrates
   #macroinvert abundance density
macroinvert_abundance_abundensity_corr.l <- make_spearman_cor_plot(taxa = "Macroinvertebrate",metric_input = "",metric = "abundance",metric_name = "total_abundance_depthzone_site",data.table = dat_macroinvert_total_abundances.env)

   #macroinvert richness
macroinvert_richness_corr.l <- make_spearman_cor_plot(taxa = "Macroinvertebrate",metric_input = "",metric = "richness",metric_name = "total_richness_depthzone_site",data.table = dat_macroinvert_richness_depthzone_site.env)

   #macroinvert diversity density
macroinvert_simpson_abundensity_corr.l <- make_spearman_cor_plot(taxa = "Macroinvertebrate",metric_input = "",metric = "Simpson diversity",metric_name = "Simpson_index_density",data.table = dat_fish_simpson.env)


#kelp
#kelp abundance density
kelp_abundance_abundensity_corr.l <- make_spearman_cor_plot(taxa = "Macroalgae",metric_input = "",metric = "abundance",metric_name = "total_abundance_depthzone_site",data.table = dat_kelp_total_abundances.env)

#kelp richness
kelp_richness_corr.l <- make_spearman_cor_plot(taxa = "Macroalgae",metric_input = "",metric = "richness",metric_name = "total_richness_depthzone_site",data.table = dat_kelp_richness_depthzone_site.env)

#kelp diversity density
kelp_simpson_abundensity_corr.l <- make_spearman_cor_plot(taxa = "Macroalgae",metric_input = "",metric = "Simpson diversity",metric_name = "Simpson_index_density",data.table = dat_kelp_simpson.env)

#Merge long data tables
all_spearman_correlations_list <- mget(ls(pattern = ".*_"))
#merge this list into one
#merge these individual objects
all_spearman_correlations <- rbindlist(all_spearman_correlations_list)

#Make facet plot
all_spearman_correlations_plot <- ggplot(all_spearman_correlations) +
  geom_tile(aes(y = DepthZone, x = env_variable, fill = coefficient)) +
#  geom_text(aes(y = DepthZone, x = env_variable, label = round(coefficient,2)), size = 2) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0, limits = c(-1,1)) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Habitat characteristic", y = "Depth zone", fill = "Spearman's\ncorrelation\ncoefficient", title = plot_title) +
  facet_wrap(~ taxa + metric + metric_input) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = -0.1), axis.title.x = element_blank(), plot.title = element_text(face = "bold"))


#MODEL FISH ABUNDANCE/DENSITY BY DEPTH ZONE

#scale all variables
# Identify numeric columns
numeric_cols <- names(dat_fish_total_abundances.env)[sapply(dat_fish_total_abundances.env, is.numeric) & !(names(dat_fish_total_abundances.env) %in% c("ID", "true_label"))]

# Scale numeric columns and create new columns with .s appended
dat_fish_total_abundances.env[, (paste0(numeric_cols, ".s")) := lapply(.SD, scale), .SDcols = numeric_cols]

#exclude two sites for which we don't have data for chlorophyll or sst
dat_fish_total_abundances.env <- dat_fish_total_abundances.env[complete.cases(dat_fish_total_abundances.env),]

options(na.action = "na.fail")

#DEEP
dat_fish_total_abundances.env.deep <- dat_fish_total_abundances.env[DepthZone == "Deep",]

#dredge
global.fish.model_deep <- lm(data = dat_fish_total_abundances.env.deep, total_abundance_depthzone_site ~ 
                      type + 
                      mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                      mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                      dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                      Substrate_index.s + Substrate_SD.s + Latitude.s)

#subset to only include one temperature variable and one chlorophyll variable
dredge_fish_density_deep <- dredge(global.fish.model_deep, m.lim = c(NA,1) ,subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
         !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
        !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
         !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
         !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
         !(max_sst_C.s && min_sst_C.s) &
         !(mean_sst_C.s && max_sst_C.s) &
         !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_fish_density_deep.dt <- data.table(dredge_fish_density_deep)

dredge_fish_density_deep.dt <- dredge_fish_density_deep.dt[delta < 10]

#model average
dredge_fish_density_deep.modavg <- model.avg(dredge_fish_density_deep, subset = delta < 10)
dredge_fish_density_deep.modavg.coefs <- transpose(data.table(dredge_fish_density_deep.modavg$coefficients[2,]))
colnames(dredge_fish_density_deep.modavg.coefs) <- colnames(dredge_fish_density_deep.modavg$coefficients)

#output row of data table that will make correlation matrix
fish_density_deep_output <- data.table(metric_category = "abundance", taxa = "fish", metric = "density", DepthZone = "Deep")
fish_density_deep_output <- cbind(fish_density_deep_output, dredge_fish_density_deep.modavg.coefs)

#Fish OUTER

dat_fish_total_abundances.env.outer <- dat_fish_total_abundances.env[DepthZone == "Outer",]

#dredge
global.fish.model_outer <- lm(data = dat_fish_total_abundances.env.outer, total_abundance_depthzone_site ~ 
                          type + 
                          mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                          mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                          dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                          Substrate_index.s + Substrate_SD.s + Latitude.s)

#subset to only include one temperature variable and one chlorophyll variable
dredge_fish_density_outer <- dredge(global.fish.model_outer, m.lim = c(NA,1) ,subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                     !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                     !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                     !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                     !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                     !(max_sst_C.s && min_sst_C.s) &
                                     !(mean_sst_C.s && max_sst_C.s) &
                                     !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_fish_density_outer.dt <- data.table(dredge_fish_density_outer)

dredge_fish_density_outer.dt <- dredge_fish_density_outer.dt[delta < 10]

#model average
dredge_fish_density_outer.modavg <- model.avg(dredge_fish_density_outer, subset = delta < 10)
dredge_fish_density_outer.modavg.coefs <- transpose(data.table(dredge_fish_density_outer.modavg$coefficients[2,]))
colnames(dredge_fish_density_outer.modavg.coefs) <- colnames(dredge_fish_density_outer.modavg$coefficients)

#output row of data table that will make correlation matrix
fish_density_outer_output <- data.table(metric_category = "abundance", taxa = "fish", metric = "density", DepthZone = "Outer")
fish_density_outer_output <- cbind(fish_density_outer_output, dredge_fish_density_outer.modavg.coefs)

#Fish MIDDLE

dat_fish_total_abundances.env.middle <- dat_fish_total_abundances.env[DepthZone == "Middle",]

#dredge
global.fish.model_middle <- lm(data = dat_fish_total_abundances.env.middle, total_abundance_depthzone_site ~ 
                           type + 
                           mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                           mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                           dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                           Substrate_index.s + Substrate_SD.s + Latitude.s)

#subset to only include one temperature variable and one chlorophyll variable
dredge_fish_density_middle <- dredge(global.fish.model_middle, m.lim = c(NA,1) ,subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                      !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                      !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                      !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                      !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                      !(max_sst_C.s && min_sst_C.s) &
                                      !(mean_sst_C.s && max_sst_C.s) &
                                      !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_fish_density_middle.dt <- data.table(dredge_fish_density_middle)

dredge_fish_density_middle.dt <- dredge_fish_density_middle.dt[delta < 10]

#model average
dredge_fish_density_middle.modavg <- model.avg(dredge_fish_density_middle, subset = delta < 10)
dredge_fish_density_middle.modavg.coefs <- transpose(data.table(dredge_fish_density_middle.modavg$coefficients[2,]))
colnames(dredge_fish_density_middle.modavg.coefs) <- colnames(dredge_fish_density_middle.modavg$coefficients)

#output row of data table that will make correlation matrix
fish_density_middle_output <- data.table(metric_category = "abundance", taxa = "fish", metric = "density", DepthZone = "Middle")
fish_density_middle_output <- cbind(fish_density_middle_output, dredge_fish_density_middle.modavg.coefs)

#Fish INNER

dat_fish_total_abundances.env.inner <- dat_fish_total_abundances.env[DepthZone == "Inner",]

#dredge
global.fish.model_inner <- lm(data = dat_fish_total_abundances.env.inner, total_abundance_depthzone_site ~ 
                            type + 
                            mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                            mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                            dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                            Substrate_index.s + Substrate_SD.s + Latitude.s)

#subset to only include one temperature variable and one chlorophyll variable
dredge_fish_density_inner <- dredge(global.fish.model_inner, m.lim = c(NA,1) ,subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                       !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                       !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                       !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                       !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                       !(max_sst_C.s && min_sst_C.s) &
                                       !(mean_sst_C.s && max_sst_C.s) &
                                       !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_fish_density_inner.dt <- data.table(dredge_fish_density_inner)

dredge_fish_density_inner.dt <- dredge_fish_density_inner.dt[delta < 10]

#model average
dredge_fish_density_inner.modavg <- model.avg(dredge_fish_density_inner, subset = delta < 10)
dredge_fish_density_inner.modavg.coefs <- transpose(data.table(dredge_fish_density_inner.modavg$coefficients[2,]))
colnames(dredge_fish_density_inner.modavg.coefs) <- colnames(dredge_fish_density_inner.modavg$coefficients)

#output row of data table that will make correlation matrix
fish_density_inner_output <- data.table(metric_category = "abundance", taxa = "fish", metric = "density", DepthZone = "Inner")
fish_density_inner_output <- cbind(fish_density_inner_output, dredge_fish_density_inner.modavg.coefs)

#merge density output
fish_density_output_list <- mget(ls(pattern = "fish_density_.*_output"))

#merge these individual objects
fish_density_output_final <- rbindlist(fish_density_output_list, use.names = T, fill = T)

#################
#FISH BIOMASS
#################

#DEEP
#dredge
global.fish.model_deep_biomass <- lm(data = dat_fish_total_abundances.env.deep, total_biomass_depthzone_site ~ 
                          type + 
                          mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                          mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                          dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                          Substrate_index.s + Substrate_SD.s + Latitude.s)

#subset to only include one temperature variable and one chlorophyll variable
dredge_fish_density_deep_biomass <- dredge(global.fish.model_deep_biomass, m.lim = c(NA,1) ,subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                     !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                     !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                     !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                     !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                     !(max_sst_C.s && min_sst_C.s) &
                                     !(mean_sst_C.s && max_sst_C.s) &
                                     !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_fish_density_deep_biomass.dt <- data.table(dredge_fish_density_deep_biomass)

dredge_fish_density_deep_biomass.dt <- dredge_fish_density_deep_biomass.dt[delta < 20]

#model average
dredge_fish_density_deep_biomass.modavg <- model.avg(dredge_fish_density_deep_biomass, subset = delta < 20)
dredge_fish_density_deep_biomass.modavg.coefs <- transpose(data.table(dredge_fish_density_deep_biomass.modavg$coefficients[2,]))
colnames(dredge_fish_density_deep_biomass.modavg.coefs) <- colnames(dredge_fish_density_deep_biomass.modavg$coefficients)

#output row of data table that will make correlation matrix
fish_density_deep_biomass_output <- data.table(metric_category = "abundance", taxa = "fish", metric = "biomass", DepthZone = "Deep")
fish_density_deep_biomass_output <- cbind(fish_density_deep_biomass_output, dredge_fish_density_deep_biomass.modavg.coefs)

#Fish OUTER

#dredge
global.fish.model_outer_biomass <- lm(data = dat_fish_total_abundances.env.outer, total_biomass_depthzone_site ~ 
                           type + 
                           mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                           mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                           dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                           Substrate_index.s + Substrate_SD.s + Latitude.s)

#subset to only include one temperature variable and one chlorophyll variable
dredge_fish_density_outer_biomass <- dredge(global.fish.model_outer_biomass, m.lim = c(NA,1) ,subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                      !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                      !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                      !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                      !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                      !(max_sst_C.s && min_sst_C.s) &
                                      !(mean_sst_C.s && max_sst_C.s) &
                                      !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_fish_density_outer_biomass.dt <- data.table(dredge_fish_density_outer_biomass)

dredge_fish_density_outer_biomass.dt <- dredge_fish_density_outer_biomass.dt[delta < 20]

#model average
dredge_fish_density_outer_biomass.modavg <- model.avg(dredge_fish_density_outer_biomass, subset = delta < 20)
dredge_fish_density_outer_biomass.modavg.coefs <- transpose(data.table(dredge_fish_density_outer_biomass.modavg$coefficients[2,]))
colnames(dredge_fish_density_outer_biomass.modavg.coefs) <- colnames(dredge_fish_density_outer_biomass.modavg$coefficients)

#output row of data table that will make correlation matrix
fish_density_outer_biomass_output <- data.table(metric_category = "abundance", taxa = "fish", metric = "biomass", DepthZone = "Outer")
fish_density_outer_biomass_output <- cbind(fish_density_outer_biomass_output, dredge_fish_density_outer_biomass.modavg.coefs)

#Fish MIDDLE
#dredge
global.fish.model_middle_biomass <- lm(data = dat_fish_total_abundances.env.middle, total_biomass_depthzone_site ~ 
                            type + 
                            mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                            mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                            dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                            Substrate_index.s + Substrate_SD.s + Latitude.s)

#subset to only include one temperature variable and one chlorophyll variable
dredge_fish_density_middle_biomass <- dredge(global.fish.model_middle_biomass, m.lim = c(NA,1) ,subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                       !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                       !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                       !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                       !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                       !(max_sst_C.s && min_sst_C.s) &
                                       !(mean_sst_C.s && max_sst_C.s) &
                                       !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_fish_density_middle_biomass.dt <- data.table(dredge_fish_density_middle_biomass)

dredge_fish_density_middle_biomass.dt <- dredge_fish_density_middle_biomass.dt[delta < 20]

#model average
dredge_fish_density_middle_biomass.modavg <- model.avg(dredge_fish_density_middle_biomass, subset = delta < 20)
dredge_fish_density_middle_biomass.modavg.coefs <- transpose(data.table(dredge_fish_density_middle_biomass.modavg$coefficients[2,]))
colnames(dredge_fish_density_middle_biomass.modavg.coefs) <- colnames(dredge_fish_density_middle_biomass.modavg$coefficients)

#output row of data table that will make correlation matrix
fish_density_middle_biomass_output <- data.table(metric_category = "abundance", taxa = "fish", metric = "biomass", DepthZone = "Middle")
fish_density_middle_biomass_output <- cbind(fish_density_middle_biomass_output, dredge_fish_density_middle_biomass.modavg.coefs)

#Fish INNER

#dredge
global.fish.model_inner_biomass <- lm(data = dat_fish_total_abundances.env.inner, total_biomass_depthzone_site ~ 
                           type + 
                           mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                           mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                           dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                           Substrate_index.s + Substrate_SD.s + Latitude.s)

#subset to only include one temperature variable and one chlorophyll variable
dredge_fish_density_inner_biomass <- dredge(global.fish.model_inner_biomass, m.lim = c(NA,1) ,subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                      !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                      !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                      !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                      !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                      !(max_sst_C.s && min_sst_C.s) &
                                      !(mean_sst_C.s && max_sst_C.s) &
                                      !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_fish_density_inner_biomass.dt <- data.table(dredge_fish_density_inner_biomass)

dredge_fish_density_inner_biomass.dt <- dredge_fish_density_inner_biomass.dt[delta < 20]

#model average
dredge_fish_density_inner_biomass.modavg <- model.avg(dredge_fish_density_inner_biomass, subset = delta < 20)
dredge_fish_density_inner_biomass.modavg.coefs <- transpose(data.table(dredge_fish_density_inner_biomass.modavg$coefficients[2,]))
colnames(dredge_fish_density_inner_biomass.modavg.coefs) <- colnames(dredge_fish_density_inner_biomass.modavg$coefficients)

#output row of data table that will make correlation matrix
fish_density_inner_biomass_output <- data.table(metric_category = "abundance", taxa = "fish", metric = "biomass", DepthZone = "Inner")
fish_density_inner_biomass_output <- cbind(fish_density_inner_biomass_output, dredge_fish_density_inner_biomass.modavg.coefs)

#merge density output
fish_density_biomass_output_list <- mget(ls(pattern = "fish_density_.*biomass_output"))

#merge these individual objects
fish_density_biomass_output_final <- rbindlist(fish_density_biomass_output_list, use.names = T, fill = T)

#flip wide to long
fish_density_biomass_output_final.l <- melt(fish_density_biomass_output_final, id.vars = c("metric_category","taxa","metric","DepthZone"), variable.name = "model_parameter",value.name = "coefficient")

fish_density_biomass_output_final.l[,DepthZone := factor(DepthZone, levels = c("Deep","Outer","Middle","Inner"))]

ggplot(fish_density_biomass_output_final.l) +
  geom_tile(aes(y = DepthZone, x = model_parameter, fill = coefficient)) +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  scale_x_discrete(position = "top", labels = c("Intercept: Island","Intercept: Mainland","Mean chlorophyll","Relief","Maximum chlorophyll","Substrate SD","Relief SD","Distance from 200m isobath","Minimum chlorophyll","Substrate","Giant kelp stipe density","Minimum SST","Latitude","Mean SST","Maximum SST"), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = -0.1))

#################################
#MODEL MACROINVERT ABUNDANCE/DENSITY
#################################
#scale all variables
# Identify numeric columns
numeric_cols <- names(dat_macroinvert_total_abundances.env)[sapply(dat_macroinvert_total_abundances.env, is.numeric) & !(names(dat_macroinvert_total_abundances.env) %in% c("ID", "true_label"))]

# Scale numeric columns and create new columns with .s appended
dat_macroinvert_total_abundances.env[, (paste0(numeric_cols, ".s")) := lapply(.SD, scale), .SDcols = numeric_cols]

#exclude two sites for which we don't have data for chlorophyll or sst
dat_macroinvert_total_abundances.env <- dat_macroinvert_total_abundances.env[complete.cases(dat_macroinvert_total_abundances.env),]

#dredge
global.model <- lme(data = dat_macroinvert_total_abundances.env, total_abundance_depthzone_site ~ 
                      #  DepthZone + 
                      type + 
                      mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                      mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                      dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                      Substrate_index.s + Substrate_SD.s + Latitude.s, random = ~1|Site)

#subset to only include one temperature variable and one chlorophyll variable
dredge_macroinvert_density <- dredge(global.model, subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                !(max_sst_C.s && min_sst_C.s) &
                                !(mean_sst_C.s && max_sst_C.s) &
                                !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_macroinvert_density.dt <- data.table(dredge_macroinvert_density)

dredge_macroinvert_density.dt <- dredge_macroinvert_density.dt[delta < 2]

#best model
best_macroinvert_density_lme <- lme(data = dat_macroinvert_total_abundances.env, total_abundance_depthzone_site ~ 
                                      type +
                                      dist_200m_bath.s +
                                      mean_sst_C.s +
                                      giantkelp_stipe_density_m2 +
                                      Relief_SD.s +
                                      Substrate_index.s, random = ~1|Site)

summary(best_macroinvert_density_lme)      
r.squaredGLMM(best_macroinvert_density_lme)

#R^2 = 0.30, with random effect = 0.76
#density increases with distance to shelf break
#density increases with substrate SD
#density increases with SD of relief
#density decreases with giant kelp stipe density
#density increases wtih latitude

#MODEL KELP ABUNDANCE/DENSITY

#scale all variables
# Identify numeric columns
numeric_cols <- names(dat_kelp_total_abundances.env)[sapply(dat_kelp_total_abundances.env, is.numeric) & !(names(dat_kelp_total_abundances.env)  %in% c("ID", "true_label"))]

# Scale numeric columns and create new columns with .s appended
dat_kelp_total_abundances.env[, (paste0(numeric_cols, ".s")) := lapply(.SD, scale), .SDcols = numeric_cols]

#exclude two sites for which we don't have data for chlorophyll or sst
dat_kelp_total_abundances.env <- dat_kelp_total_abundances.env[complete.cases(dat_kelp_total_abundances.env),]

#dredge
global.model <- lme(data = dat_kelp_total_abundances.env, total_abundance_depthzone_site ~ 
                      #  DepthZone + 
                      type + 
                      mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                      mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                      dist_200m_bath.s +  Relief_index.s + Relief_SD.s +
                      Substrate_index.s + Substrate_SD.s + Latitude.s, random = ~1|Site)

#subset to only include one temperature variable and one chlorophyll variable
dredge_kelp_density <- dredge(global.model, subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                       !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                       !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                       !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                       !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                       !(max_sst_C.s && min_sst_C.s) &
                                       !(mean_sst_C.s && max_sst_C.s) &
                                       !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_kelp_density.dt <- data.table(dredge_kelp_density)

dredge_kelp_density.dt <- dredge_kelp_density.dt[delta < 2]

#best model
best_kelp_density_lme <- lme(data = dat_kelp_total_abundances.env, total_abundance_depthzone_site ~ 
                               type +
                               Relief_index.s +
                               Substrate_index.s
                               , random = ~1|Site)

summary(best_kelp_density_lme)      
r.squaredGLMM(best_kelp_density_lme)

#R^2 = 0.14, with random effect = 0.47
#density decreases with mean cholorphyl
#density increases with mean temperature
#density decreases with relief SD
#density increases with substrate size (index)

#not intuitive...
#################
######MODEL FISH RICHNESS

#scale all variables
# Identify numeric columns
numeric_cols <- names(dat_fish_richness_depthzone_site.env)[sapply(dat_fish_richness_depthzone_site.env, is.numeric) & !(names(dat_fish_richness_depthzone_site.env) %in% c("ID", "true_label"))]

# Scale numeric columns and create new columns with .s appended
dat_fish_richness_depthzone_site.env[, (paste0(numeric_cols, ".s")) := lapply(.SD, scale), .SDcols = numeric_cols]

#exclude two sites for which we don't have data for chlorophyll or sst
dat_fish_richness_depthzone_site.env <- dat_fish_richness_depthzone_site.env[complete.cases(dat_fish_richness_depthzone_site.env),]

options(na.action = "na.fail")

#dredge
global.model <- lme(data = dat_fish_richness_depthzone_site.env, total_richness_depthzone_site ~ 
                      #  DepthZone + 
                      type + 
                      mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                      mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                      dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                      Substrate_index.s + Substrate_SD.s + Latitude.s, random = ~1|Site)

#subset to only include one temperature variable and one chlorophyll variable
dredge_fish_richness <- dredge(global.model, subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                !(max_sst_C.s && min_sst_C.s) &
                                !(mean_sst_C.s && max_sst_C.s) &
                                !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_fish_richness.dt <- data.table(dredge_fish_richness)

dredge_fish_richness.dt <- dredge_fish_richness.dt[delta < 2]

#best model
best_fish_richness_lme <- lme(data = dat_fish_richness_depthzone_site.env, total_richness_depthzone_site ~ 
                               type + 
                              Latitude.s +
                               giantkelp_stipe_density_m2 +
                                mean_chl_mg_m3.s +
                                Relief_SD.s, 
                                random = ~1|Site)

summary(best_fish_richness_lme)
r.squaredGLMM(best_fish_richness_lme)

#R^2 = 0.38, with random effect = 0.52
#baseline natural mainland has higher richness
#richness decreases with substrate size
#richness increases with relief and with SD of relief
#richness increases with SST
#richness decreases with giant kelp stipe density
#richness increases wtih latitude

#MODEL MACROINVERT RICHNESS

#scale all variables
# Identify numeric columns
numeric_cols <- names(dat_macroinvert_richness_depthzone_site.env)[sapply(dat_macroinvert_richness_depthzone_site.env, is.numeric) & !(names(dat_macroinvert_richness_depthzone_site.env) %in% c("ID", "true_label"))]

# Scale numeric columns and create new columns with .s appended
dat_macroinvert_richness_depthzone_site.env[, (paste0(numeric_cols, ".s")) := lapply(.SD, scale), .SDcols = numeric_cols]

#exclude two sites for which we don't have data for chlorophyll or sst
dat_macroinvert_richness_depthzone_site.env <- dat_macroinvert_richness_depthzone_site.env[complete.cases(dat_macroinvert_richness_depthzone_site.env),]

#dredge
global.model <- lme(data = dat_macroinvert_richness_depthzone_site.env, total_richness_depthzone_site ~ 
                      #  DepthZone + 
                      type + 
                      mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                      mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                      dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                      Substrate_index.s + Substrate_SD.s + Latitude.s, random = ~1|Site)

#subset to only include one temperature variable and one chlorophyll variable
dredge_macroinvert_richness <- dredge(global.model, subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                       !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                       !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                       !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                       !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                       !(max_sst_C.s && min_sst_C.s) &
                                       !(mean_sst_C.s && max_sst_C.s) &
                                       !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_macroinvert_richness.dt <- data.table(dredge_macroinvert_richness)

dredge_macroinvert_richness.dt <- dredge_macroinvert_richness.dt[delta < 2]

#best model
best_macroinvert_richness_lme <- lme(data = dat_macroinvert_richness_depthzone_site.env, total_richness_depthzone_site ~ 
                                      type +
                                      giantkelp_stipe_density_m2 +
                                      Latitude.s +
                                      mean_sst_C.s +
                                      min_chl_mg_m3.s +
                                      Substrate_index.s +
                                      Relief_SD.s +
                                      Substrate_SD.s, random = ~1|Site)

summary(best_macroinvert_richness_lme)      
r.squaredGLMM(best_macroinvert_richness_lme)

#R^2 = 0.48, with random effect = 0.57
#richness higher at mainland than on islands
#richness decreases with giant kelp density
#richness increases with latitude
#richness decreases with mean temp
#richness decreases with minimum chlorophyll
#richness decreases with substrate size (index)
#richness increases with relief SD
#richness increases with substrate SD

#MODEL KELP RICHNESS

#scale all variables
# Identify numeric columns
numeric_cols <- names(dat_kelp_richness_depthzone_site.env)[sapply(dat_kelp_richness_depthzone_site.env, is.numeric) & !(names(dat_kelp_richness_depthzone_site.env) %in% c("ID", "true_label"))]

# Scale numeric columns and create new columns with .s appended
dat_kelp_richness_depthzone_site.env[, (paste0(numeric_cols, ".s")) := lapply(.SD, scale), .SDcols = numeric_cols]

#exclude two sites for which we don't have data for chlorophyll or sst
dat_kelp_richness_depthzone_site.env <- dat_kelp_richness_depthzone_site.env[complete.cases(dat_kelp_richness_depthzone_site.env),]

#dredge
global.model <- lme(data = dat_kelp_richness_depthzone_site.env, total_richness_depthzone_site ~ 
                      #  DepthZone + 
                      type + 
                      mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                      mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                      dist_200m_bath.s + Relief_index.s + Relief_SD.s +
                      Substrate_index.s + Substrate_SD.s + Latitude.s, random = ~1|Site)

#subset to only include one temperature variable and one chlorophyll variable
dredge_kelp_richness <- dredge(global.model, subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                        !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                        !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                        !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                        !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                        !(max_sst_C.s && min_sst_C.s) &
                                        !(mean_sst_C.s && max_sst_C.s) &
                                        !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_kelp_richness.dt <- data.table(dredge_kelp_richness)

dredge_kelp_richness.dt <- dredge_kelp_richness.dt[delta < 2]

#best model
best_kelp_richness_lme <- lme(data = dat_kelp_richness_depthzone_site.env, total_richness_depthzone_site ~ 
                                       type +
                                       giantkelp_stipe_density_m2 +
                                       Latitude.s +
                                       mean_sst_C.s +
                                       min_chl_mg_m3.s +
                                       Substrate_index.s +
                                       Relief_SD.s +
                                       Substrate_SD.s, random = ~1|Site)

summary(best_kelp_richness_lme)      
r.squaredGLMM(best_kelp_richness_lme)

#R^2 = 0.11, with random effect = 0.44

#################
######MODEL FISH Simpson Diversity

#scale all variables
# Identify numeric columns
numeric_cols <- names(dat_fish_simpson.env)[sapply(dat_fish_simpson.env, is.numeric) & !(names(dat_fish_simpson.env)  %in% c("ID", "true_label"))]

# Scale numeric columns and create new columns with .s appended
dat_fish_simpson.env[, (paste0(numeric_cols, ".s")) := lapply(.SD, scale), .SDcols = numeric_cols]

#exclude two sites for which we don't have data for chlorophyll or sst
dat_fish_simpson.env <- dat_fish_simpson.env[complete.cases(dat_fish_simpson.env),]

#dredge
global.model <- lme(data = dat_fish_simpson.env, Simpson_index_density ~ 
                      #  DepthZone + 
                      type + 
                      mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                      mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                      dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                      Substrate_index.s + Substrate_SD.s + Latitude.s, random = ~1|Site)

#subset to only include one temperature variable and one chlorophyll variable
dredge_fish_diversity <- dredge(global.model, subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                 !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                 !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                 !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                 !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                 !(max_sst_C.s && min_sst_C.s) &
                                 !(mean_sst_C.s && max_sst_C.s) &
                                 !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_fish_diversity.dt <- data.table(dredge_fish_diversity)

dredge_fish_diversity.dt <- dredge_fish_diversity.dt[delta < 2]

#best model
best_fish_diversity_lme <- lme(data = dat_fish_simpson.env, Simpson_index_density ~ 
                               giantkelp_stipe_density_m2 + 
                              type,
                              random = ~1|Site)

summary(best_fish_diversity_lme)
r.squaredGLMM(best_fish_diversity_lme)

#R^2 = 0.27, with random effect = 0.40
#baseline natural mainland has higher diversity
#diversity decreases with substrate size
#diversity increases with relief and with SD of relief
#diversity increases with SST
#diversity decreases with giant kelp stipe density
#diversity increases wtih latitude

#MODEL MACROINVERT diversity

#scale all variables
# Identify numeric columns
numeric_cols <- names(dat_macroinvert_simpson.env)[sapply(dat_macroinvert_simpson.env, is.numeric) & !(names(dat_macroinvert_simpson.env) %in% c("ID", "true_label"))]

# Scale numeric columns and create new columns with .s appended
dat_macroinvert_simpson.env[, (paste0(numeric_cols, ".s")) := lapply(.SD, scale), .SDcols = numeric_cols]

#exclude two sites for which we don't have data for chlorophyll or sst
dat_macroinvert_simpson.env <- dat_macroinvert_simpson.env[complete.cases(dat_macroinvert_simpson.env),]

#dredge
global.model <- lme(data = dat_macroinvert_simpson.env, Simpson_index_density ~ 
                      #  DepthZone + 
                      type + 
                      mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                      mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                      dist_200m_bath.s + giantkelp_stipe_density_m2 + Relief_index.s + Relief_SD.s +
                      Substrate_index.s + Substrate_SD.s + Latitude.s, random = ~1|Site)

#subset to only include one temperature variable and one chlorophyll variable
dredge_macroinvert_diversity <- dredge(global.model, subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                        !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                        !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                        !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                        !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                        !(max_sst_C.s && min_sst_C.s) &
                                        !(mean_sst_C.s && max_sst_C.s) &
                                        !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_macroinvert_diversity.dt <- data.table(dredge_macroinvert_diversity)

dredge_macroinvert_diversity.dt <- dredge_macroinvert_diversity.dt[delta < 2]

#best model
best_macroinvert_diversity_lme <- lme(data = dat_macroinvert_simpson.env, Simpson_index_density ~ 
                                        giantkelp_stipe_density_m2 +
                                        mean_sst_C.s
                                        , random = ~1|Site)

summary(best_macroinvert_diversity_lme)      
r.squaredGLMM(best_macroinvert_diversity_lme)

#R^2 = 0.13, with random effect = 0.56
#diversity higher at mainland than on islands
#diversity decreases with giant kelp density
#diversity increases with latitude
#diversity decreases with mean temp
#diversity decreases with minimum chlorophyll
#diversity decreases with substrate size (index)
#diversity increases with relief SD
#diversity increases with substrate SD

#MODEL KELP diversity

#scale all variables
# Identify numeric columns
numeric_cols <- names(dat_kelp_simpson.env)[sapply(dat_kelp_simpson.env, is.numeric) & !(names(dat_kelp_simpson.env)  %in% c("ID", "true_label"))]

# Scale numeric columns and create new columns with .s appended
dat_kelp_simpson.env[, (paste0(numeric_cols, ".s")) := lapply(.SD, scale), .SDcols = numeric_cols]

#exclude two sites for which we don't have data for chlorophyll or sst
dat_kelp_simpson.env <- dat_kelp_simpson.env[complete.cases(dat_kelp_simpson.env),]

#dredge
global.model <- lme(data = dat_kelp_simpson.env, Simpson_index_density ~ 
                      #  DepthZone + 
                      type + 
                      mean_chl_mg_m3.s + max_chl_mg_m3.s + min_chl_mg_m3.s +
                      mean_sst_C.s + max_sst_C.s + min_sst_C.s +
                      dist_200m_bath.s +  Relief_index.s + Relief_SD.s +
                      Substrate_index.s + Substrate_SD.s + Latitude.s, random = ~1|Site)

#subset to only include one temperature variable and one chlorophyll variable
dredge_kelp_diversity <- dredge(global.model, subset = !(mean_chl_mg_m3.s && max_chl_mg_m3.s && min_chl_mg_m3.s ) &
                                 !(mean_chl_mg_m3.s && max_chl_mg_m3.s) &
                                 !(max_chl_mg_m3.s && min_chl_mg_m3.s) &
                                 !(mean_chl_mg_m3.s && min_chl_mg_m3.s) &
                                 !(mean_sst_C.s && max_sst_C.s && min_sst_C.s) &
                                 !(max_sst_C.s && min_sst_C.s) &
                                 !(mean_sst_C.s && max_sst_C.s) &
                                 !(mean_sst_C.s && min_sst_C.s))

#table with top models
dredge_kelp_diversity.dt <- data.table(dredge_kelp_diversity)

dredge_kelp_diversity.dt <- dredge_kelp_diversity.dt[delta < 2]

#best model
best_kelp_diversity_lme <- lme(data = dat_kelp_simpson.env, Simpson_index_density ~ 
                                Substrate_index.s, random = ~1|Site)

summary(best_kelp_diversity_lme)      
r.squaredGLMM(best_kelp_diversity_lme)

#R^2 = 0.07, with random effect = 0.33
#diversity decreases with distance from shelf break
#diversity increases with max chlorophyl
#diversity decreases with SD of substrate
#diversity is baseline higher at islands
