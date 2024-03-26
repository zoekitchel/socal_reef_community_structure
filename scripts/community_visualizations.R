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
library(cowplot)
library(ggnewscale)

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

 ##################################################
#Long data to wide data for vegan analyses
##################################################

#melt long to wide
dat_fish_averages_bysite.wide <- dcast(dat_fish_long_density, Region + Site + DepthZone ~ Species, value.var = "mean_density_m2", fun = mean)

dat_macroinvert_averages_bysite.wide <- dcast(dat_macroinvert_long_density, Region + Site + DepthZone ~ BenthicReefSpecies, value.var = "mean_density_m2", fun = mean)

dat_kelp_averages_bysite.wide <- dcast(dat_kelp_long_density, Region + Site + DepthZone ~ BenthicReefSpecies, value.var = "mean_density_m2", fun = mean)

#check that all ordered in same way
stopifnot(dat_fish_averages_bysite.wide[,1:3]==dat_kelp_averages_bysite.wide[,1:3])
stopifnot(dat_fish_averages_bysite.wide[,1:3]==dat_macroinvert_averages_bysite.wide[,1:3])

#Merge for all species visualization
dat_averages_bysite.wide <- cbind(dat_fish_averages_bysite.wide,
                                  dat_kelp_averages_bysite.wide[,4:ncol(dat_kelp_averages_bysite.wide)],
                                  dat_macroinvert_averages_bysite.wide[,4:ncol(dat_macroinvert_averages_bysite.wide)])

####Add variables!!
#dat_averages_bysite.wide.envir <- VRG_lat_lon_only[dat_averages_bysite.wide, on = c("Site","DepthZone")]
#dat_averages_bysite.wide.envir <- distance_200mbathy_bysite[dat_averages_bysite.wide.envir, on = c("Site","DepthZone")]
#dat_averages_bysite.wide.envir <- macro_density_bysite[dat_averages_bysite.wide.envir, on = c("Site","DepthZone")]

#####################
#PERMANOVA, Permutational Multivariate Analysis of Variance (perMANOVA)
#####################
#FULL COMMUNITY

#Start with root transformation (high abundance spp less pull in ordination)
dat_averages_bysite.wide.rt <- sqrt(dat_averages_bysite.wide[,c(4:ncol(dat_averages_bysite.wide)), with = FALSE])


#add back reference columns
dat_averages_bysite.wide.rt <- cbind(dat_averages_bysite.wide[,c(1:3), with = FALSE], dat_averages_bysite.wide.rt)

    #alternative version without artificial reefs
    dat_averages_bysite_NOARM.wide.rt <- dat_averages_bysite.wide.rt[DepthZone != "ARM",]

#only keep Species
dat_averages_bysite.wide.rt.trim <- dat_averages_bysite.wide.rt[,c(4:ncol(dat_averages_bysite.wide.rt)), with = FALSE]

  #alternative version without artificial reefs
  dat_averages_bysite_NOARM.wide.rt.trim <- dat_averages_bysite_NOARM.wide.rt[,c(4:ncol(dat_averages_bysite_NOARM.wide.rt)), with = FALSE]

#permanova including ARM as depth zone
permanova_allspp <- adonis2(
  dat_averages_bysite.wide.rt.trim ~ dat_averages_bysite.wide.rt$DepthZone,
  method = "bray"
)

permanova_allspp

#permanova not including ARM as depth zone
permanova_allspp_noarm <- adonis2(
  dat_averages_bysite_NOARM.wide.rt.trim ~ dat_averages_bysite_NOARM.wide.rt$DepthZone,
  method = "bray"
)
permanova_allspp_noarm

permanova_allspp_noarm <- adonis2(
  dat_averages_bysite_NOARM.wide.rt.trim ~ dat_averages_bysite_NOARM.wide.rt$DepthZone +
    dat_averages_bysite_NOARM.wide.rt$macro_mean_density_m2 +
    dat_averages_bysite_NOARM.wide.rt$dist_200m_bath +
    dat_averages_bysite_NOARM.wide.rt$BO_sstmean +
    dat_averages_bysite_NOARM.wide.rt$Region,
  method = "bray"
)

permanova_allspp_noarm

#depth zone accounts for 19% of variation
#depth zone PLUS AR accounts for 23% of variation

#ONLY FISH

#Start with root transformation (high abundance spp less pull in ordination)
dat_fish_averages_bysite.wide.rt <- sqrt(dat_fish_averages_bysite.wide[,c(4:ncol(dat_fish_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_fish_averages_bysite.wide.rt <- cbind(dat_fish_averages_bysite.wide[,c(1:3), with = FALSE], dat_fish_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_fish_averages_bysite.wide.rt.trim <- dat_fish_averages_bysite.wide.rt[,c(4:ncol(dat_fish_averages_bysite.wide.rt)), with = FALSE]

permanova_fish <- adonis2(
  dat_fish_averages_bysite.wide.rt.trim ~ dat_fish_averages_bysite.wide.rt$DepthZone,
  method = "bray"
)

permanova_fish

#fish names
fish_names <- colnames(dat_fish_averages_bysite.wide.rt.trim)

#depth zone affects fish community composition, accounting for 23% of variation in composition

#ONLY MACRO

#Start with root transformation (high abundance spp less pull in ordination)
dat_macroinvert_averages_bysite.wide.rt <- sqrt(dat_macroinvert_averages_bysite.wide[,c(4:ncol(dat_macroinvert_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_macroinvert_averages_bysite.wide.rt <- cbind(dat_macroinvert_averages_bysite.wide[,c(1:3), with = FALSE], dat_macroinvert_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_macroinvert_averages_bysite.wide.rt.trim <- dat_macroinvert_averages_bysite.wide.rt[,c(4:ncol(dat_macroinvert_averages_bysite.wide.rt)), with = FALSE]

permanova_macroinvert <- adonis2(
  dat_macroinvert_averages_bysite.wide.rt.trim ~ dat_macroinvert_averages_bysite.wide.rt$DepthZone,
  method = "bray"
)

permanova_macroinvert

#macro names
macro_names <- colnames(dat_macroinvert_averages_bysite.wide.rt.trim)

#depth zone affects macroinvert community composition, but only accounting for 18% of variation in composition


#ONLY KELP

#sometimes, no kelp at all, need to delete these rows
kelp_row_sums <- rowSums(dat_kelp_averages_bysite.wide[,5:ncol(dat_kelp_wide_density)])

dat_kelp_averages_bysite.wide[,rowSums := kelp_row_sums]

dat_kelp_averages_bysite.wide.r <- dat_kelp_averages_bysite.wide[rowSums>0]

#Start with root transformation (high abundance spp less pull in ordination)
dat_kelp_averages_bysite.wide.rt <- sqrt(dat_kelp_averages_bysite.wide.r[,c(4:ncol(dat_kelp_averages_bysite.wide)), with = FALSE])

#add back reference columns
dat_kelp_averages_bysite.wide.rt <- cbind(dat_kelp_averages_bysite.wide.r[,c(1:3), with = FALSE], dat_kelp_averages_bysite.wide.rt)

#only keep Site and DepthZone
dat_kelp_averages_bysite.wide.rt.trim <- dat_kelp_averages_bysite.wide.rt[,c(4:ncol(dat_kelp_averages_bysite.wide.rt)), with = FALSE]

permanova_kelp <- adonis2(
  dat_kelp_averages_bysite.wide.rt.trim ~ dat_kelp_averages_bysite.wide.rt$DepthZone,
  method = "bray"
)

permanova_kelp

#kelp names
kelp_names <- colnames(dat_kelp_averages_bysite.wide.rt.trim[,1:20])

#whether or not a site is an artificial reef affects fish community composition, accounting for 21% of variation in composition

#Take away:
#


###################
#PLOT NMDS
###################
#NMDS “preserves the rank order of the inter-point dissimilarities as well as possible within the constraints of a small number of dimensions” (Anderson et al. 2008, p. 105).
#It has the least restrictive assumptions (any distance measure, no assumptions about shape of species response curves).

#all species

#convert to distance matrix
dat_averages_bysite.dist <- vegdist(dat_averages_bysite.wide[,4:ncol(dat_averages_bysite.wide)], distance = "bray")

full_nmds <- metaMDS(dat_averages_bysite.wide[,4:ncol(dat_averages_bysite.wide)], #metaMDS automatically does this for us
                     autotransform = FALSE,
                     distance = "bray",
                     engine = "monoMDS",
                     wascores = T,
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
#         k = 2, stress = 0.21
#         k = 3, stress = 0.14
#         k = 4, stress = 0.10, good ordination with no real risk of drawing false inferences (Clark 1993 pg 126)
#         k = 5, stress = 0.09
#         k = 6, stress = 0.08
#         k = 7, stress = 0.07

#The fit of a NMDS ordination can be assessed by plotting the original dissimilarities (z$diss) against the (Euclidean) ordination distances (z$dist) 

stressplot(full_nmds)

#link NMDS points to site details
dat_averages_bysite.nmds <- cbind(dat_averages_bysite.wide[,1:3], full_nmds$points)

#For some reason, Marina Del Rey coded as Region = Artificial reef, change to Santa Monica Bay

dat_averages_bysite.nmds[Region == "Artificial Reef",Region := "Santa Monica Bay"]

#Pull out species scores
species_scores <- data.table(vegan::scores(full_nmds,"species"))
species_scores[,Species := rownames(vegan::scores(full_nmds,"species"))]

#add category
species_scores[,`Spp category`:= c(rep("Fish",length(fish_names)), rep("Kelp",length(kelp_names)), rep("Macroinvert",length(macro_names)))]

#put through spp function to get common name
spp_list <- get_taxa(taxon_list = unique(species_scores$Species))

#top 5 values
top_5 <- c("Pugettia producta", "Micrometers minimus", "Hypterprosopon argenteum", "Syngnathus leptorhynchus", "Triopha catalinae",
           "Phanerodon atripes", "Cadlina flavomaculata", "Sebastes paucispinis", "Holothuria zacae", "Peltodoris mullineri", "Sargassum Palmer",
           "Undaria pinnatifid", "Alaria marginata", "Flabellina pricei", "Sebastes constellatus", "Calliostoma annulatum", "Janolus barbarensis",
           "Rathbunella hypoplecta", "Acanthodoris rhodoceras", "Cancer productus","Taliepus nuttallii","Cymatogaster aggregata","Peltodoris mullineri",
           "Sargassum horneri","Sargassum palmeri","Macrocystis pyrifera")

full_nmds_plot_bysitetype <- ggplot() +
  geom_point(data = dat_averages_bysite.nmds, aes(x = MDS1, y = MDS2, color = DepthZone, shape = DepthZone), size = 3) +
  stat_ellipse(data = dat_averages_bysite.nmds, aes(x = MDS1, y = MDS2, color = DepthZone)) + #95% confidence level for a mlutivariate t-distribution
  scale_color_manual(values = c("cornflowerblue","coral1","chartreuse3","darkorchid2","black")) +
  scale_shape_manual(values = c(15,16,17,18,12)) +
 # new_scale_color() +
 # ggrepel::geom_text_repel(data = species_scores[Species %in% top_5], aes(x = NMDS1, y = NMDS2, label = Species, color = `Spp category`), size = 2.5, fontface ="bold.italic") +
 # scale_color_manual(values = c("slategrey","black","brown")) +
  annotate(geom = "text", x = 1.4, y = 1.8, label = "4d Stress: 0.11", fontface = "bold") +
  theme_classic() +
#  lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme(
    axis.title = element_blank(),
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        legend.position = c(.1, .83),
        legend.background = element_blank(), legend.key = element_blank()
     )


ggsave(full_nmds_plot_bysitetype, path = "figures", filename = "full_nmds_plot_bysitetype.jpg", height = 6, width = 8)

#species labels with color
full_nmds_plot_bysitetype_spplabel_color <- ggplot() +
  geom_point(data = dat_averages_bysite.nmds, aes(x = MDS1, y = MDS2, color = DepthZone, shape = DepthZone), size = 3) +
  stat_ellipse(data = dat_averages_bysite.nmds, aes(x = MDS1, y = MDS2, color = DepthZone)) + #95% confidence level for a mlutivariate t-distribution
  scale_color_manual(values = c("cornflowerblue","coral1","chartreuse3","darkorchid2","black")) +
  scale_shape_manual(values = c(15,16,17,18,12)) +
  new_scale_color() +
  ggrepel::geom_text_repel(data = species_scores[Species %in% top_5], aes(x = NMDS1, y = NMDS2, label = Species, color = `Spp category`), size = 2.5, fontface ="bold.italic") +
  scale_color_manual(values = c("slategrey","black","brown")) +
  annotate(geom = "text", x = 1.8, y = 1.8, label = "4d Stress: 0.11", fontface = "bold") +
  theme_classic() +
  #  lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    #     legend.position = c(.2, .83),
    legend.background = element_blank(), legend.key = element_blank()
  )


ggsave(full_nmds_plot_bysitetype_spplabel_color, path = "figures", filename = "full_nmds_plot_bysitetype_spplabel_color.jpg", height = 6, width = 8)

#species label black and white
full_nmds_plot_bysitetype_spplabel_bw <- ggplot() +
  geom_point(data = dat_averages_bysite.nmds, aes(x = MDS1, y = MDS2, color = DepthZone, shape = DepthZone), size = 3) +
  stat_ellipse(data = dat_averages_bysite.nmds, aes(x = MDS1, y = MDS2, color = DepthZone)) + #95% confidence level for a mlutivariate t-distribution
  scale_color_manual(values = c("cornflowerblue","coral1","chartreuse3","darkorchid2","black")) +
  scale_shape_manual(values = c(15,16,17,18,12)) +
  ggrepel::geom_text_repel(data = species_scores[Species %in% top_5], aes(x = NMDS1, y = NMDS2, label = Species), size = 2.5, fontface ="bold.italic") +
  annotate(geom = "text", x = 1.8, y = 1.8, label = "4d Stress: 0.11", fontface = "bold") +
  theme_classic() +
  #  lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    #     legend.position = c(.2, .83),
    legend.background = element_blank(), legend.key = element_blank()
  )


ggsave(full_nmds_plot_bysitetype_spplabel_bw, path = "figures", filename = "full_nmds_plot_bysitetype_spplabel_bw.jpg", height = 6, width = 8)


    full_nmds_plot_byregion <- ggplot(data = dat_averages_bysite.nmds) +
  geom_point(aes(x = MDS1, y = MDS2, color = Region), size = 2) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black", "#A65628" ,"#F781BF")) +
 # lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.background = element_blank(), legend.key = element_blank())


ggsave(full_nmds_plot_byregion, path = "figures", filename = "full_nmds_plot_byregion.jpg", height = 7, width = 7)

full_nmds_plot_byregion_zone <- ggplot(data = dat_averages_bysite.nmds) +
  geom_point(aes(x = MDS1, y = MDS2, color = Region, shape = DepthZone), size = 2) +
  scale_shape_manual(guide = "none", values = c(15,16,17,18,12)) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black", "#A65628" ,"#F781BF")) +
  # lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.background = element_blank(), legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=4)))

#extract legend for shape only
depthzone_legend <- get_legend(ggplot(data = dat_averages_bysite.nmds) +
  geom_point(aes(x = MDS1, y = MDS2, shape = DepthZone), size = 2) +
    scale_shape_manual(values = c(15,16,17,18,12)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.background = element_blank(), legend.key = element_blank()) +
         guides(shape = guide_legend(override.aes = list(size=4))))

#add shape legend to plot
full_nmds_plot_byregion_zone_legend <- ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  draw_plot(full_nmds_plot_byregion_zone, x = 0, y = 0, height = 1, width = 1) +
  draw_plot(depthzone_legend, x = 0.18, y = 0.7, height = 0.2, width = 0.2)

ggsave(full_nmds_plot_byregion_zone_legend, path = "figures", filename = "full_nmds_plot_byregion_zone_legend.jpg", height = 7, width = 7)

#Visualize by kelp density
#link kelp density
dat_averages_bysite.nmds_pluskelp <- dat_averages_bysite.nmds[macro_density_bysite, on = c("Site","DepthZone")]

full_nmds_plot_byregion_withkelp <- ggplot(data = dat_averages_bysite.nmds_pluskelp) +
  geom_point(aes(x = MDS1, y = MDS2, color = Region, size = macro_mean_density_m2*100)) +
  scale_size_continuous(guide = "none") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black", "#A65628" ,"#F781BF")) +
  # lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.background = element_blank(), legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=4)))

full_nmds_plot_byregion_withkelp 

#extract legend for shape only
macrocystis_density_legend <- get_legend(ggplot(data = dat_averages_bysite.nmds_pluskelp) +
                                           geom_point(aes(x = MDS1, y = MDS2,size = macro_mean_density_m2*100)) +
                                           scale_size_continuous(expression(paste("Average giant kelp density per 100m" ^2))) +
                                           # lims(x= c(-1.4,1.2), y = c(-2,2)) +
                                           theme_classic() +
                                           guides(size = guide_legend(title.position = "top")) +
                                           theme(axis.title = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 axis.text = element_blank(),
                                                 legend.position = "top",
                                                 legend.direction = "horizontal",
                                                 legend.background = element_blank(), legend.key = element_blank()))

#add shape legend to plot
full_nmds_plot_byregion_withkelp_legend <- ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  draw_plot(full_nmds_plot_byregion_withkelp, x = 0, y = 0, height = 1, width = 1) +
  draw_plot(macrocystis_density_legend, x = 0.17, y = 0.7, height = 0.2, width = 0.2)

ggsave(full_nmds_plot_byregion_withkelp_legend, path = "figures", filename = "full_nmds_plot_byregion_withkelp_legend.jpg", height = 7, width = 7)



#Visualize by distance to 200m isobath
#link kelp density
dat_averages_bysite.nmds_dist_200m <- distance_200mbathy_bysite[dat_averages_bysite.nmds, on = c("Site","DepthZone")]

full_nmds_plot_byregion_dist_200m <- ggplot(data = dat_averages_bysite.nmds_dist_200m) +
  geom_point(aes(x = MDS1, y = MDS2, color = Region, size = dist_200m_bath/1000)) +
  scale_size_continuous(guide = "none") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black", "#A65628" ,"#F781BF")) +
  # lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.background = element_blank(), legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=4)))

full_nmds_plot_byregion_dist_200m 

#extract legend for shape only
dist_200m_legend <- get_legend(ggplot(data = dat_averages_bysite.nmds_dist_200m) +
                                           geom_point(aes(x = MDS1, y = MDS2,size = dist_200m_bath/1000)) +
                                           scale_size_continuous("Kilometers to 200m isobath") +
                                           # lims(x= c(-1.4,1.2), y = c(-2,2)) +
                                           theme_classic() +
                                           guides(size = guide_legend(title.position = "top")) +
                                           theme(axis.title = element_blank(),
                                                 axis.ticks = element_blank(),
                                                 axis.text = element_blank(),
                                                 legend.position = "top",
                                                 legend.direction = "horizontal",
                                                 legend.background = element_blank(), legend.key = element_blank()))

#add shape legend to plot
full_nmds_plot_byregion_dist_200m_legend <- ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  draw_plot(full_nmds_plot_byregion_dist_200m, x = 0, y = 0, height = 1, width = 1) +
  draw_plot(dist_200m_legend, x = 0.17, y = 0.7, height = 0.2, width = 0.2)

ggsave(full_nmds_plot_byregion_dist_200m_legend, path = "figures", filename = "full_nmds_plot_byregion_dist_200m_legend.jpg", height = 7, width = 7)


#Visualize by temp
#link kelp density
dat_averages_bysite.nmds_200m_color <- ggplot(data = dat_averages_bysite.nmds_dist_200m) +
  geom_point(aes(x = MDS1, y = MDS2, color = dist_200m_bath)) +
  # scale_size_continuous(guide = "none") +
  scale_color_viridis() +
  # lims(x= c(-1.4,1.2), y = c(-2,2)) +
  labs(color = "Meters to 200m Isobath") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        #   legend.title = element_blank(),
        legend.background = element_blank(), legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=4)))

full_nmds_plot_byregion_meantemp 



#Visualize by temp
#link kelp density
dat_averages_bysite.nmds_meantemp <- VRG_lat_lon_only[dat_averages_bysite.nmds, on = c("Site","DepthZone")]

full_nmds_plot_byregion_meantemp <- ggplot(data = dat_averages_bysite.nmds_meantemp) +
  geom_point(aes(x = MDS1, y = MDS2, color = BO_sstmean)) +
 # scale_size_continuous(guide = "none") +
  scale_colour_gradientn(colors = c("blue","lightblue","pink","red")) +
  # lims(x= c(-1.4,1.2), y = c(-2,2)) +
  labs(color = "Temperature") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
     #   legend.title = element_blank(),
        legend.background = element_blank(), legend.key = element_blank()) +
  guides(color = guide_legend(override.aes = list(size=4)))

full_nmds_plot_byregion_meantemp 


#extract legend for shape only
dist_200m_legend <- get_legend(ggplot(data = dat_averages_bysite.nmds_meantemp) +
                                 geom_point(aes(x = MDS1, y = MDS2,size = dist_200m_bath/1000)) +
                                 scale_size_continuous("Kilometers to 200m isobath") +
                                 # lims(x= c(-1.4,1.2), y = c(-2,2)) +
                                 theme_classic() +
                                 guides(size = guide_legend(title.position = "top")) +
                                 theme(axis.title = element_blank(),
                                       axis.ticks = element_blank(),
                                       axis.text = element_blank(),
                                       legend.position = "top",
                                       legend.direction = "horizontal",
                                       legend.background = element_blank(), legend.key = element_blank()))

#add shape legend to plot
full_nmds_plot_byregion_meantemp_legend <- ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  draw_plot(full_nmds_plot_byregion_meantemp, x = 0, y = 0, height = 1, width = 1) +
  draw_plot(dist_200m_legend, x = 0.17, y = 0.7, height = 0.2, width = 0.2)

ggsave(full_nmds_plot_byregion_meantemp_legend, path = "figures", filename = "full_nmds_plot_byregion_meantemp_legend.jpg", height = 7, width = 7)







#PermANOVA for FISH with kelp as predictor
dat_fish_averages_bysite.wide.rt_pluskelp <- macro_density_bysite[dat_fish_averages_bysite.wide.rt, on = c("Site","DepthZone")]


permanova_fish_kelp_predictor <- adonis2(
  dat_fish_averages_bysite.wide.rt_pluskelp[,5:ncol(dat_fish_averages_bysite.wide.rt_pluskelp)] ~ dat_fish_averages_bysite.wide.rt_pluskelp$macro_mean_density_m2 + dat_fish_averages_bysite.wide.rt_pluskelp$DepthZone,
  method = "bray"
)

permanova_fish_kelp_predictor

#Macrocystis density explains 6% of fish community structure, depth zone explains 19% of community structure

















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





