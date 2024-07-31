# CREATION DATE 28 Jan 2024
# MODIFICATION DATE 29 July 2024

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Community visualizations and analyses

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
dat_fish_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]
dat_macroinvert_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]
dat_kelp_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]


##################################################
#Long data to wide data for vegan analyses
##################################################
#environmental predictors
#pull environmental data
all_env_lat_lon <- fread(file.path("data","enviro_predictors","all_env_lat_lon.csv"))

keep <- colnames(all_env_lat_lon)
keep <- keep[!(keep %in% c("Region"))]

#remove region, and lat lon (FOR SOME REASON EXTRA ROWS FOR UNIQUE DEPTHZONE and SITE COMBOS, NEED TO FIX in pull_enviro_data once CDFW bathymetry link is back up)
all_env_lat_lon <- all_env_lat_lon[,..keep]

nrow(unique(all_env_lat_lon[,.(DepthZone, Site)]))

nrow(unique(all_env_lat_lon[,.(DepthZone, Site, Substrate_SD, Relief_SD, mean_sst_C, dist_200m_bath)]))

all_env_lat_lon <- all_env_lat_lon[, lapply(.SD, mean, na.rm=TRUE), by=c("Site","DepthZone")] 
#temporary fix above

#fish
dat_fish_long_density <- dat_fish_site_averages[,.(taxa,Region,Site,DepthZone, mean_density_m2)]
dat_fish_long_density[,mean_density_m2 := mean(mean_density_m2),.(taxa,Region,Site,DepthZone)]
dat_fish_long_density <- unique(dat_fish_long_density[,.(taxa,Region,Site,DepthZone,mean_density_m2)])
#merge with environmental
dat_fish_long_density <- dat_fish_long_density[all_env_lat_lon, on = c("Site","DepthZone")] #dropped rows, not sure why

dat_fish_long_biomass <- dat_fish_site_averages[,.(taxa,Region,Site,DepthZone, mean_wt_density_g_m2)]
dat_fish_long_biomass[,mean_wt_density_g_m2 := mean(mean_wt_density_g_m2),.(taxa,Region,Site,DepthZone)]
dat_fish_long_biomass <- unique(dat_fish_long_biomass[,.(taxa,Region,Site,DepthZone,mean_wt_density_g_m2)])
#merge with environmental
dat_fish_long_biomass <- dat_fish_long_biomass[all_env_lat_lon, on = c("Site","DepthZone")] #dropped rows, not sure why

#macroinvert
dat_macroinvert_long_density <- dat_macroinvert_site_averages[,.(taxa,Region,Site,DepthZone, mean_density_m2)]
dat_macroinvert_long_density[,mean_density_m2 := mean(mean_density_m2),.(taxa,Region,Site,DepthZone)]
dat_macroinvert_long_density <- unique(dat_macroinvert_long_density[,.(taxa,Region,Site,DepthZone,mean_density_m2)])
#merge with environmental
dat_macroinvert_long_density <- dat_macroinvert_long_density[all_env_lat_lon, on = c("Site","DepthZone")] #dropped rows, not sure why



#kelp
dat_kelp_long_density <- dat_kelp_site_averages[,.(taxa,Region,Site,DepthZone, mean_density_m2)]
dat_kelp_long_density[,mean_density_m2 := mean(mean_density_m2),.(taxa,Region,Site,DepthZone)]
dat_kelp_long_density <- unique(dat_kelp_long_density[,.(taxa,Region,Site,DepthZone,mean_density_m2)])
#merge with environmental
dat_kelp_long_density <- dat_kelp_long_density[all_env_lat_lon, on = c("Site","DepthZone")] #dropped rows, not sure why, check later



#melt long to wide
dat_fish_wide_density <- dcast(dat_fish_long_density, Region + Site + DepthZone + mean_chl_mg_m3 + max_chl_mg_m3 + min_chl_mg_m3 + mean_sst_C + max_sst_C + min_sst_C + dist_200m_bath +  giantkelp_stipe_density_m2 + Relief_index + Relief_SD + Substrate_index + Substrate_SD + Latitude ~ taxa, value.var = "mean_density_m2", fun = mean)

dat_fish_wide_biomass <- dcast(dat_fish_long_biomass, Region + Site + DepthZone + mean_chl_mg_m3 + max_chl_mg_m3 + min_chl_mg_m3 + mean_sst_C + max_sst_C + min_sst_C + dist_200m_bath +  giantkelp_stipe_density_m2 + Relief_index + Relief_SD + Substrate_index + Substrate_SD + Latitude ~ taxa, value.var = "mean_wt_density_g_m2", fun = mean)

dat_macroinvert_wide_density <- dcast(dat_macroinvert_long_density, Region + Site + DepthZone + mean_chl_mg_m3 + max_chl_mg_m3 + min_chl_mg_m3 + mean_sst_C + max_sst_C + min_sst_C + dist_200m_bath +  giantkelp_stipe_density_m2 + Relief_index + Relief_SD + Substrate_index + Substrate_SD + Latitude ~ taxa, value.var = "mean_density_m2", fun = mean)

dat_kelp_wide_density <- dcast(dat_kelp_long_density, Region + Site + DepthZone + mean_chl_mg_m3 + max_chl_mg_m3 + min_chl_mg_m3 + mean_sst_C + max_sst_C + min_sst_C + dist_200m_bath +  giantkelp_stipe_density_m2 + Relief_index + Relief_SD + Substrate_index + Substrate_SD + Latitude ~ taxa, value.var = "mean_density_m2", fun = mean)

#spp per depthzone
dat_fish_long_density_removezeros <- dat_fish_long_density[mean_density_m2>0,]
dat_fish_long_biomass_removezeros <- dat_fish_long_biomass[mean_wt_density_g_m2>0,]
dat_macroinvert_long_density_removezeros <- dat_macroinvert_long_density[mean_density_m2>0,]
dat_kelp_long_density_removezeros <- dat_kelp_long_density[mean_density_m2>0,]


#number of taxa per depth zone?
fish_depthzone_site_richness_abun <- dat_fish_long_density_removezeros[,.(count_spp = .N,sum_abun = sum(mean_density_m2)),.(Site, DepthZone, Region)]
fish_depthzone_site_richness_abun[,category:="fish"]
fish_depthzone_site_biomass <- dat_fish_long_biomass_removezeros[,.(sum_biomass = sum(mean_wt_density_g_m2)),.(Site, DepthZone, Region)]
fish_depthzone_site_richness_abun_biomass <- fish_depthzone_site_richness_abun[fish_depthzone_site_biomass, on = c("Site","DepthZone","Region")]

macroinvert_depthzone_site_richness <- dat_macroinvert_long_density_removezeros[,.(count_spp = .N, sum_abun = sum(mean_density_m2)),.(Site, DepthZone, Region)]
macroinvert_depthzone_site_richness[,category:="macroinvert"]
kelp_depthzone_site_richness <- dat_kelp_long_density_removezeros[,.(count_spp = .N, sum_abun = sum(mean_density_m2)),.(Site, DepthZone, Region)]
kelp_depthzone_site_richness[,category:="kelp"]

#rbind these
depthzone_site_richness <- rbind(fish_depthzone_site_richness_abun_biomass, macroinvert_depthzone_site_richness, kelp_depthzone_site_richness, fill = TRUE)

##################################################################################################
#Long data to wide data for multivariate community visualizations and analyses
##################################################################################################

#melt long to wide
dat_fish_averages_bysite.wide <- dcast(dat_fish_long_density_removezeros, Region + Site + DepthZone + mean_chl_mg_m3 + max_chl_mg_m3 + min_chl_mg_m3 + mean_sst_C + max_sst_C + min_sst_C + dist_200m_bath +  giantkelp_stipe_density_m2 + Relief_index + Relief_SD + Substrate_index + Substrate_SD + Latitude ~ taxa, value.var = "mean_density_m2", fun = mean, fill = 0)

dat_fish_biomass_averages_bysite.wide <- dcast(dat_fish_long_biomass_removezeros, Region + Site + DepthZone ~ taxa, value.var = "mean_wt_density_g_m2", fun = mean, fill = 0)

dat_macroinvert_averages_bysite.wide <- dcast(dat_macroinvert_long_density_removezeros, Region + Site + DepthZone ~ taxa, value.var = "mean_density_m2", fun = mean, fill = 0)

dat_kelp_averages_bysite.wide <- dcast(dat_kelp_long_density_removezeros, Region + Site + DepthZone ~ taxa, value.var = "mean_density_m2", fun = mean, fill = 0)

#Merge for all density data tables for species visualization (biomass of fish kept separate)
dat_averages_bysite.wide <- dat_fish_averages_bysite.wide[dat_macroinvert_averages_bysite.wide, on = c("Region","Site","DepthZone")]
dat_averages_bysite.wide <- dat_kelp_averages_bysite.wide[dat_averages_bysite.wide, on = c("Region","Site","DepthZone")]

#there are some sites that have no kelp, and these come up as NAs, change to 0 
dat_averages_bysite.wide[is.na(dat_averages_bysite.wide)] <- 0


##################################################################################################
#Visualize full PcoA (all species, colored by Depth Zone, no ARM yet)
##################################################################################################
#CHECK FOR NUMBERS IN 165
stopifnot(colnames(dat_averages_bysite.wide[,c(1:3,24,36)]) == c("Region","Site","DepthZone","mean_chl_mg_m3","Latitude"))
          
#species only
dat_averages_bysite.wide.spp <- dat_averages_bysite.wide[,c(4:23,37:210)]

#log base 2 transform (Jonathan F. Jupke, Ralf B. Schäfer)
dat_averages_bysite.wide.spp.l <- log2(1+dat_averages_bysite.wide.spp)

#distance matrix (distance = bray curtis)
dat_averages_bysite.dist <- vegdist(dat_averages_bysite.wide.spp.l, distance = "bray")

#PCoA
allspp_PCoA <- wcmdscale(dat_averages_bysite.dist, eig = TRUE)

allspp_PCoA.pnts <- data.table(allspp_PCoA$points)

#add env and group variables
PCoA_logtrans_env <- cbind(dat_averages_bysite.wide[,c(1:3,24:36)],allspp_PCoA.pnts[,1:2])

#adjust factor order
PCoA_logtrans_env[,DepthZone := factor(DepthZone,
                                       levels = c("Inner","Middle","Outer","Deep","ARM"),
                                       labels = c("Inner","Middle","Outer","Deep","AR"))]

#add mainland versus island designation
PCoA_logtrans_env[,type := factor(ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")))]


#visualize
PCoA_allspp_allsite <- ggplot(PCoA_logtrans_env) +
  stat_ellipse(geom = "polygon", aes(Dim1, Dim2, fill = DepthZone), color = NA,alpha = 0.4) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100","black"), guide = guide_legend(reverse = TRUE)) +
  geom_point(aes(Dim1, Dim2, color = DepthZone, shape = DepthZone), size = 4, fill = "#FE6100") +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100","black"), guide = guide_legend(reverse = TRUE)) +
 # lims(x = c(-0.5,0.45),y = c(-0.5,0.35)) +
  scale_shape_manual(values = c(15,17,19,23,7), guide = guide_legend(reverse = TRUE)) +
  theme_classic()

#central point for each
PCoA_logtrans_env.centroid <- PCoA_logtrans_env[, lapply(.SD, mean, na.rm=TRUE), by=c("DepthZone"), .SDcols = c("Dim1","Dim2")] 

PCoA_allspp_centroids <- ggplot() +
  stat_ellipse(data = PCoA_logtrans_env, geom = "polygon", aes(Dim1, Dim2, fill = DepthZone), color = "white",alpha = 0.2) +
  scale_fill_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100","black"), guide = guide_legend(reverse = TRUE)) +
  geom_point(data = PCoA_logtrans_env.centroid, aes(Dim1, Dim2, color = DepthZone, shape = DepthZone), size = 6, fill = "#FE6100") +
  scale_color_manual(values = c("#015AB5", "#785EF0","#DC277F","#FE6100","black"), guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(values = c(15,17,19,23,7), guide = guide_legend(reverse = TRUE)) +
  labs(x = "PCoA 1",y = "PCoA 2")+
  theme_classic() +
  theme(legend.position = "top", legend.justification = "center", legend.title = element_blank(),
        legend.text = element_text(size = 18, face = "bold"),   # Increase legend text size
        legend.key.size = unit(1.5, "lines"))     # Increase legend key size

#plot each depth zone independently, as this is how we will do variable analysis
#DEEP
PCoA_allspp_deep <- ggplot(PCoA_logtrans_env[DepthZone == "Deep"]) +
  geom_point(aes(Dim1, Dim2, shape = type, fill = type),color = "darkgrey", size = 4) +
  scale_fill_manual(values = c("#F3C393","#FE6100")) +
  scale_shape_manual(values = c(24,25)) +
  lims(y = c(-0.6,0.4),x = c(-0.65,0.35)) +
  labs(fill = "Reef type", shape = "Reef type", x = "PCoA 1",y = "PCoA 2")+
  theme_classic()

#Outer
PCoA_allspp_outer <- ggplot(PCoA_logtrans_env[DepthZone == "Outer"]) +
  geom_point(aes(Dim1, Dim2, shape = type, fill = type),color = "darkgrey", size = 4) +
  scale_fill_manual(values = c("#E0D1F1","#785EF0")) +
  scale_shape_manual(values = c(24,25)) +
  lims(y = c(-0.65,0.35),x = c(-0.6,0.4)) +
  labs(fill = "Reef type", shape = "Reef type", x = "PCoA 1",y = "PCoA 2")+
  theme_classic()

#Middle
PCoA_allspp_middle <- ggplot(PCoA_logtrans_env[DepthZone == "Middle"]) +
  geom_point(aes(Dim1, Dim2, shape = type, fill = type),color = "darkgrey", size = 4) +
  scale_fill_manual(values = c("#F3E0F3","#DC277F")) +
  scale_shape_manual(values = c(24,25)) +
  lims(y = c(-0.6,0.4),x = c(-0.5,0.5)) +
  labs(fill = "Reef type", shape = "Reef type", x = "PCoA 1",y = "PCoA 2")+
  theme_classic()

#Inner
PCoA_allspp_inner <- ggplot(PCoA_logtrans_env[DepthZone == "Inner"]) +
  geom_point(aes(Dim1, Dim2, shape = type, fill = type),color = "darkgrey", size = 4) +
  scale_fill_manual(values = c("#ADD4F5","#015AB5")) +
  scale_shape_manual(values = c(24,25)) +
  lims(y = c(-0.7,0.3),x = c(-0.4,0.6)) +
  labs(fill = "Reef type", shape = "Reef type", x = "PCoA 1",y = "PCoA 2")+
  theme_classic()

#Artificial reef
PCoA_allspp_AR <- ggplot(PCoA_logtrans_env[DepthZone == "AR"]) +
  geom_point(aes(Dim1, Dim2),color = "darkgrey", size = 4, shape = 25, fill = "black") +
  scale_shape_manual(values = c(24,25)) +
  lims(y = c(-0.45,0.55),x = c(-0.7,0.3)) +
  labs(fill = "Reef type", shape = "Reef type", x = "PCoA 1",y = "PCoA 2")+
  theme_classic()

#DUMMY FOR LEGEND
PCoA_allspp_legend <- get_legend(ggplot(PCoA_logtrans_env[DepthZone == "Inner"]) +
  geom_point(aes(Dim1, Dim2, shape = type, fill = type),color = "darkgrey", size = 8) +
  scale_fill_manual(values = c("lightgrey","black")) +
  scale_shape_manual(values = c(24,25)) +
  labs(fill = "Reef type", shape = "Reef type", x = "PCoA 1",y = "PCoA 2")+
  theme_classic() +
    theme(  legend.title = element_text(size = 20),  # Increase legend title size
            legend.text = element_text(size = 16),   # Increase legend text size
            legend.key.size = unit(1.5, "lines")     # Increase legend key size
    ))

#merge all plots
PCoA_allspp_zones <- plot_grid(PCoA_allspp_inner+theme(legend.position = "null"),
          PCoA_allspp_outer+theme(legend.position = "null"),
          PCoA_allspp_middle+theme(legend.position = "null"),
          PCoA_allspp_deep+theme(legend.position = "null"),
          PCoA_allspp_AR + theme(legend.position = "null"),
          PCoA_allspp_legend,  ncol = 2, labels = c(" Inner","Middle"," Outer"," Deep","   AR"), label_x = 0.09)

PCoA_allspp_full <- plot_grid(PCoA_allspp_centroids, PCoA_allspp_zones, ncol = 1, rel_heights = c(0.8,1))

ggsave(PCoA_allspp_full, path = "figures", filename = "PCoA_allspp_full.jpg",height = 11, width = 7, unit = "in" )

##################################################################################################
#Distance-based ReDundancy Analysis (dbRDA)
##################################################################################################

#base 2 log-transformed data

#first step of a dbRDA is to apply a PCoA to the distance matrix and to keep all principal coordinates

#












#SKIP FOR NOW
########################
#Latent variable modeling
#########################
#response = counts (poisson or negative binomial; log link; method = VA/LA)
#response = biomass (tweedie distribution, log link; method = LA)
#row.eff = defining type of row effects (none, fixed, random)
#offset = offsets
#power = defining power parameter of tweedie distribution
#starting.val = judicious chocie of starting values for latent variables

#GLLVMs can be used as a model-based approach to unconstrained ordination by including two latent variables in the model but NO predictors
  #corresponding ordination plot then provides a graphical representation of which sites are similar in terms of thei rspecies composition


#model based ordination
y <- as.matrix(site_spp)
x <- as.matrix(site_env.f)

fitp <- gllvm(y, family = poisson()) #takes a few minutes
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1)) #set up full plot space
plot(fitp, var.colors = 1) #visualize model performance
saveRDS(fitp, file.path("model_output","fitp.rds"))
#fitp <- readRDS(file.path("model_output","fitp.rds"))
#Warning message:
#  In gllvm(y, family = poisson()) : There are responses full of zeros. (not sure what that means)

fitvm <- gllvm(y, family = "negative.binomial")
par(mfrow = c(3, 2), mar = c(4, 4, 2, 1)) #set up full plot space
plot(fitvm, var.colors = 1) #visualize model performance
saveRDS(fitvm, file.path("model_output","fitvm.rds"))


#check AIC to see which is a better model
#poisson is better model!

#note that latent variables provide some capacity to account for overdispersion, so overdispersed counts do not always require us to move beyond Poisson

#visualize!
ordiplot(fitp, biplot = T, ind.spp = 20, xlim = c(-3,4.5), ylim = c(-0.0005,0.0005)) #for some reason the ordination can't be visualized right

#Including environmental variables

#example
library(mvabund)
data(antTraits)


#####################
#PERMANOVA, Permutational Multivariate Analysis of Variance (perMANOVA)
#####################
#FULL COMMUNITY

#Square root transformation (high abundance spp less pull in ordination)
dat_averages_bysite.wide.rt <- sqrt(dat_averages_bysite.wide[,c(4:ncol(dat_averages_bysite.wide)), with = FALSE])


#add back reference columns
dat_averages_bysite.wide.rt <- cbind(dat_averages_bysite.wide[,c(1:3), with = FALSE], dat_averages_bysite.wide.rt)

    #alternative version without artificial reefs
    dat_averages_bysite_NOARM.wide.rt <- dat_averages_bysite.wide.rt[DepthZone != "ARM",]

#only keep Species
dat_averages_bysite.wide.rt.trim <- dat_averages_bysite.wide.rt[,c(4:ncol(dat_averages_bysite.wide.rt)), with = FALSE]

  #alternative version without artificial reefs
  dat_averages_bysite_NOARM.wide.rt.trim <- dat_averages_bysite_NOARM.wide.rt[,c(4:ncol(dat_averages_bysite_NOARM.wide.rt)), with = FALSE]

#perMANOVA including ARM as depth zone
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

#depth zone accounts for 19% of variation for natural reef sites
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
  dat_kelp_averages_bysite.wide.rt.trim[,1:20] ~ dat_kelp_averages_bysite.wide.rt$DepthZone,
  method = "bray"
)

permanova_kelp

#kelp names
kelp_names <- colnames(dat_kelp_averages_bysite.wide.rt.trim[,1:20])

#whether or not a site is an artificial reef affects fish community composition, accounting for 21% of variation in composition

###################
#PLOT NMDS
###################
#NMDS “preserves the rank order of the inter-point dissimilarities as well as possible within the constraints of a small number of dimensions” (Anderson et al. 2008, p. 105).
#It has the least restrictive assumptions (any distance measure, no assumptions about shape of species response curves).

#We may replace these with latent variable plots

#all species

#convert to distance matrix
dat_averages_bysite.dist <- vegdist(dat_averages_bysite.wide[,4:ncol(dat_averages_bysite.wide)], distance = "bray")

full_nmds <- metaMDS(dat_averages_bysite.wide[,4:ncol(dat_averages_bysite.wide)], #metaMDS automatically does this for us
                     autotransform = FALSE,
                     distance = "bray",
                     engine = "monoMDS",
                     wascores = T, #weighted average scores
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
    #SEE LINES BELOW
    #output: species_scores_top

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
  ggrepel::geom_text_repel(data = species_scores_top, aes(x = MDS1, y = MDS2, label = Species, color = `Spp category`), size = 2.5, fontface ="bold.italic") +
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
  ggrepel::geom_text_repel(data = species_scores_top, aes(x = MDS1, y = MDS2, label = Species), size = 2.5, fontface ="bold.italic") +
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
  geom_point(aes(x = MDS1, y = MDS2, color = Region), size = 2.5) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "black", "#A65628" ,"#F781BF")) +
  new_scale_color() +
  stat_ellipse(data = dat_averages_bysite.nmds, aes(x = MDS1, y = MDS2, color = DepthZone)) + #95% confidence level for a mlutivariate t-distribution
  scale_color_manual(values = c("cornflowerblue","coral1","chartreuse3","darkorchid2","black")) +
 # guides(color = guide_legend(nrow = 2)) +
   # lims(x= c(-1.4,1.2), y = c(-2,2)) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
    #    legend.position = "top",
      #  legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.background = element_blank(), legend.key = element_blank())


ggsave(full_nmds_plot_byregion, path = "figures", filename = "full_nmds_plot_byregion.jpg", height = 7, width = 9)

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

#extract legend for color ellipse only
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

#######NEED TO ADD IN HABITAT/ENVIRONMENTAL VARIABLES AGAIN
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
















#######################
########SPECIES SCORES
#######################
#For species scores, you can do it in a few different ways. See this page for a helpful explanation of the differences: https://stackoverflow.com/questions/60937869/using-envfit-vegan-to-calculate-species-scores
# "wascores() are best thought of as optima, where the abundance declines as you move away from the weighted centroid"
# "wascores() takes the coordinates of the points in the nmds space and computes the mean on each dimension, weighting observations by the abundance of the species at each point. Hence the species score returned by wascores() is a weighted centroid in the NMDS space for each species, where the weights are the abundances of the species."
# "envfit() in an NMDS is not necessarily a good thing as we wouldn't expect abundances to vary linearly over the ordination space."


#I believe we should be calculating species scores from raw data, not from transformed data as are automatically in the NMDS

species_scores <- data.table(vegan::wascores(full_nmds$points, #ordination scores
                                             dat_averages_bysite.wide[,4:ncol(dat_averages_bysite.wide)]#weights; species abundances
                                             ))

spp_list <- colnames(dat_averages_bysite.wide[,4:ncol(dat_averages_bysite.wide)])
species_scores[,Species := spp_list]

#add category
species_scores[,`Spp category`:= c(rep("Fish",length(fish_names)), rep("Kelp",length(kelp_names)), rep("Macroinvert",length(macro_names)))]
#add category
species_scores[,`Spp category`:= c(rep("Fish",length(fish_names)), rep("Kelp",length(kelp_names)), rep("Macroinvert",length(macro_names)))]


#top and bottom three for NMDS1 & 2
#sort by MDS1
setkey(species_scores, MDS1)
mds1_full_spp_scores1 <- head(species_scores[complete.cases(species_scores),],5)[,Species]
mds1_full_spp_scores2 <- tail(species_scores[complete.cases(species_scores),],5)[,Species]
setkey(species_scores, MDS2)
mds1_full_spp_scores3 <- head(species_scores[complete.cases(species_scores),],5)[,Species]
mds1_full_spp_scores4 <- tail(species_scores[complete.cases(species_scores),],5)[,Species]

top_all_spp <- sort(unique(c(mds1_full_spp_scores1,
                 mds1_full_spp_scores2,
                 mds1_full_spp_scores3,
                 mds1_full_spp_scores4
                 )))

species_scores_top <- species_scores[Species %in% top_all_spp]

#add to simple NMDS plot

full_nmds_plot_speciesscores <- ggplot(data = dat_averages_bysite.nmds) +
  geom_point(aes(x = MDS1, y = MDS2,shape = `Site type`), size = 3, alpha = 0.3) +
  geom_text_repel(data = species_scores_top, aes(x = MDS1, y = MDS2, label = Species), size = 2) +
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

#START HERE
#GENERALIZED LINEAR LATENT VARIABLE MODELS (glvvm)
library(gllvm)
col_keep <- c(4:23,37:ncol(dat_averages_bysite.wide))

dat_averages_bysite.wide.bio <- dat_averages_bysite.wide[,..col_keep][complete.cases(dat_averages_bysite.wide)]
dat_averages_bysite.wide.envnum <- dat_averages_bysite.wide[,24:36][complete.cases(dat_averages_bysite.wide)]
dat_averages_bysite.wide.envcat <- dat_averages_bysite.wide[,1:3][complete.cases(dat_averages_bysite.wide)]

dat_averages_bysite.wide.envnum <- scale(dat_averages_bysite.wide.envnum) #centered and scaled

#fit model
fitp <- gllvm(dat_averages_bysite.wide.bio, family = poisson())
fit_ord <- gllvm(dat_averages_bysite.wide.bio, family = poisson())

#equally good fit for both, stick with Poisson
ordiplot(fitp, biplot = T, ind.spp = 15, 
         xlim = c(-5,3), ylim = c(-0.0000012,0.0000020))
points(fitp, biplot = T, col = dat_averages_bysite.wide.envcat$DepthZone)

rbPal <- c("#83752D", "#EAC4AA", "#D38EB6", "#F6D798","#76A6B8", "#2F624E")

for(i in 1:ncol(X)){
  Col <- rbPal[as.numeric(cut(X[,i], breaks = 20))]
  ordiplot(fitnb, symbols = T, s.colors = Col, main = colnames(X)[i], 
           biplot = TRUE)

cr <- getResidualCor(fitp)

ggsave(cr, path = file.path("figures"), file = "correlation.pdf", width = 20, height = 20, unit = "in")

library(corrplot); library(gclus)

pdf(file = file.path("figures","correlation_gllvm.pdf"))
corrplot(cr, diag = FALSE, type = "lower", method = "square", tl.srt = 25, tl.cex = 0.2, addgrid.col = NA)
dev.off()


data("antTraits")
antTraits[1]
#can take long data
gllvm()

ordiplot(fitp, biplot = TRUE)
abline(h = 0, v = 0, lty=2)

#assess model fit
par(mfrow = c(1,2))
plot(fitp, which = 1:2)

#
