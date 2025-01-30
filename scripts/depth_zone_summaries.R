# CREATION DATE 28 Jan 2024
# MODIFIED DATE 20 Dec 2024

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Create depth zone summaries (feed into first figure in paper)

#############################
##Setup
#############################
library(ggplot2)
library(data.table)
library(cowplot)
library(ggpattern) #patterned bars
library(RColorBrewer)

########################
##Load data
########################
dat_event.r <- readRDS(file.path("data","processed_crane", "dat_event.r.rds"))
dat_fish_site_averages <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages.rds"))
dat_macroinvert_site_averages <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages.rds"))
dat_kelp_site_averages <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages.rds"))

#Load hex codes for all top spp
spp_color_key <- fread(file.path("keys","spp_color_key.csv"), strip.white = T)

#included giant kelp stipes for in situ habitat data, but DELETE for community analyses
dat_kelp_site_averages <- dat_kelp_site_averages[BenthicReefSpecies != "Macrocystis pyrifera stipes",] 

#pull in spp taxonomy info
species_key <- fread(file.path("keys","species_key.csv"))
#capitalize California
species_key[, common_name_final := gsub("california sheephead", "California sheephead", common_name_final)]
#add common name for Chestnut whelk
species_key[, common_name_final := ifelse(taxa == "Neobernaya spadicea","chestnut cowry",common_name_final)]
#add common name for elk kelp
species_key[, common_name_final := ifelse(taxa == "Pelagophycus porra","elk kelp",common_name_final)]


#find and replace Semicossyphus pulcher with Bodianus pulcher in VRG data
dat_fish_site_averages[, Species := gsub("Semicossyphus pulcher", "Bodianus pulcher", Species)]

#find and replace Hermosilla azurea  with Kyphosus azureus in VRG data
dat_fish_site_averages[, Species := gsub("Hermosilla azurea", "Kyphosus azureus", Species)]

#change genus for Neoagarum fimbriatum in VRG data
dat_kelp_site_averages[, BenthicReefSpecies := gsub("Agarum fimbriatum", "Neoagarum fimbriatum", BenthicReefSpecies)]

#change genus for Neobernaya spadicea (seems like a spelling error?) in VRG data
dat_macroinvert_site_averages[, BenthicReefSpecies := gsub("Neobemaya spadicea", "Neobernaya spadicea", BenthicReefSpecies)]

#link site averaged data with species key
dat_fish_site_averages <- species_key[dat_fish_site_averages, on = c("taxa" = "Species")]
dat_macroinvert_site_averages <- species_key[dat_macroinvert_site_averages, on = c("taxa" = "BenthicReefSpecies")]
dat_kelp_site_averages <- species_key[dat_kelp_site_averages, on = c("taxa" = "BenthicReefSpecies")]

#Manually add common name for chainbladder kelp
dat_kelp_site_averages[,common_name_final := ifelse(taxa == "Stephanocystis spp.", "chain bladder kelp",common_name_final)]

#Manually add common name for sargassum sp
dat_kelp_site_averages[,common_name_final := ifelse(taxa == "Sargassum sp", "",common_name_final)]

#Split depth zone for artificial reefs in PV and SM
dat_fish_site_averages[,DepthZone := ifelse(DepthZone %in% c("ARM","Module") & AR_Complex == "Palos Verdes Reef","AR_PVR",ifelse(DepthZone %in% c("ARM","Module") & AR_Complex != "Palos Verdes", "AR_SM",as.character(DepthZone)))]
dat_kelp_site_averages[,DepthZone := ifelse(DepthZone %in% c("ARM","Module") & AR_Complex == "Palos Verdes Reef","AR_PVR",ifelse(DepthZone %in% c("ARM","Module") & AR_Complex != "Palos Verdes", "AR_SM",as.character(DepthZone)))]
dat_macroinvert_site_averages[,DepthZone := ifelse(DepthZone %in% c("ARM","Module") & AR_Complex == "Palos Verdes Reef","AR_PVR",ifelse(DepthZone %in% c("ARM","Module") & AR_Complex != "Palos Verdes", "AR_SM",as.character(DepthZone)))]

#Summary of avg site counts and biomass ####
dat_fish_site_averages.summary <- dat_fish_site_averages[,.(density_count100m2_summed = round(sum(mean_density_m2*100),1),biomass_kg100m2_summed = round(sum(mean_wt_density_g_m2)/1000*100,1)),.(DepthZone,Site)]
dat_kelp_site_averages.summary <- dat_kelp_site_averages[,.(density_count100m2_summed = round(sum(mean_density_m2)*100,1)), .(DepthZone,Site)]
dat_macroinvert_site_averages.summary <- dat_macroinvert_site_averages[,.(density_count100m2_summed = round(sum(mean_density_m2)*100,1)),.(DepthZone,Site)]
########################
##Averaged across all sites, top species per depth zone
########################

#fish
dat_fish_site_averages[,mean_depthzone_density_m2 := mean(mean_density_m2),.(taxa, DepthZone)] 
dat_fish_site_averages[,mean_depthzone_wt_density_g_m2 := mean(mean_wt_density_g_m2),.(taxa, DepthZone)] 

col_keep_fish <- colnames(dat_fish_site_averages[,c(1:10,14,17,18)])

dat_fish_averages <- unique(dat_fish_site_averages[,..col_keep_fish])

#invert
dat_macroinvert_site_averages[,mean_depthzone_density_m2 := mean(mean_density_m2),.(taxa, DepthZone)] 

col_keep_macroinvert <- colnames(dat_macroinvert_site_averages[,c(1:10,15,17)])

dat_macroinvert_averages <- unique(dat_macroinvert_site_averages[,..col_keep_macroinvert])
  
#kelp
dat_kelp_site_averages[,mean_depthzone_density_m2 := mean(mean_density_m2),.(taxa, DepthZone)] 

col_keep_kelp <- colnames(dat_kelp_site_averages[,c(1:10,15,17)])

dat_kelp_averages <- unique(dat_kelp_site_averages[,..col_keep_kelp])

#Add column for whether species is in top 5 species (by density) for that zone
#fish density
dat_fish_averages[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"), .(DepthZone)]
dat_fish_averages[, common_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,""), .(DepthZone)]
dat_fish_averages[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
dat_fishdensity_averages.u <- unique(dat_fish_averages[,.(Species_top5,common_top5, DepthZone, summed_mean_depthzone_density_m2)])
dat_fishdensity_averages.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]
dat_fishdensity_averages.u[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]
dat_fishdensity_averages.u[,rel_abun := round(summed_mean_depthzone_density_m2/sum(summed_mean_depthzone_density_m2),2),.(DepthZone)]

#macroinvert density
dat_macroinvert_averages[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"), .(DepthZone)]
dat_macroinvert_averages[, common_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,""), .(DepthZone)]
dat_macroinvert_averages[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
dat_macroinvertdensity_averages.u <- unique(dat_macroinvert_averages[,.(Species_top5,common_top5, DepthZone, summed_mean_depthzone_density_m2)])
dat_macroinvertdensity_averages.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]
dat_macroinvertdensity_averages.u[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]
dat_macroinvertdensity_averages.u[,rel_abun := round(summed_mean_depthzone_density_m2/sum(summed_mean_depthzone_density_m2),2),.(DepthZone)]


#Kelp density
dat_kelp_averages[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"), .(DepthZone)]
dat_kelp_averages[, common_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,""), .(DepthZone)]
dat_kelp_averages[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
dat_kelp_averages.u <- unique(dat_kelp_averages[,.(Species_top5, common_top5, DepthZone, summed_mean_depthzone_density_m2)])
dat_kelp_averages.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]
dat_kelp_averages.u[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]
dat_kelp_averages.u[,rel_abun := round(summed_mean_depthzone_density_m2/sum(summed_mean_depthzone_density_m2),2),.(DepthZone)]


#Add column for whether species is in top 5 species (by biomass) for that zone, and then sum biomass for all others
#fish biomass
dat_fish_averages[, Species_top5_biomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,taxa,"Other"), .(DepthZone)]
dat_fish_averages[, common_top5_biomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,common_name_final,""), .(DepthZone)]
dat_fish_averages[, summed_mean_depthzone_biomass_m2 := sum(mean_depthzone_wt_density_g_m2), .(DepthZone, Species_top5_biomass)]
dat_fishbiomass_averages.u <- unique(dat_fish_averages[,.(Species_top5_biomass,common_top5_biomass, DepthZone, summed_mean_depthzone_biomass_m2)])
dat_fishbiomass_averages.u[,full_label := ifelse(Species_top5_biomass == "Other","Other",paste0(Species_top5_biomass,"\n", common_top5_biomass))]
dat_fishbiomass_averages.u[,DepthZone := factor(DepthZone, levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"))]


#manually set factor order for plotting
dat_fishbiomass_averages.u[,full_label := factor(full_label, levels = c(
"Anisotremus davidsonii\nxantic sargo",       "Bodianus pulcher\nCalifornia sheephead",     "Chromis punctipinnis\nblacksmith",          
"Girella nigricans\nopaleye",                 "Hypsypops rubicundus\nGaribaldi damselfish",                                   
"Paralabrax clathratus\nkelp bass",           "Paralabrax nebulifer\nbarred sand bass",   "Rhacochilus toxotes\nrubberlip seaperch",  "Stereolepis gigas\ngiant sea bass",  "Other"     ))]

#Zoom into individual artificial reef complexes


#Fish density
dat_fish_site_averages.AR <- dat_fish_site_averages[DepthZone %in% c("AR_SM","AR_PVR")]
dat_fish_site_averages.AR[,mean_AR_complex_wt_density_g_m2 := mean(mean_wt_density_g_m2),.(taxa, AR_Complex)] 
dat_fish_site_averages.AR[,mean_AR_complex_density_m2 := mean(mean_density_m2),.(taxa, AR_Complex)] 
dat_fish_averages.AR <- unique(dat_fish_site_averages.AR[,.(worms_id, taxa, common_name_final, kingdom, phylum, class, order, family, genus, rank, AR_Complex,mean_AR_complex_density_m2, mean_AR_complex_wt_density_g_m2)])
dat_fish_averages.AR[, Species_top5 := ifelse(frank(-mean_AR_complex_density_m2)<=5,taxa,"Other"), .(AR_Complex)]
dat_fish_averages.AR[, common_top5 := ifelse(frank(-mean_AR_complex_density_m2)<=5,common_name_final,""), .(AR_Complex)]
dat_fish_averages.AR[, summed_mean_AR_complex_density_m2 := sum(mean_AR_complex_density_m2), .(AR_Complex, Species_top5)]
dat_fishdensity_averages_AR.u <- unique(dat_fish_averages.AR[,.(Species_top5,common_top5, AR_Complex, summed_mean_AR_complex_density_m2)])
dat_fishdensity_averages_AR.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]

#Fish biomass
#AR only, and average by AR Complex
dat_fish_averages.AR[, Species_top5_biomass := ifelse(frank(-mean_AR_complex_wt_density_g_m2)<=5,taxa,"Other"), .(AR_Complex)]
dat_fish_averages.AR[, common_top5_biomass := ifelse(frank(-mean_AR_complex_wt_density_g_m2)<=5,common_name_final,""), .(AR_Complex)]
dat_fish_averages.AR[, summed_mean_AR_complex_density_wt_m2 := sum(mean_AR_complex_wt_density_g_m2), .(AR_Complex, Species_top5)]
dat_fishbiomass_averages_AR.u <- unique(dat_fish_averages.AR[,.(Species_top5_biomass,common_top5_biomass, AR_Complex, summed_mean_AR_complex_density_wt_m2)])
dat_fishbiomass_averages_AR.u[,full_label_biomass := ifelse(Species_top5_biomass == "Other","Other",paste0(Species_top5_biomass,"\n", common_top5_biomass))]

#Kelp density
dat_kelp_site_averages.AR <- dat_kelp_site_averages[DepthZone %in% c("AR_SM","AR_PVR")]
dat_kelp_site_averages.AR[,mean_AR_complex_density_m2 := mean(mean_density_m2),.(taxa, AR_Complex)] 
dat_kelp_averages.AR <- unique(dat_kelp_site_averages.AR[,.(worms_id, taxa, common_name_final, kingdom, phylum, class, order, family, genus, rank, AR_Complex,mean_AR_complex_density_m2)])
dat_kelp_averages.AR[, Species_top5 := ifelse(frank(-mean_AR_complex_density_m2)<=5,taxa,"Other"), .(AR_Complex)]
dat_kelp_averages.AR[, common_top5 := ifelse(frank(-mean_AR_complex_density_m2)<=5,common_name_final,""), .(AR_Complex)]
dat_kelp_averages.AR[, summed_mean_AR_complex_density_m2 := sum(mean_AR_complex_density_m2), .(AR_Complex, Species_top5)]
dat_kelpdensity_averages_AR.u <- unique(dat_kelp_averages.AR[,.(Species_top5,common_top5, AR_Complex, summed_mean_AR_complex_density_m2)])
dat_kelpdensity_averages_AR.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]

#Macroinvert density
dat_macroinvert_site_averages.AR <- dat_macroinvert_site_averages[DepthZone %in% c("AR_SM","AR_PVR")]
dat_macroinvert_site_averages.AR[,mean_AR_complex_density_m2 := mean(mean_density_m2),.(taxa, AR_Complex)] 
dat_macroinvert_averages.AR <- unique(dat_macroinvert_site_averages.AR[,.(worms_id, taxa, common_name_final, kingdom, phylum, class, order, family, genus, rank, AR_Complex,mean_AR_complex_density_m2)])
dat_macroinvert_averages.AR[, Species_top5 := ifelse(frank(-mean_AR_complex_density_m2)<=5,taxa,"Other"), .(AR_Complex)]
dat_macroinvert_averages.AR[, common_top5 := ifelse(frank(-mean_AR_complex_density_m2)<=5,common_name_final,""), .(AR_Complex)]
dat_macroinvert_averages.AR[, summed_mean_AR_complex_density_m2 := sum(mean_AR_complex_density_m2), .(AR_Complex, Species_top5)]
dat_macroinvertdensity_averages_AR.u <- unique(dat_macroinvert_averages.AR[,.(Species_top5,common_top5, AR_Complex, summed_mean_AR_complex_density_m2)])
dat_macroinvertdensity_averages_AR.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]


#how many sites in each complex?

unique(dat_fish_site_averages.AR[,.(Site, AR_Complex)])

#PVR = 18
#Santa Monica = 3
#Marina del Rey = 3
#Hermosa Beach = 1

############################################################
#Stacked barplots for all sites together, if not in top 5 species, summed into OTHER category
############################################################

#Make color palettes
color_palette_fishdensity <- c( 
                   "#D178A4", #rubberlip sea perch
                   "#8DD3C7",#opaleye
                   "#BEBADA" ,#rock wrasse
                   "#CCCC8F",#California sheephead
                   "#FDB462", #garibaldi
                   "#FCCDE5",#senorita
                   "forestgreen",#kelp bass
                   "#80B1D3",#barred sand bass
                   "#FB8072" ,#kelp perch
                   "#BC80BD" ,#bluebanded goby
                   "#CCEBC5" ,#blacksmith
                    "#D9D9D9")

color_palette_fishbiomass <- c("deepskyblue2",#sargo
                               "#D178A4", #rubberlip sea perch
                               "#FDB462" ,#garibaldi
                               "#E16A86",#giant sea bass
                               "#CCCC8F",#california sheephead
                               "#8DD3C7" ,#opaleye
                               "forestgreen",#kelp bass
                               "#80B1D3",#barred sand bass
                               "#CCEBC5",#blacksmith
                               "#D9D9D9" #other
                               )


color_palette_fishall <- c("#D178A4", #rubberlip sea perch
                           "#FB8072",#kelp perch
                           "#CCEBC5",#blacksmith
                           "#8DD3C7" ,#opaleye
                           "#BEBADA",#rock wrasse
                           "#CCCC8F",#california sheephead
                           "#FDB462" ,#garibaldi
                           "#FCCDE5" ,#senorita
                           "forestgreen",#kelp bass
                           "#80B1D3",#barred sand bass
                           "#BC80BD" ,#bluebanded goby
                           "deepskyblue2",#xantic sargo
                           "#E16A86",#giant sea bass
                           "#D9D9D9")#other


color_palette_macroinvert <- c("#8DD3C7",#stalked tunicate
                               "#CCCC8F",#starburst anenome
                               "blueviolet", #warty sea cucumber
                               "#BEBADA" ,#bat star
                               "indianred3",#brown gorgonian
                               "#80B1D3" ,#crowned urchin
                               "#FCCDE5" ,#golden gorgonian
                               "#BC80BD" ,#wavy turban snail
                               "#FDB462" ,#purple sea urchin
                               "#B3DE69" ,#red urchin
                               "lightskyblue",#tube dwelling anenome
                               "#FB8072",#red gorgonian
                               "#D9D9D9")


color_palette_kelp <- c("#8DD3C7",#southern sea palm
                        "#CCCC8F",#s. palmeri
                        "#BEBADA" ,#stalked kelp
                        "#FB8072",#chain bladder kelp
                        "#80B1D3" ,#fringed sea kelp
                        "#FDB462" ,#s. horneri
                        "#B3DE69" ,#golden kombu
                        "#FCCDE5" ,#giant kelp
                        "#D9D9D9")

#add new rank column
dat_fishdensity_averages.u[,rank := ifelse(full_label == "Other",400,frank(summed_mean_depthzone_density_m2))]

#fish density
fish_density_top5_stacked <- ggplot(dat_fishdensity_averages.u, aes(fill=reorder(full_label,rank), y=(summed_mean_depthzone_density_m2*100), x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "                   Depth zone                                            AR Complex", y = expression(paste("Average density per 100m" ^2)), fill = "Fish species") +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos Verdes", "Santa Monica\nBay")) +
  scale_fill_manual(values = color_palette_fishdensity) +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(aes(xintercept = 4.5), linewidth = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5))

ggsave(fish_density_top5_stacked, path = file.path("figures"), filename = "fish_density_top5_stacked.jpg", height = 6, width = 8, units = "in")

#add new rank column
dat_fishbiomass_averages.u[,summed_mean_depthzone_biomass_kg_100m := (summed_mean_depthzone_biomass_m2*100/1000)]
dat_fishbiomass_averages.u[,rank := ifelse(full_label == "Other",400,frank(summed_mean_depthzone_biomass_m2))]

#fish biomass
fish_biomass_top5_stacked <- ggplot(dat_fishbiomass_averages.u, aes(fill=reorder(full_label,rank), y=(summed_mean_depthzone_biomass_m2*100/1000), x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos Verdes", "Santa Monica\nBay")) +
  labs(x = "                   Depth zone                                            AR Complex", y = expression(paste("Average biomass in kg per 100m" ^2)), fill = "Fish species") +
  scale_fill_manual(values = color_palette_fishbiomass) +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(aes(xintercept = 4.5), linewidth = 1) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5))

ggsave(fish_biomass_top5_stacked, path = file.path("figures"), filename = "fish_biomass_top5_stacked.jpg", height = 6, width = 8, units = "in")

#merge fish plots
#dummy data.table

all_spp <-c("Rhacochilus toxotes\nrubberlip seaperch",
            "Brachyistius frenatus\nkelp perch"    ,
            "Chromis punctipinnis\nblacksmith",
"Girella nigricans\nopaleye",
"Halichoeres semicinctus\nrock wrasse"  ,
"Bodianus pulcher\nCalifornia sheephead",
"Hypsypops rubicundus\nGaribaldi damselfish",
"Oxyjulis californica\nsenorita"        ,
"Paralabrax clathratus\nkelp bass"    ,
"Paralabrax nebulifer\nbarred sand bass",
"Lythrypnus dalli\nbluebanded goby",
"Anisotremus davidsonii\nxantic sargo"     ,
"Stereolepis gigas\ngiant sea bass"  ,
"Other")

dummy_fish_dt <- data.table(`Fish species` = factor(all_spp, levels = c("Rhacochilus toxotes\nrubberlip seaperch",
                                                                        "Brachyistius frenatus\nkelp perch"    ,
                                                                        "Chromis punctipinnis\nblacksmith",
                                                                        "Girella nigricans\nopaleye",
                                                                        "Halichoeres semicinctus\nrock wrasse"  ,
                                                                        "Bodianus pulcher\nCalifornia sheephead",
                                                                        "Hypsypops rubicundus\nGaribaldi damselfish",
                                                                        "Oxyjulis californica\nsenorita"        ,
                                                                        "Paralabrax clathratus\nkelp bass"    ,
                                                                        "Paralabrax nebulifer\nbarred sand bass",
                                                                        "Lythrypnus dalli\nbluebanded goby",
                                                                        "Anisotremus davidsonii\nxantic sargo"     ,
                                                                        "Stereolepis gigas\ngiant sea bass"  ,
                                                                        "Other")), num = rep(1,14))

fish_leg <- get_legend(ggplot(dummy_fish_dt) +
                         geom_col(aes(x = `Fish species`, y = num, fill = `Fish species`)) +
                         scale_fill_manual(values = color_palette_fishall) +
                         theme_classic() +
                         theme(legend.position = "top",legend.direction = "horizontal"))
#merge w/o legends
fish_abundance_top5_stacked <- plot_grid(fish_density_top5_stacked + theme(legend.position = "none"), fish_biomass_top5_stacked + theme(legend.position = "none"), ncol = 2, labels = c("a.","b."))

#merge with legend
fish_abundance_top5_stacked.l <- plot_grid(fish_leg,fish_abundance_top5_stacked,ncol = 1, rel_heights = c(2,10))

ggsave(fish_abundance_top5_stacked.l, path = file.path("figures"), filename = "fish_abundance_top5_stacked.l.jpg", height = 6, width = 12, units = "in")


#Surprised giant sea bass shows up as top! look into this quickly (all from Santa Monica ARs)
ggplot(dat_fish_site_averages[taxa == "Stereolepis gigas"]) +
  geom_boxplot(aes(x = DepthZone, y = mean_density_m2)) +
  theme_classic()

ggplot(dat_fish_site_averages[taxa == "Stereolepis gigas"]) +
  geom_boxplot(aes(x = DepthZone, y = mean_wt_density_g_m2)) +
  theme_classic()
  


#macroinvert density
#add new rank column
dat_macroinvertdensity_averages.u[,rank := ifelse(full_label == "Other",400,frank(summed_mean_depthzone_density_m2))]

macroinvert_density_top5_stacked <- ggplot(dat_macroinvertdensity_averages.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "                   Depth zone                                            AR Complex", y = expression(paste("Average density per 100m" ^2)), fill = "Macroinvertebrate species") +
  scale_fill_manual(values = color_palette_macroinvert) +
    scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos Verdes","Santa Monica\nBay")) +
  geom_vline(aes(xintercept = 4.5), linewidth = 1) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5))

ggsave(macroinvert_density_top5_stacked, path = file.path("figures"), filename = "macroinvert_density_top5_stacked.jpg", height = 6, width = 8, units = "in")

#kelp density
#add new rank column
dat_kelp_averages.u[,rank := ifelse(full_label == "Other",400,frank(summed_mean_depthzone_density_m2))]

kelp_density_top5_stacked <- ggplot(dat_kelp_averages.u, aes(fill=reorder(full_label, rank), y=summed_mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "                   Depth zone                                            AR Complex", y = expression(paste("Average density per 100m" ^2)), fill = "Macroalgae species") +
  scale_fill_manual(values = color_palette_kelp) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","Palos Verdes","Santa Monica\nBay")) +
  geom_vline(aes(xintercept = 4.5), linewidth = 1) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5))

ggsave(kelp_density_top5_stacked, path = file.path("figures"), filename = "kelp_density_top5_stacked.jpg", height = 6, width = 8, units = "in")

#Stacked bar plots for Artificial reefs only, compare AR complexes
#fish density

#rank species by abundance
dat_fishdensity_averages_AR.u[,abun_total := sum(summed_mean_AR_complex_density_m2),Species_top5]
#add new rank column
dat_fishdensity_averages_AR.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]
#Change order of factors for AR complex to be latitudinal
dat_fishdensity_averages_AR.u[,AR_Complex := factor(AR_Complex,
                                                    levels = c("Santa Monica","Marina del Rey 2","Hermosa Beach","Palos Verdes Reef"),
                                                    labels = c("Santa Monica","Marina del Rey","Hermosa Beach","Palos Verdes Reef"))]

fish_AR_density_palette <- c("#BEBADA" ,#rock wrasse
                              "#D178A4", #rubberlip sea perch
                             "#FB8072" ,#kelp perch
                              "#8DD3C7",#opaleye
                              "#CCCC8F",#California sheephead
                             "#80B1D3",#barred sand bass,
                              "forestgreen" ,#kelp bass
                              "#CCEBC5" ,#blacksmith
                              "#D9D9D9") #other

fish_density_top5_stacked_AR <- ggplot(dat_fishdensity_averages_AR.u, aes(fill=reorder(full_label, rank), y=(summed_mean_AR_complex_density_m2*100), x=AR_Complex)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "AR complex", y = expression(paste("Average density per 100m" ^2)), fill = "Fish species") +
  scale_fill_manual(values = fish_AR_density_palette) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5))

#fish biomass
#rank species by abundance
#total abundance
dat_fishbiomass_averages_AR.u[,abun_total := sum(summed_mean_AR_complex_density_wt_m2),Species_top5_biomass]
#add new rank column
dat_fishbiomass_averages_AR.u[,rank := ifelse(full_label_biomass == "Other",400,frank(abun_total))]
#Change order of factors for AR complex to be latitudinal
dat_fishbiomass_averages_AR.u[,AR_Complex := factor(AR_Complex,
                                                    levels = c("Santa Monica","Marina del Rey 2","Hermosa Beach","Palos Verdes Reef"),
                                                    labels = c("Santa Monica","Marina del Rey","Hermosa Beach","Palos Verdes Reef"))]

fish_AR_biomass_palette <- c("#D178A4",#seaperch
                              "#8DD3C7",#opaleye
                             "forestgreen" ,#kelp bass
                              "#CCEBC5" ,#blacksmith
                              "#CCCC8F",#California sheephead
                             "#DC5725", #broomtail grouper
                              "#80B1D3",#barred sand bass
                             "#E16A86",#giant sea bass
                              "#D9D9D9") #other

fish_biomass_top5_stacked_AR <- ggplot(dat_fishbiomass_averages_AR.u, aes(fill=reorder(full_label_biomass, rank), y=(summed_mean_AR_complex_density_wt_m2/100*1000), x=AR_Complex)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "AR complex", y = expression(paste("Average biomass in kg per 100m" ^2)), fill = "Fish species") +
  scale_fill_manual(values = fish_AR_biomass_palette) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5))


#dummy data.table
all_spp_AR <-c( "Halichoeres semicinctus\nrock wrasse"    ,
                "Rhacochilus toxotes\nrubberlip seaperch"    ,
                "Girella nigricans\nopaleye",
               "Chromis punctipinnis\nblacksmith",
               "Bodianus pulcher\nCalifornia sheephead",
               "Paralabrax clathratus\nkelp bass"  ,
               "Paralabrax nebulifer\nbarred sand bass",
               "Stereolepis gigas\ngiant sea bass" ,
               "Mycteroperca xenarcha\nbroomtail grouper" ,"Other")

fishall_AR_palette <- c("#BEBADA" ,#rock wrasse
                        "#D178A4", #rubberlip sea perch
                        "#8DD3C7",#opaleye
                             "#CCEBC5" ,#blacksmith
                             "#CCCC8F",#California sheephead
                             "forestgreen" ,#kelp bass
                             "#80B1D3",#barred sand bass,
                             "#E16A86",#giant sea bass
                             "#DC5725", #broomtail grouper
                             "#D9D9D9") #other

dummy_fish_AR_dt <- data.table(`Fish species` = factor(all_spp_AR, levels = c( "Halichoeres semicinctus\nrock wrasse"    ,
                                                                            "Rhacochilus toxotes\nrubberlip seaperch"    ,
                                                                            "Girella nigricans\nopaleye",
                                                                            "Chromis punctipinnis\nblacksmith",
                                                                            "Bodianus pulcher\nCalifornia sheephead",
                                                                            "Paralabrax clathratus\nkelp bass"  ,
                                                                            "Paralabrax nebulifer\nbarred sand bass",
                                                                            "Stereolepis gigas\ngiant sea bass" ,
                                                                            "Mycteroperca xenarcha\nbroomtail grouper" ,
                                                                            "Other")), num = rep(1,10))

fish_AR_leg <- get_legend(ggplot(dummy_fish_AR_dt) +
                         geom_col(aes(x = `Fish species`, y = num, fill = `Fish species`)) +
                         scale_fill_manual(values = fishall_AR_palette) +
                         theme_classic() +
                           theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"), legend.position = "top") +
                           guides(fill = guide_legend(keywidth = 1, keyheight = 1, title.position = "top",title.hjust = 0.5)))

#merge artificial reef fish plots

fish_abundance_top5_stacked_AR <- plot_grid(fish_density_top5_stacked_AR+theme(legend.position = "null", axis.title.x = element_blank()), fish_biomass_top5_stacked_AR+theme(legend.position = "null", axis.title.x = element_blank()), ncol = 2, labels = c("a.","b."))
fish_abundance_top5_stacked_AR.l <- plot_grid(fish_AR_leg,fish_abundance_top5_stacked_AR, ncol = 1,rel_heights = c(2,10))

#kelp density
#rank species by abundance
#total abundance
dat_kelpdensity_averages_AR.u[,abun_total := sum(summed_mean_AR_complex_density_m2),Species_top5]
#add new rank column
dat_kelpdensity_averages_AR.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]
#Change order of factors for AR complex to be latitudinal
dat_kelpdensity_averages_AR.u[,AR_Complex := factor(AR_Complex,
                                                    levels = c("Santa Monica","Marina del Rey 2","Hermosa Beach","Palos Verdes Reef"),
                                                    labels = c("Santa Monica","Marina del Rey","Hermosa Beach","Palos Verdes Reef"))]

kelp_AR_density_palette  <- c(
                              "#80B1D3" ,#fringed sea kelp,
                              "#8DD3C7",#southern sea palm
                              "#BEBADA" ,#stalked kelp
                              "#FB8072",#chain bladder kelp
                              "#FCCDE5" ,#giant kelp
                              "#B3DE69" ,#golden kombu
                              "#D9D9D9")

kelp_density_top5_stacked_AR <- ggplot(dat_kelpdensity_averages_AR.u, aes(fill=reorder(full_label, rank), y=(summed_mean_AR_complex_density_m2*100), x=AR_Complex)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "AR complex", y = expression(paste("Average density per 100m" ^2)), fill = "Macroalgae species") +
  scale_fill_manual(values = kelp_AR_density_palette) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"), legend.position = "top") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, title.position = "top",title.hjust = 0.5))

#macroinvert density
#rank species by abundance
#total abundance
dat_macroinvertdensity_averages_AR.u[,abun_total := sum(summed_mean_AR_complex_density_m2),Species_top5]
#add new rank column
dat_macroinvertdensity_averages_AR.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]
#Change order of factors for AR complex to be latitudinal
dat_macroinvertdensity_averages_AR.u[,AR_Complex := factor(AR_Complex,
                                                    levels = c("Santa Monica","Marina del Rey 2","Hermosa Beach","Palos Verdes Reef"),
                                                    labels = c("Santa Monica","Marina del Rey","Hermosa Beach","Palos Verdes Reef"))]

macroinvert_AR_density_palette <- c( "#BC80BD" ,#wavy turban snail
                                     "#8DD3C7",#stalked tunicate
                                     "#FDB462" ,#purple sea urchin
                                     "blueviolet", #warty sea cucumber
                                     "cyan4", #chestnut cowry
                                     "indianred3",#brown gorgonian
                                     "burlywood", #Kellets whelk
                                     "lightskyblue",#tube dwelling anenome
                                     "#B3DE69" ,#red urchin
                                     "#FB8072",#red gorgonian
                                     "#FCCDE5" ,#golden gorgonian
                                     "#D9D9D9")

macroinvert_density_top5_stacked_AR <- ggplot(dat_macroinvertdensity_averages_AR.u, aes(fill=reorder(full_label, rank), y=(summed_mean_AR_complex_density_m2*100), x=AR_Complex)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "AR complex", y = expression(paste("Average density per 100m" ^2)), fill = "Macroinvertebrate species") +
  scale_fill_manual(values = macroinvert_AR_density_palette) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"), legend.position = "top") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, title.position = "top",title.hjust = 0.5))

#merge these all together (zoomed in on artificial reefs))

#kelp and macroinvert
kelp_macroinvert_density_top5_stacked_AR <- plot_grid(kelp_density_top5_stacked_AR,
                                                      macroinvert_density_top5_stacked_AR,
                                                      ncol = 2, labels = c("c.","d."))

#fish, kelp, macroinverts on ARs
all_spp_top5_stacked_AR <- plot_grid(fish_abundance_top5_stacked_AR.l, kelp_macroinvert_density_top5_stacked_AR, nrow = 2, rel_heights = c(1,1.1))

ggsave(all_spp_top5_stacked_AR, path = "figures", filename = "all_spp_top5_stacked_AR.jpg", height = 9, width = 14, unit ="in")

############################################################
#Stacked barplots for all sites together, all species grouped by taxonomy
############################################################


#some strange order identification, will specify for following species
dat_fish_averages[taxa == "Anisotremus davidsonii", order := "Perciformes"]
dat_fish_averages[taxa == "Caulolatilus princeps", order := "Perciformes"]
dat_fish_averages[taxa == "Cheilotrema saturnum", order := "Acanthuriformes"]
dat_fish_averages[taxa == "Halichoeres semicinctus", order := "Labriformes"]
dat_fish_averages[taxa == "Oxyjulis californica", order := "Labriformes"]
dat_fish_averages[taxa == "Pristigenys serrula", order := "Perciformes"]
dat_fish_averages[taxa == "Bodianus pulcher", order := "Labriformes"]
dat_fish_averages[taxa == "Brachyistius frenatus", order := "Perciformes"]
dat_fish_averages[taxa == "Chromis punctipinnis", order := "Perciformes"]
dat_fish_averages[taxa == "Cymatogaster aggregata", order := "Perciformes"]
dat_fish_averages[taxa == "Embiotoca jacksoni", order := "Perciformes"]
dat_fish_averages[taxa == "Embiotoca lateralis", order := "Perciformes"]
dat_fish_averages[taxa == "Hyperprosopon argenteum", order := "Perciformes"]
dat_fish_averages[taxa == "Hypsurus caryi", order := "Perciformes"]
dat_fish_averages[taxa == "Hypsypops rubicundus", order := "Perciformes"]
dat_fish_averages[taxa == "Micrometrus minimus", order := "Perciformes"]
dat_fish_averages[taxa == "Phanerodon atripes", order := "Perciformes"]
dat_fish_averages[taxa == "Phanerodon furcatus", order := "Perciformes"]
dat_fish_averages[taxa == "Rhacochilus toxotes", order := "Perciformes"]
dat_fish_averages[taxa == "Zalembius rosaceus", order := "Perciformes"]
dat_fish_averages[taxa == "Damalichthys vacca", order := "Perciformes"]
dat_fish_averages[taxa == "Hermosilla azurea", order := "Centrarchiformes"] 
dat_fish_averages[taxa == "Lythrypnus sp.", order := "Gobiiformes"] 
dat_fish_averages[taxa == "Syngnathus leptorhynchus", order := "Syngnathiformes"]


#Sum density across orders
dat_fish_summed_order <- dat_fish_averages[,.(mean_depthzone_density_m2_summed_by_order = sum(mean_depthzone_density_m2),
                                               mean_depthzone_wt_density_g_m2_summed_by_order=sum(mean_depthzone_wt_density_g_m2)),
                                            .(order,DepthZone)] 

#generate 16 colors
palette1 <- brewer.pal(8, "Set1")
palette2 <- brewer.pal(8, "Set2")
combined_palette <- c(palette1, palette2)

#need to rearrange colors

combined_palette <- c("#E41A1C", "#377EB8", "#4DAF4A","#F781BF", "#984EA3" ,"#FF7F00" ,"#FFFF33" ,"#66C2A5" ,"#FC8D62" ,"#8DA0CB", "#E78AC3" ,"#B3B3B3" ,"#A65628","#FFD92F", "#E5C494" ,"#A6D854")

#fish density
 fish_density_order_stacked <- ggplot(dat_fish_summed_order,
                                     aes(fill=order, y=mean_depthzone_density_m2_summed_by_order*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Summed average density per 100m" ^2)), fill = "Order") +
  scale_fill_manual(values = combined_palette) +
   scale_y_continuous(expand = c(0,0)) +
  theme_classic()

ggsave(fish_density_order_stacked, path = file.path("figures"), filename = "fish_density_order_stacked.jpg",
       height = 4.5, width = 5.5, units = "in")

#fish biomass
fish_biomass_order_stacked <- ggplot(dat_fish_summed_order, aes(fill=order, y=mean_depthzone_wt_density_g_m2_summed_by_order*100/1000, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = bquote("Summed average biomass (kg per 100m"^2*")"), fill = "Order") +
  scale_fill_manual(values = combined_palette) +
  theme_classic()

ggsave(fish_biomass_order_stacked, path = file.path("figures"), filename = "fish_biomass_order_stacked.jpg",
       height = 4.5, width = 5.5, units = "in")


############################################################
#Now, average biomass and density by Island vs ARM vs Coast ####
############################################################

#New Column Identifying Island versus Mainland
dat_event.r[,type := ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")]
dat_fish_site_averages[,type := ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")]
dat_macroinvert_site_averages[,type := ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")]
dat_kelp_site_averages[,type := ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")]

#average by type and depth zone
dat_fish_averages_sitetype <- dat_fish_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2),
                                                        mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                                     .(taxa, common_name_final, DepthZone, type)] 
dat_macroinvert_averages_sitetype <- dat_macroinvert_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                                   .(taxa, common_name_final, DepthZone, type)]  
dat_kelp_averages_sitetype <- dat_kelp_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                     .(taxa, common_name_final, DepthZone, type)]  

#What fish species unique to island?
setdiff(unique(dat_fish_averages_sitetype[type == "Island" & mean_depthzone_density_m2>0]$taxa), unique(dat_fish_averages_sitetype[type == "Mainland" & mean_depthzone_density_m2>0]$taxa))
#What fish species unique to mainland?
setdiff(unique(dat_fish_averages_sitetype[type == "Mainland" & mean_depthzone_density_m2>0]$taxa), unique(dat_fish_averages_sitetype[type == "Island" & mean_depthzone_density_m2>0]$taxa))

#What macroinvert species unique to island?
setdiff(unique(dat_macroinvert_averages_sitetype[type == "Island" & mean_depthzone_density_m2>0]$taxa), unique(dat_macroinvert_averages_sitetype[type == "Mainland" & mean_depthzone_density_m2>0]$taxa))
#What macroinvert species unique to mainland?
setdiff(unique(dat_macroinvert_averages_sitetype[type == "Mainland" & mean_depthzone_density_m2>0]$taxa), unique(dat_macroinvert_averages_sitetype[type == "Island" & mean_depthzone_density_m2>0]$taxa))

#What macroinvert species unique to island?
setdiff(unique(dat_kelp_averages_sitetype[type == "Island" & mean_depthzone_density_m2>0]$taxa), unique(dat_kelp_averages_sitetype[type == "Mainland" & mean_depthzone_density_m2>0]$taxa))
#What macroinvert species unique to mainland?
setdiff(unique(dat_kelp_averages_sitetype[type == "Mainland" & mean_depthzone_density_m2>0]$taxa), unique(dat_kelp_averages_sitetype[type == "Island" & mean_depthzone_density_m2>0]$taxa))



#Identify top 5 species, sum other into 'other' category
###FISH DENSITY#####

dat_fish_averages_sitetype[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, common_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, summed_mean_depthzone_sitetype_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, type, Species_top5, common_top5)]
dat_fishdensity_averages_sitetype.u <- unique(dat_fish_averages_sitetype[,.(Species_top5, type, DepthZone, summed_mean_depthzone_sitetype_density_m2, common_top5)])
dat_fishdensity_averages_sitetype.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]
dat_fishdensity_averages_sitetype.u[,DepthZone := factor(DepthZone,levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"), labels = c("Depth zone\nInner","Depth zone\nMiddle","Depth zone\nOuter","Depth zone\nDeep","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
dat_fishdensity_averages_sitetype.u <- spp_color_key[dat_fishdensity_averages_sitetype.u, on = c("species_name" = "Species_top5")]
dat_fishdensity_averages_sitetype.u[,rel_abun := round(summed_mean_depthzone_sitetype_density_m2/sum(summed_mean_depthzone_sitetype_density_m2),2),.(DepthZone, type)]


####KELP######
dat_kelp_averages_sitetype[, Species_top5_kelp := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(DepthZone,type)]
dat_kelp_averages_sitetype[, common_top5_kelp := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(DepthZone,type)]
dat_kelp_averages_sitetype[, summed_mean_depthzone_sitetype_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, type, Species_top5_kelp, common_top5_kelp)]
dat_kelpdensity_averages_sitetype.u <- unique(dat_kelp_averages_sitetype[,.(Species_top5_kelp, type, DepthZone, summed_mean_depthzone_sitetype_density_m2, common_top5_kelp)])
dat_kelpdensity_averages_sitetype.u[,full_label := ifelse(Species_top5_kelp == "Other","Other",paste0(Species_top5_kelp,"\n", common_top5_kelp))]
dat_kelpdensity_averages_sitetype.u[,DepthZone := factor(DepthZone,levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"), labels = c("Depth zone\nInner","Depth zone\nMiddle","Depth zone\nOuter","Depth zone\nDeep","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
dat_kelpdensity_averages_sitetype.u <- spp_color_key[dat_kelpdensity_averages_sitetype.u, on = c("species_name" = "Species_top5_kelp")]


#####MACRO######
dat_macroinvert_averages_sitetype[, Species_top5_macroinvert := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(DepthZone,type)]
dat_macroinvert_averages_sitetype[, common_top5_macroinvert := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(DepthZone,type)]
dat_macroinvert_averages_sitetype[, summed_mean_depthzone_sitetype_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, type, Species_top5_macroinvert, common_top5_macroinvert)]
dat_macroinvertdensity_averages_sitetype.u <- unique(dat_macroinvert_averages_sitetype[,.(Species_top5_macroinvert, type, DepthZone, summed_mean_depthzone_sitetype_density_m2, common_top5_macroinvert)])
dat_macroinvertdensity_averages_sitetype.u[,full_label := ifelse(Species_top5_macroinvert == "Other","Other",paste0(Species_top5_macroinvert,"\n", common_top5_macroinvert))]
dat_macroinvertdensity_averages_sitetype.u[,DepthZone := factor(DepthZone,levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"), labels = c("Depth zone\nInner","Depth zone\nMiddle","Depth zone\nOuter","Depth zone\nDeep","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
dat_macroinvertdensity_averages_sitetype.u <- spp_color_key[dat_macroinvertdensity_averages_sitetype.u, on = c("species_name" = "Species_top5_macroinvert")]

######FISHBIOMASS######
dat_fish_averages_sitetype[, Species_top5_fishbiomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,taxa,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, common_top5_fishbiomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,common_name_final,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, summed_mean_depthzone_sitetype_biomass_m2 := sum(mean_depthzone_wt_density_g_m2), .(DepthZone, type, Species_top5_fishbiomass, common_top5_fishbiomass)]
dat_fishbiomass_averages_sitetype.u <- unique(dat_fish_averages_sitetype[,.(Species_top5_fishbiomass, type, DepthZone, summed_mean_depthzone_sitetype_biomass_m2, common_top5_fishbiomass)])
dat_fishbiomass_averages_sitetype.u[,full_label := ifelse(Species_top5_fishbiomass == "Other","Other",paste0(Species_top5_fishbiomass,"\n", common_top5_fishbiomass))]
dat_fishbiomass_averages_sitetype.u[,DepthZone := factor(DepthZone,levels = c("Inner","Middle","Outer","Deep","AR_PVR","AR_SM"), labels = c("Depth zone\nInner","Depth zone\nMiddle","Depth zone\nOuter","Depth zone\nDeep","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
dat_fishbiomass_averages_sitetype.u <- spp_color_key[dat_fishbiomass_averages_sitetype.u, on = c("species_name" = "Species_top5_fishbiomass")]

############################################################
#Stacked barplots for average biomass and density by site location (island/mainland) and depth, only top 5 species
############################################################
#fish density
#rank to order plot correctly
#total abundance
dat_fishdensity_averages_sitetype.u[,abun_total := sum(summed_mean_depthzone_sitetype_density_m2),species_name]

#add new rank column
dat_fishdensity_averages_sitetype.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]

#set order of hex by rank
dat_fishdensity_averages_sitetype.u[,hex := reorder(hex,rank)]

fish_density_top5_sitetype_stacked <- ggplot(dat_fishdensity_averages_sitetype.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_density_m2*100, x=type)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "", y = expression(paste("Average density per 100m" ^2)), fill = "Fish\nspecies") +
  scale_fill_manual(values = levels(dat_fishdensity_averages_sitetype.u$hex)) +
  scale_y_continuous(expand = c(0,0.1)) +
  facet_grid(~DepthZone, scales = "free_x", space = "free")+
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))

ggsave(fish_density_top5_sitetype_stacked, path = file.path("figures"), filename = "fish_density_top5_sitetype_stacked.jpg", height = 5, width = 10, units = "in")

#fish biomass
#rank to order plot correctly
#total abundance
dat_fishbiomass_averages_sitetype.u[,abun_total := sum(summed_mean_depthzone_sitetype_biomass_m2),species_name]
#add new rank column
dat_fishbiomass_averages_sitetype.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]

#set order of hex by rank
dat_fishbiomass_averages_sitetype.u[,hex := reorder(hex,rank)]

fish_biomass_top5_sitetype_stacked <- ggplot(dat_fishbiomass_averages_sitetype.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_biomass_m2*100/1000, x=type)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "", y = bquote("Average biomass (kg per 100m"^2*")"), fill = "Fish\nspecies") +
  scale_fill_manual(values =levels(dat_fishbiomass_averages_sitetype.u$hex)) +
  scale_y_continuous(expand = c(0,0.1)) +
  facet_grid(~DepthZone, scales = "free_x", space = "free") +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))

ggsave(fish_biomass_top5_sitetype_stacked, path = file.path("figures"), filename = "fish_biomass_top5_sitetype_stacked.jpg", height = 5, width = 10, units = "in")

#Merge fish
fish_top5_sitetype_stacked <- cowplot::plot_grid(fish_density_top5_sitetype_stacked, fish_biomass_top5_sitetype_stacked, ncol = 1, labels = c("a.","b."))

ggsave(fish_top5_sitetype_stacked, path = file.path("figures"), filename = "fish_top5_sitetype_stacked.jpg", height = 10, width = 10, units = "in")


#macroinvert density

#rank to order plot correctly
#total abundance
dat_macroinvertdensity_averages_sitetype.u[,abun_total := sum(summed_mean_depthzone_sitetype_density_m2),species_name]
#add new rank column
dat_macroinvertdensity_averages_sitetype.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]

#set order of hex by rank
dat_macroinvertdensity_averages_sitetype.u[,hex := reorder(hex,rank)]


macroinvert_density_top5_sitetype_stacked <- ggplot(dat_macroinvertdensity_averages_sitetype.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_density_m2*100, x=type)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "", y = expression(paste("Average density per 100m" ^2)), fill = "Macroinvertebrate\nspecies") +
  scale_fill_manual(values =levels(dat_macroinvertdensity_averages_sitetype.u$hex)) +
  scale_y_continuous(expand = c(0,0.1)) +
  facet_grid(~DepthZone, scales = "free_x", space = "free") +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))


ggsave(macroinvert_density_top5_sitetype_stacked, path = file.path("figures"), filename = "macroinvert_top5_sitetype_stacked.jpg", height = 5, width = 10, units = "in")

#kelp density
#rank to order plot correctly
#total abundance
dat_kelpdensity_averages_sitetype.u[,abun_total := sum(summed_mean_depthzone_sitetype_density_m2),species_name]
#add new rank column
dat_kelpdensity_averages_sitetype.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]


#set order of hex by rank
dat_kelpdensity_averages_sitetype.u[,hex := reorder(hex,rank)]

#kelp density
kelp_density_top5_sitetype_stacked <- ggplot(dat_kelpdensity_averages_sitetype.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_density_m2*100, x=type)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "", y = expression(paste("Average density per 100m" ^2)), fill = "Macroalgae\nspecies") +
  scale_fill_manual(values = levels(dat_kelpdensity_averages_sitetype.u$hex)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(~DepthZone, scales = "free_x", space = "free") +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))

ggsave(kelp_density_top5_sitetype_stacked, path = file.path("figures"), filename = "kelp_density_top5_sitetype_stacked.jpg", height = 5, width = 10, units = "in")

############################################################
#Now, average biomass and density by Island vs ARM but ***NO depth zones***
############################################################

#New Column Identifying Island versus Mainland
dat_fish_site_averages[,type_ARsplit := ifelse(DepthZone %in% c("AR_PVR", "AR_SM"),DepthZone,as.character(type))]
dat_macroinvert_site_averages[,type_ARsplit := ifelse(DepthZone %in% c("AR_PVR", "AR_SM"),DepthZone,as.character(type))]
dat_kelp_site_averages[,type_ARsplit := ifelse(DepthZone %in% c("AR_PVR", "AR_SM"),DepthZone,as.character(type))]

#average by type and depth zone
dat_fish_averages_sitetype_only_OUTERDEEP <- dat_fish_site_averages[!(DepthZone %in% c("Inner","Middle")), .(mean_depthzone_density_m2 = mean(mean_density_m2),
                                                                                                   mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                                               .(taxa, common_name_final, type_ARsplit)]

dat_fish_averages_sitetype_only <- dat_fish_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2),
                                                        mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                                     .(taxa, common_name_final, type_ARsplit)] 

dat_macroinvert_averages_sitetype_only_OUTERDEEP <- dat_macroinvert_site_averages[!(DepthZone %in% c("Inner","Middle")),.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                                             .(taxa, common_name_final, type_ARsplit)]  
dat_macroinvert_averages_sitetype_only <- dat_macroinvert_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                                   .(taxa, common_name_final, type_ARsplit)]  

dat_kelp_averages_sitetype_only_OUTERDEEP <- dat_kelp_site_averages[!(DepthZone %in% c("Inner","Middle")),.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                               .(taxa, common_name_final, type_ARsplit)]  
dat_kelp_averages_sitetype_only <- dat_kelp_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                     .(taxa, common_name_final, type_ARsplit)]  


#Identify top 5 species, sum other into 'other' category
###FISH DENSITY#####

dat_fish_averages_sitetype_only[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(type_ARsplit)]
dat_fish_averages_sitetype_only[, common_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(type_ARsplit)]
dat_fish_averages_sitetype_only[, summed_mean_depthzone_sitetype_only_density_m2 := sum(mean_depthzone_density_m2), .(type_ARsplit, Species_top5, common_top5)]
dat_fishdensity_averages_sitetype_only.u <- unique(dat_fish_averages_sitetype_only[,.(Species_top5, type_ARsplit, summed_mean_depthzone_sitetype_only_density_m2, common_top5)])
dat_fishdensity_averages_sitetype_only.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]
dat_fishdensity_averages_sitetype_only.u[,type_ARsplit := factor(type_ARsplit,levels = c("Island","Mainland","AR_PVR","AR_SM"), labels = c("Natural island","Mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
dat_fishdensity_averages_sitetype_only.u <- spp_color_key[dat_fishdensity_averages_sitetype_only.u, on = c("species_name" = "Species_top5")]

  #Only outer and deep natural sites
  dat_fish_averages_sitetype_only_OUTERDEEP[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(type_ARsplit)]
  dat_fish_averages_sitetype_only_OUTERDEEP[, common_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(type_ARsplit)]
  dat_fish_averages_sitetype_only_OUTERDEEP[, summed_mean_depthzone_sitetype_only_density_m2 := sum(mean_depthzone_density_m2), .(type_ARsplit, Species_top5, common_top5)]
  dat_fishdensity_averages_sitetype_only_OUTERDEEP.u <- unique(dat_fish_averages_sitetype_only_OUTERDEEP[,.(Species_top5, type_ARsplit, summed_mean_depthzone_sitetype_only_density_m2, common_top5)])
  dat_fishdensity_averages_sitetype_only_OUTERDEEP.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]
  dat_fishdensity_averages_sitetype_only_OUTERDEEP.u[,type_ARsplit := factor(type_ARsplit,levels = c("Island","Mainland","AR_PVR","AR_SM"), labels = c("Natural island\nOuter/Deep","Mainland\nOuter/Deep","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
  dat_fishdensity_averages_sitetype_only_OUTERDEEP.u <- spp_color_key[dat_fishdensity_averages_sitetype_only_OUTERDEEP.u, on = c("species_name" = "Species_top5")]
  

####KELP######
dat_kelp_averages_sitetype_only[, Species_top5_kelp := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(type_ARsplit)]
dat_kelp_averages_sitetype_only[, common_top5_kelp := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(type_ARsplit)]
dat_kelp_averages_sitetype_only[, summed_mean_depthzone_sitetype_only_density_m2 := sum(mean_depthzone_density_m2), .( type_ARsplit, Species_top5_kelp, common_top5_kelp)]
dat_kelpdensity_averages_sitetype_only.u <- unique(dat_kelp_averages_sitetype_only[,.(Species_top5_kelp, type_ARsplit,  summed_mean_depthzone_sitetype_only_density_m2, common_top5_kelp)])
dat_kelpdensity_averages_sitetype_only.u[,full_label := ifelse(Species_top5_kelp == "Other","Other",paste0(Species_top5_kelp,"\n", common_top5_kelp))]
dat_kelpdensity_averages_sitetype_only.u[,type_ARsplit := factor(type_ARsplit,levels = c("Island","Mainland","AR_PVR","AR_SM"), labels = c("Natural island","Mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
dat_kelpdensity_averages_sitetype_only.u <- spp_color_key[dat_kelpdensity_averages_sitetype_only.u, on = c("species_name" = "Species_top5_kelp")]

  #Only outer and deep natural sites
  dat_kelp_averages_sitetype_only_OUTERDEEP[, Species_top5_kelp := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(type_ARsplit)]
  dat_kelp_averages_sitetype_only_OUTERDEEP[, common_top5_kelp := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(type_ARsplit)]
  dat_kelp_averages_sitetype_only_OUTERDEEP[, summed_mean_depthzone_sitetype_only_density_m2 := sum(mean_depthzone_density_m2), .( type_ARsplit, Species_top5_kelp, common_top5_kelp)]
  dat_kelpdensity_averages_sitetype_only_OUTERDEEP.u <- unique(dat_kelp_averages_sitetype_only_OUTERDEEP[,.(Species_top5_kelp, type_ARsplit,  summed_mean_depthzone_sitetype_only_density_m2, common_top5_kelp)])
  dat_kelpdensity_averages_sitetype_only_OUTERDEEP.u[,full_label := ifelse(Species_top5_kelp == "Other","Other",paste0(Species_top5_kelp,"\n", common_top5_kelp))]
  dat_kelpdensity_averages_sitetype_only_OUTERDEEP.u[,type_ARsplit := factor(type_ARsplit,levels = c("Island","Mainland","AR_PVR","AR_SM"), labels = c("Natural island","Mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
  dat_kelpdensity_averages_sitetype_only_OUTERDEEP.u <- spp_color_key[dat_kelpdensity_averages_sitetype_only_OUTERDEEP.u, on = c("species_name" = "Species_top5_kelp")]

#####MACRO######
dat_macroinvert_averages_sitetype_only[, Species_top5_macroinvert := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(type_ARsplit)]
dat_macroinvert_averages_sitetype_only[, common_top5_macroinvert := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(type_ARsplit)]
dat_macroinvert_averages_sitetype_only[, summed_mean_depthzone_sitetype_only_density_m2 := sum(mean_depthzone_density_m2), .(type_ARsplit, Species_top5_macroinvert, common_top5_macroinvert)]
dat_macroinvertdensity_averages_sitetype_only.u <- unique(dat_macroinvert_averages_sitetype_only[,.(Species_top5_macroinvert, type_ARsplit,summed_mean_depthzone_sitetype_only_density_m2, common_top5_macroinvert)])
dat_macroinvertdensity_averages_sitetype_only.u[,full_label := ifelse(Species_top5_macroinvert == "Other","Other",paste0(Species_top5_macroinvert,"\n", common_top5_macroinvert))]
dat_macroinvertdensity_averages_sitetype_only.u[,type_ARsplit := factor(type_ARsplit,levels = c("Island","Mainland","AR_PVR","AR_SM"), labels = c("Natural island","Mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
dat_macroinvertdensity_averages_sitetype_only.u <- spp_color_key[dat_macroinvertdensity_averages_sitetype_only.u, on = c("species_name" = "Species_top5_macroinvert")]

  #Only outer and deep natural reef zones
  dat_macroinvert_averages_sitetype_only_OUTERDEEP[, Species_top5_macroinvert := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(type_ARsplit)]
  dat_macroinvert_averages_sitetype_only_OUTERDEEP[, common_top5_macroinvert := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(type_ARsplit)]
  dat_macroinvert_averages_sitetype_only_OUTERDEEP[, summed_mean_depthzone_sitetype_only_density_m2 := sum(mean_depthzone_density_m2), .(type_ARsplit, Species_top5_macroinvert, common_top5_macroinvert)]
  dat_macroinvertdensity_averages_sitetype_only_OUTERDEEP.u <- unique(dat_macroinvert_averages_sitetype_only_OUTERDEEP[,.(Species_top5_macroinvert, type_ARsplit,summed_mean_depthzone_sitetype_only_density_m2, common_top5_macroinvert)])
  dat_macroinvertdensity_averages_sitetype_only_OUTERDEEP.u[,full_label := ifelse(Species_top5_macroinvert == "Other","Other",paste0(Species_top5_macroinvert,"\n", common_top5_macroinvert))]
  dat_macroinvertdensity_averages_sitetype_only_OUTERDEEP.u[,type_ARsplit := factor(type_ARsplit,levels = c("Island","Mainland","AR_PVR","AR_SM"), labels = c("Natural island","Mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
  dat_macroinvertdensity_averages_sitetype_only_OUTERDEEP.u <- spp_color_key[dat_macroinvertdensity_averages_sitetype_only_OUTERDEEP.u, on = c("species_name" = "Species_top5_macroinvert")]


######FISHBIOMASS######
dat_fish_averages_sitetype_only[, Species_top5_fishbiomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,taxa,"Other"),.(type_ARsplit)]
dat_fish_averages_sitetype_only[, common_top5_fishbiomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,common_name_final,"Other"),.(type_ARsplit)]
dat_fish_averages_sitetype_only[, summed_mean_depthzone_sitetype_only_biomass_m2 := sum(mean_depthzone_wt_density_g_m2), .(type_ARsplit, Species_top5_fishbiomass, common_top5_fishbiomass)]
dat_fishbiomass_averages_sitetype_only.u <- unique(dat_fish_averages_sitetype_only[,.(Species_top5_fishbiomass, type_ARsplit, summed_mean_depthzone_sitetype_only_biomass_m2, common_top5_fishbiomass)])
dat_fishbiomass_averages_sitetype_only.u[,full_label := ifelse(Species_top5_fishbiomass == "Other","Other",paste0(Species_top5_fishbiomass,"\n", common_top5_fishbiomass))]
dat_fishbiomass_averages_sitetype_only.u[,type_ARsplit := factor(type_ARsplit,levels = c("Island","Mainland","AR_PVR","AR_SM"), labels = c("Natural island","Mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
dat_fishbiomass_averages_sitetype_only.u <- spp_color_key[dat_fishbiomass_averages_sitetype_only.u, on = c("species_name" = "Species_top5_fishbiomass")]

    #ONLY OUTER AND DEEP NATURAL REEF ZONES
    dat_fish_averages_sitetype_only_OUTERDEEP[, Species_top5_fishbiomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,taxa,"Other"),.(type_ARsplit)]
    dat_fish_averages_sitetype_only_OUTERDEEP[, common_top5_fishbiomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,common_name_final,"Other"),.(type_ARsplit)]
    dat_fish_averages_sitetype_only_OUTERDEEP[, summed_mean_depthzone_sitetype_only_biomass_m2 := sum(mean_depthzone_wt_density_g_m2), .(type_ARsplit, Species_top5_fishbiomass, common_top5_fishbiomass)]
    dat_fishbiomass_averages_sitetype_only_OUTERDEEP.u <- unique(dat_fish_averages_sitetype_only_OUTERDEEP[,.(Species_top5_fishbiomass, type_ARsplit, summed_mean_depthzone_sitetype_only_biomass_m2, common_top5_fishbiomass)])
    dat_fishbiomass_averages_sitetype_only_OUTERDEEP.u[,full_label := ifelse(Species_top5_fishbiomass == "Other","Other",paste0(Species_top5_fishbiomass,"\n", common_top5_fishbiomass))]
    dat_fishbiomass_averages_sitetype_only_OUTERDEEP.u[,type_ARsplit := factor(type_ARsplit,levels = c("Island","Mainland","AR_PVR","AR_SM"), labels = c("Natural island","Mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica\nBay"))]
    dat_fishbiomass_averages_sitetype_only_OUTERDEEP.u <- spp_color_key[dat_fishbiomass_averages_sitetype_only_OUTERDEEP.u, on = c("species_name" = "Species_top5_fishbiomass")]


############################################################
#Stacked barplots for average biomass and density by site location (island/mainland) and depth, only top 5 species
############################################################
#fish density
#rank to order plot correctly
#total abundance
dat_fishdensity_averages_sitetype_only.u[,abun_total := sum(summed_mean_depthzone_sitetype_only_density_m2),species_name]

#add new rank column
dat_fishdensity_averages_sitetype_only.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]

#set order of hex by rank
dat_fishdensity_averages_sitetype_only.u[,hex := reorder(hex,rank)]

fish_density_top5_sitetype_only_stacked <- ggplot(dat_fishdensity_averages_sitetype_only.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_only_density_m2*100, x=type_ARsplit)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "", y = expression(paste("Average density per 100m" ^2)), fill = "Fish\nspecies") +
  scale_fill_manual(values = levels(dat_fishdensity_averages_sitetype_only.u$hex)) +
  scale_x_discrete(labels = c("Natural island","Natural mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica Bay")) +
  scale_y_continuous(expand = c(0,0.1)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"), legend.position = "top") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, title.position = "top",title.hjust = 0.5))

ggsave(fish_density_top5_sitetype_only_stacked, path = file.path("figures"), filename = "fish_density_top5_sitetype_only_stacked.jpg", height = 5, width = 10, units = "in")


    #OUTER DEEP ONLY
#rank to order plot correctly
#total abundance
dat_fishdensity_averages_sitetype_only_OUTERDEEP.u[,abun_total := sum(summed_mean_depthzone_sitetype_only_density_m2),species_name]

#add new rank column
dat_fishdensity_averages_sitetype_only_OUTERDEEP.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]

#set order of hex by rank
dat_fishdensity_averages_sitetype_only_OUTERDEEP.u[,hex := reorder(hex,rank)]

fish_density_top5_sitetype_only_stacked_OUTERDEEP <- ggplot(dat_fishdensity_averages_sitetype_only_OUTERDEEP.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_only_density_m2*100, x=type_ARsplit)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "", y = expression(paste("Average density per 100m" ^2)), fill = "Fish\nspecies") +
  scale_fill_manual(values = levels(dat_fishdensity_averages_sitetype_only_OUTERDEEP.u$hex)) +
  scale_x_discrete(labels = c("Natural island\nOuter/Deep","Natural mainland\nOuter/Deep","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica Bay")) +
  scale_y_continuous(expand = c(0,0.1)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"), legend.position = "top") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, title.position = "top",title.hjust = 0.5))

ggsave(fish_density_top5_sitetype_only_stacked_OUTERDEEP, path = file.path("figures"), filename = "fish_density_top5_sitetype_only_stacked_OUTERDEEP.jpg", height = 5, width = 10, units = "in")


#fish biomass
#rank to order plot correctly
#total abundance
dat_fishbiomass_averages_sitetype_only.u[,abun_total := sum(summed_mean_depthzone_sitetype_only_biomass_m2),species_name]
#add new rank column
dat_fishbiomass_averages_sitetype_only.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]

#set order of hex by rank
dat_fishbiomass_averages_sitetype_only.u[,hex := reorder(hex,rank)]

fish_biomass_top5_sitetype_only_stacked <- ggplot(dat_fishbiomass_averages_sitetype_only.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_only_biomass_m2*100/1000, x=type_ARsplit)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "", y = bquote("Average biomass (kg per 100m"^2*")"), fill = "Fish\nspecies") +
  scale_fill_manual(values =levels(dat_fishbiomass_averages_sitetype_only.u$hex)) +
    scale_x_discrete(labels = c("Natural island","Natural mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica Bay")) +
  scale_y_continuous(expand = c(0,0.1)) +
  theme_classic() +
  theme(legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, title.position = "top",title.hjust = 0.5))

ggsave(fish_biomass_top5_sitetype_only_stacked, path = file.path("figures"), filename = "fish_biomass_top5_sitetype_only_stacked.jpg", height = 5, width = 10, units = "in")

      #OUTER DEEP ONLY
      #total abundance
      dat_fishbiomass_averages_sitetype_only_OUTERDEEP.u[,abun_total := sum(summed_mean_depthzone_sitetype_only_biomass_m2),species_name]
      #add new rank column
      dat_fishbiomass_averages_sitetype_only_OUTERDEEP.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]
      
      #set order of hex by rank
      dat_fishbiomass_averages_sitetype_only_OUTERDEEP.u[,hex := reorder(hex,rank)]
      
      fish_biomass_top5_sitetype_only_stacked_OUTERDEEP <- ggplot(dat_fishbiomass_averages_sitetype_only_OUTERDEEP.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_only_biomass_m2*100/1000, x=type_ARsplit)) + 
        geom_bar(position="stack", stat="identity") +
        labs(x = "", y = bquote("Average biomass (kg per 100m"^2*")"), fill = "Fish\nspecies") +
        scale_fill_manual(values =levels(dat_fishbiomass_averages_sitetype_only_OUTERDEEP.u$hex)) +
        scale_x_discrete(labels = c("Natural island","Natural mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica Bay")) +
        scale_y_continuous(expand = c(0,0.1)) +
        theme_classic() +
        theme(legend.position = "top", legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
        guides(fill = guide_legend(keywidth = 1, keyheight = 1, title.position = "top",title.hjust = 0.5))
      
      ggsave(fish_biomass_top5_sitetype_only_stacked_OUTERDEEP, path = file.path("figures"), filename = "fish_biomass_top5_sitetype_only_stacked_OUTERDEEP.jpg", height = 5, width = 10, units = "in")


#Merge fish
fish_top5_sitetype_only_stacked <- cowplot::plot_grid(fish_density_top5_sitetype_only_stacked, fish_biomass_top5_sitetype_only_stacked, ncol = 2, labels = c("a.","b."), align = "hv")

ggsave(fish_top5_sitetype_only_stacked, path = file.path("figures"), filename = "fish_top5_sitetype_only_stacked.jpg", height = 10, width = 10, units = "in")

#Merge fish outerdeep only
fish_top5_sitetype_only_stacked_OUTERDEEP <- cowplot::plot_grid(fish_density_top5_sitetype_only_stacked_OUTERDEEP, fish_biomass_top5_sitetype_only_stacked_OUTERDEEP, ncol = 2, labels = c("a.","b."), align = "hv")

ggsave(fish_top5_sitetype_only_stacked_OUTERDEEP, path = file.path("figures"), filename = "fish_top5_sitetype_only_stacked_OUTERDEEP.jpg", height = 10, width = 10, units = "in")


#macroinvert density

#rank to order plot correctly
#total abundance
dat_macroinvertdensity_averages_sitetype_only.u[,abun_total := sum(summed_mean_depthzone_sitetype_only_density_m2),species_name]
#add new rank column
dat_macroinvertdensity_averages_sitetype_only.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]

#set order of hex by rank
dat_macroinvertdensity_averages_sitetype_only.u[,hex := reorder(hex,rank)]


macroinvert_density_top5_sitetype_only_stacked <- ggplot(dat_macroinvertdensity_averages_sitetype_only.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_only_density_m2*100, x=type_ARsplit)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "", y = expression(paste("Average density per 100m" ^2)), fill = "Macroinvertebrate\nspecies") +
  scale_fill_manual(values =levels(dat_macroinvertdensity_averages_sitetype_only.u$hex)) +
    scale_x_discrete(labels = c("Natural island","Natural mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica Bay")) +
  scale_y_continuous(expand = c(0,0.1)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"), legend.position = "top") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, title.position = "top",title.hjust = 0.5))

ggsave(macroinvert_density_top5_sitetype_only_stacked, path = file.path("figures"), filename = "macroinvert_top5_sitetype_only_stacked.jpg", height = 5, width = 10, units = "in")

    #OUTER DEEP ONLY
    #total abundance
    dat_macroinvertdensity_averages_sitetype_only_OUTERDEEP.u[,abun_total := sum(summed_mean_depthzone_sitetype_only_OUTERDEEP_density_m2),species_name]
    #add new rank column
    dat_macroinvertdensity_averages_sitetype_only_OUTERDEEP.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]
    
    #set order of hex by rank
    dat_macroinvertdensity_averages_sitetype_only_OUTERDEEP.u[,hex := reorder(hex,rank)]
    
    
    macroinvert_density_top5_sitetype_only_stacked_OUTERDEEP <- ggplot(dat_macroinvertdensity_averages_sitetype_only_OUTERDEEP.u,
                                                                       aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_only_density_m2*100, x=type_ARsplit)) + 
      geom_bar(position="stack", stat="identity") +
      labs(x = "", y = expression(paste("Average density per 100m" ^2)), fill = "Macroinvertebrate\nspecies") +
      scale_fill_manual(values =levels(dat_macroinvertdensity_averages_sitetype_only_OUTERDEEP.u$hex)) +
      scale_x_discrete(labels = c("Natural island","Natural mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica Bay")) +
      scale_y_continuous(expand = c(0,0.1)) +
      theme_classic() +
      theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"), legend.position = "top") +
      guides(fill = guide_legend(keywidth = 1, keyheight = 1, title.position = "top",title.hjust = 0.5))
    
    ggsave(macroinvert_density_top5_sitetype_only_stacked_OUTERDEEP, path = file.path("figures"), filename = "macroinvert_top5_sitetype_only_stacked_OUTERDEEP.jpg", height = 5, width = 10, units = "in")


#kelp density
#rank to order plot correctly
#total abundance
dat_kelpdensity_averages_sitetype_only.u[,abun_total := sum(summed_mean_depthzone_sitetype_only_density_m2),species_name]
#add new rank column
dat_kelpdensity_averages_sitetype_only.u[,rank := ifelse(full_label == "Other",400,frank(abun_total))]


#set order of hex by rank
dat_kelpdensity_averages_sitetype_only.u[,hex := reorder(hex,rank)]

#kelp density
kelp_density_top5_sitetype_only_stacked <- ggplot(dat_kelpdensity_averages_sitetype_only.u, aes(fill=reorder(full_label,rank), y=summed_mean_depthzone_sitetype_only_density_m2*100, x=type_ARsplit)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "", y = expression(paste("Average density per 100m" ^2)), fill = "Macroalgae\nspecies") +
  scale_x_discrete(labels = c("Natural island","Natural mainland","Artificial reef\nPalos Verdes","Artificial reef\nSanta Monica Bay")) +
  scale_fill_manual(values = levels(dat_kelpdensity_averages_sitetype_only.u$hex)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"), legend.position = "top") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1, title.position = "top",title.hjust = 0.5))

ggsave(kelp_density_top5_sitetype_only_stacked, path = file.path("figures"), filename = "kelp_density_top5_sitetype_only_stacked.jpg", height = 5, width = 10, units = "in")

#Merge these four figures

top5_sitetype_only_stacked <- plot_grid(fish_density_top5_sitetype_only_stacked, fish_biomass_top5_sitetype_only_stacked,
                                        kelp_density_top5_sitetype_only_stacked, macroinvert_density_top5_sitetype_only_stacked, ncol = 2, nrow = 2,
                                        labels = c("a.","b.","c.","d."))

ggsave(top5_sitetype_only_stacked, path = file.path("figures"), filename = "top5_sitetype_only_stacked.jpg", height = 12, width = 16, units = "in")


