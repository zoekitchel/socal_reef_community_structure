# CREATION DATE 28 Jan 2024
# MODIFIED DATE 29 Jun 2024

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

#number of sites per depth zone
dat_event.r[,number_sites_depthzone := uniqueN(Site),.(DepthZone)]

number_sites_depthzone <- unique(dat_event.r[,.(DepthZone, number_sites_depthzone)])

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
dat_fishdensity_averages_deep_ar.u <- unique(dat_fish_averages[,.(Species_top5,common_top5, DepthZone, summed_mean_depthzone_density_m2)])
dat_fishdensity_averages_deep_ar.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]

#manually set factor order for plotting
dat_fishdensity_averages_deep_ar.u[,full_label := factor(full_label, levels = c(
  "Bodianus pulcher\nCalifornia sheephead", "Brachyistius frenatus\nkelp perch","Chromis punctipinnis\nblacksmith",
 "Girella nigricans\nopaleye", "Halichoeres semicinctus\nrock wrasse",
 "Hypsypops rubicundus\nGaribaldi damselfish" ,"Lythrypnus dalli\nbluebanded goby" ,
 "Oxyjulis californica\nsenorita",
 "Paralabrax clathratus\nkelp bass", "Paralabrax nebulifer\nbarred sand bass","Other"))]

#macroinvert density
dat_macroinvert_averages[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"), .(DepthZone)]
dat_macroinvert_averages[, common_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,""), .(DepthZone)]
dat_macroinvert_averages[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
dat_macroinvertdensity_averages_deep_ar.u <- unique(dat_macroinvert_averages[,.(Species_top5,common_top5, DepthZone, summed_mean_depthzone_density_m2)])
dat_macroinvertdensity_averages_deep_ar.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]

#manually set factor order for plotting
dat_macroinvertdensity_averages_deep_ar.u[,full_label := factor(full_label, levels = c(
"Anthopleura sola\nstarburst anenome","Centrostephanus coronatus\ncrowned urchin" ,
"Leptogorgia chilensis\nred gorgonian", "Megastraea undosa\nwavy turban snail",
"Mesocentrotus franciscanus\nred urchin", "Muricea californica\ngolden gregorian" ,
"Patiria miniata\nbat star" ,
"Strongylocentrotus purpuratus\npurple sea urchin", "Styela montereyensis\nstalked tunicate","Other"))]


#Kelp density
dat_kelp_averages[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"), .(DepthZone)]
dat_kelp_averages[, common_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,""), .(DepthZone)]
dat_kelp_averages[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
dat_kelp_averages_deep_ar.u <- unique(dat_kelp_averages[,.(Species_top5, common_top5, DepthZone, summed_mean_depthzone_density_m2)])
dat_kelp_averages_deep_ar.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]

#manually set factor order for plotting
dat_kelp_averages_deep_ar.u[,full_label := factor(full_label, levels = c(
"Agarum fimbriatum","Eisenia arborea\nsouthern sea palm", "Laminaria farlowii\ngolden kombu",
"Macrocystis pyrifera\ngiant kelp", "Pterygophora californica\nstalked kelp",
"Sargassum horneri","Sargassum palmeri","Stephanocystis spp." ,
"Other"), labels = c(
  "Agarum fimbriatum",#"#8DD3C7"
  "Eisenia arborea\nsouthern sea palm",#"#CCCC8F"
  "Laminaria farlowii\ngolden kombu",#"#BEBADA"
  "Macrocystis pyrifera\ngiant kelp",#"#FB8072"
  "Pterygophora californica\nstalked kelp",#"#80B1D3"
  "Sargassum horneri",#"#FDB462" 
  "Sargassum palmeri",#"#B3DE69" 
  "Stephanocystis spp." ,#"#FCCDE5"
  "Other"))] #"#D9D9D9"


#Add column for whether species is in top 5 species (by biomass) for that zone, and then sum biomass for all others
#fish biomass
dat_fish_averages[, Species_top5_biomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,taxa,"Other"), .(DepthZone)]
dat_fish_averages[, common_top5_biomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,common_name_final,""), .(DepthZone)]
dat_fish_averages[, summed_mean_depthzone_biomass_m2 := sum(mean_depthzone_wt_density_g_m2), .(DepthZone, Species_top5_biomass)]
dat_fishbiomass_averages_deep_ar.u <- unique(dat_fish_averages[,.(Species_top5_biomass,common_top5_biomass, DepthZone, summed_mean_depthzone_biomass_m2)])
dat_fishbiomass_averages_deep_ar.u[,full_label := ifelse(Species_top5_biomass == "Other","Other",paste0(Species_top5_biomass,"\n", common_top5_biomass))]

#manually set factor order for plotting
dat_fishbiomass_averages_deep_ar.u[,full_label := factor(full_label, levels = c(
"Anisotremus davidsonii\nxantic sargo",       "Bodianus pulcher\nCalifornia sheephead",     "Chromis punctipinnis\nblacksmith",          
"Girella nigricans\nopaleye",                 "Hypsypops rubicundus\nGaribaldi damselfish",                                   
"Paralabrax clathratus\nkelp bass",           "Paralabrax nebulifer\nbarred sand bass",     "Stereolepis gigas\ngiant sea bass",  "Other"     ))]



############################################################
#Stacked barplots for all sites together, if not in top 5 species, summed into OTHER category
############################################################

#Make color palettes
color_palette <- c("#80B1D3",#barred sand bass,
                   "#CCCC8F",#California sheephead
                   "#8DD3C7",#opaleye
                   "#BEBADA" ,#rock wrasse
                   "#FDB462", #garibaldi
                   "forestgreen" ,#kelp bass
                   "#FCCDE5",#senorita
                   "#FB8072" ,#kelp perch
                   "#BC80BD" ,#bluebanded goby
                   "#CCEBC5" ,#blacksmith
                    "#D9D9D9")
color_palette_fishbiomass <- c("deepskyblue2",#sargo
                               "#FDB462" ,#garibaldi
                                "#E16A86",#giant sea bass
                               "#80B1D3",#barred sand bass
                               "#8DD3C7" ,#opaleye
                               "#CCCC8F",#california sheephead
                               "forestgreen",#kelp bass
                               "#CCEBC5",#blacksmith
                               "#D9D9D9" #other
                               )


color_palette_fishall <- c("#CCCC8F",#california sheephead
                           "#FB8072",#kelp perch
                           "#CCEBC5",#blacksmith
                           "#8DD3C7" ,#opaleye
                           "#BEBADA",#rock wrasse
                           "#FDB462" ,#garibaldi
                           "#BC80BD" ,#bluebanded goby
                           "#FCCDE5" ,#senorita
                           "forestgreen",#kelp bass
                           "#80B1D3",#barred sand bass
                           "deepskyblue2",#xantic sargo
                           "#E16A86",#giant sea bass
                           "#D9D9D9")#other

color_palette_macroinvert <- c("#8DD3C7", "#CCCC8F", "#BEBADA" ,"#FB8072", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#FCCDE5" ,"#BC80BD" , "#D9D9D9")
color_palette_kelp <- c("#8DD3C7", "#CCCC8F", "#BEBADA" ,"#FB8072", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#FCCDE5" , "#D9D9D9")

#add new rank column
dat_fishdensity_averages_deep_ar.u[,rank := ifelse(full_label == "Other",400,frank(summed_mean_depthzone_density_m2))]

#fish density
fish_density_top5_stacked <- ggplot(dat_fishdensity_averages_deep_ar.u, aes(fill=reorder(full_label,rank), y=(summed_mean_depthzone_density_m2*100), x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Average density per 100m" ^2)), fill = "Fish species") +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","AR")) +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5))

ggsave(fish_density_top5_stacked, path = file.path("figures"), filename = "fish_density_top5_stacked.jpg", height = 4, width = 6, units = "in")

#add new rank column
dat_fishbiomass_averages_deep_ar.u[,summed_mean_depthzone_biomass_kg_100m := (summed_mean_depthzone_biomass_m2*100/1000)]
dat_fishbiomass_averages_deep_ar.u[,rank := ifelse(full_label == "Other",400,frank(summed_mean_depthzone_biomass_m2))]

#fish biomass
fish_biomass_top5_stacked <- ggplot(dat_fishbiomass_averages_deep_ar.u, aes(fill=reorder(full_label,rank), y=(summed_mean_depthzone_biomass_m2*100/1000), x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","AR")) +
  labs(x = "Depth zone", y = expression(paste("Average biomass in kg per 100m" ^2)), fill = "Fish species") +
  scale_fill_manual(values = color_palette_fishbiomass) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5))

ggsave(fish_biomass_top5_stacked, path = file.path("figures"), filename = "fish_biomass_top5_stacked.jpg", height = 4, width = 6, units = "in")

#merge fish plots
#dummy data.table
all_spp <-c("Bodianus pulcher\nCalifornia sheephead",     "Brachyistius frenatus\nkelp perch"    ,      "Chromis punctipinnis\nblacksmith",          
"Girella nigricans\nopaleye",                 "Halichoeres semicinctus\nrock wrasse"  ,     "Hypsypops rubicundus\nGaribaldi damselfish",
"Lythrypnus dalli\nbluebanded goby",          "Oxyjulis californica\nsenorita"        ,     "Paralabrax clathratus\nkelp bass"    ,      
"Paralabrax nebulifer\nbarred sand bass",  "Anisotremus davidsonii\nxantic sargo"     , 
"Stereolepis gigas\ngiant sea bass"  , "Other")

dummy_fish_dt <- data.table(`Fish species` = factor(all_spp, levels = c("Bodianus pulcher\nCalifornia sheephead",     "Brachyistius frenatus\nkelp perch"    ,      "Chromis punctipinnis\nblacksmith",          
                                                                        "Girella nigricans\nopaleye",                 "Halichoeres semicinctus\nrock wrasse"  ,     "Hypsypops rubicundus\nGaribaldi damselfish",
                                                                        "Lythrypnus dalli\nbluebanded goby",          "Oxyjulis californica\nsenorita"        ,     "Paralabrax clathratus\nkelp bass"    ,      
                                                                        "Paralabrax nebulifer\nbarred sand bass",  "Anisotremus davidsonii\nxantic sargo"     , 
                                                                        "Stereolepis gigas\ngiant sea bass"  , "Other")), num = rep(1,13))

fish_leg <- get_legend(ggplot(dummy_fish_dt) +
                         geom_col(aes(x = `Fish species`, y = num, fill = `Fish species`)) +
                         scale_fill_manual(values = color_palette_fishall) +
                         theme_classic() +
                         theme(legend.position = "top",legend.direction = "horizontal"))
#merge w/o legends
fish_abundance_top5_stacked <- plot_grid(fish_density_top5_stacked + theme(legend.position = "none"), fish_biomass_top5_stacked + theme(legend.position = "none"), ncol = 2, labels = c("a.","b."))

#merge with legend
fish_abundance_top5_stacked.l <- plot_grid(fish_leg,fish_abundance_top5_stacked,ncol = 1, rel_heights = c(1,10))

ggsave(fish_abundance_top5_stacked.l, path = file.path("figures"), filename = "fish_abundance_top5_stacked.l.jpg", height = 8.5, width = 10, units = "in")


#Surprised giant sea bass shows up as top! look into this quickly
ggplot(dat_fish_site_averages[taxa == "Stereolepis gigas"]) +
  geom_boxplot(aes(x = DepthZone, y = mean_density_m2)) +
  theme_classic()

ggplot(dat_fish_site_averages[taxa == "Stereolepis gigas"]) +
  geom_boxplot(aes(x = DepthZone, y = mean_wt_density_g_m2)) +
  theme_classic()
  

#macroinvert density
macroinvert_density_top5_stacked <- ggplot(dat_macroinvertdensity_averages_deep_ar.u, aes(fill=full_label, y=summed_mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Average density per 100m" ^2)), fill = "Macroinvertebrate species") +
  scale_fill_manual(values = color_palette_macroinvert) +
    scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","AR")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5))

ggsave(macroinvert_density_top5_stacked, path = file.path("figures"), filename = "macroinvert_density_top5_stacked.jpg", height = 4, width = 6, units = "in")

#kelp density
kelp_density_top5_stacked <- ggplot(dat_kelp_averages_deep_ar.u, aes(fill=full_label, y=summed_mean_depthzone_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Average density per 100m" ^2)), fill = "Macroalgae species") +
  scale_fill_manual(values = color_palette_kelp) +
  scale_x_discrete(labels = c("Inner","Middle","Outer","Deep","AR")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold")) +
  guides(fill = guide_legend(keywidth = 1.5, keyheight = 1.5))

ggsave(kelp_density_top5_stacked, path = file.path("figures"), filename = "kelp_density_top5_stacked.jpg", height = 4, width = 6, units = "in")


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
#Now, average biomass and density by Island vs ARM vs Coast
############################################################

#New Column Identifying ARM vs Island vs Natural Coast
dat_fish_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]
dat_macroinvert_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]
dat_kelp_site_averages[,type := ifelse(DepthZone == "ARM","ARM",ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Natural mainland"))]

#average by type and depth zone
dat_fish_averages_sitetype <- dat_fish_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2),
                                                        mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                                     .(taxa, common_name_final, DepthZone, type)] 
dat_macroinvert_averages_sitetype <- dat_macroinvert_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                                   .(taxa, common_name_final, DepthZone, type)]  
dat_kelp_averages_sitetype <- dat_kelp_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                     .(taxa, common_name_final, DepthZone, type)]  

#Identify top 5 species, sum other into 'other' category
###FISH DENSITY#####

dat_fish_averages_sitetype[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, common_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, summed_mean_depthzone_sitetype_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, type, Species_top5, common_top5)]
dat_fishdensity_averages_sitetype.u <- unique(dat_fish_averages_sitetype[,.(Species_top5, type, DepthZone, summed_mean_depthzone_sitetype_density_m2, common_top5)])
dat_fishdensity_averages_sitetype.u[,full_label := ifelse(Species_top5 == "Other","Other",paste0(Species_top5,"\n", common_top5))]

#add colors and manually set factor order for plotting
dat_fishdensity_averages_sitetype.u[, full_label := factor(full_label, levels = c( "Bodianus pulcher\nCalifornia sheephead"  ,   "Brachyistius frenatus\nkelp perch"    ,      "Chromis punctipinnis\nblacksmith"    ,      
                                                            "Girella nigricans\nopaleye"             ,    "Halichoeres semicinctus\nrock wrasse"    ,   "Hypsypops rubicundus\nGaribaldi damselfish",
                                                            "Lythrypnus dalli\nbluebanded goby"     ,     "Medialuna californiensis\nhalfmoon"     ,       
                                                            "Oxyjulis californica\nsenorita"        ,     "Paralabrax clathratus\nkelp bass"       ,    "Paralabrax nebulifer\nbarred sand bass"  ,  
                                                            "Rhinogobiops nicholsii\nblackeye goby" , "Other"                                 ))]

####KELP######
dat_kelp_averages_sitetype[, Species_top5_kelp := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(DepthZone,type)]
dat_kelp_averages_sitetype[, common_top5_kelp := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(DepthZone,type)]
dat_kelp_averages_sitetype[, summed_mean_depthzone_sitetype_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, type, Species_top5_kelp, common_top5_kelp)]
dat_kelpdensity_averages_sitetype.u <- unique(dat_kelp_averages_sitetype[,.(Species_top5_kelp, type, DepthZone, summed_mean_depthzone_sitetype_density_m2, common_top5_kelp)])
dat_kelpdensity_averages_sitetype.u[,full_label := ifelse(Species_top5_kelp == "Other","Other",paste0(Species_top5_kelp,"\n", common_top5_kelp))]

#add colors and manually set factor order for plotting
dat_kelpdensity_averages_sitetype.u[, full_label := factor(full_label, levels = c(
  "Agarum fimbriatum\nNA","Egregia menziesii\nfeather boa kelp"  ,  "Eisenia arborea\nsouthern sea palm"  ,"Laminaria farlowii\ngolden kombu",
 "Macrocystis pyrifera\ngiant kelp" , "Pelagophycus porra\n"  ,  "Pterygophora californica\nstalked kelp",
 "Sargassum horneri\nS. horneri" ,"Sargassum palmeri\nS. palmeri"  ,  "Sargassum sp\nNA",  "Stephanocystis spp.\nNA"  ,"Other" ), 
 labels = c(
   "Agarum fimbriatum","Egregia menziesii\nfeather boa kelp"  ,  "Eisenia arborea\nsouthern sea palm"  ,"Laminaria farlowii\ngolden kombu",
   "Macrocystis pyrifera\ngiant kelp" , "Pelagophycus porra\nelk Kelp"  ,  "Pterygophora californica\nstalked kelp",
   "Sargassum horneri" ,"Sargassum palmeri"  ,  "Sargassum spp.",  "Stephanocystis spp."  ,"Other" ))]

#####MACRO######
dat_macroinvert_averages_sitetype[, Species_top5_macroinvert := ifelse(frank(-mean_depthzone_density_m2)<=5,taxa,"Other"),.(DepthZone,type)]
dat_macroinvert_averages_sitetype[, common_top5_macroinvert := ifelse(frank(-mean_depthzone_density_m2)<=5,common_name_final,"Other"),.(DepthZone,type)]
dat_macroinvert_averages_sitetype[, summed_mean_depthzone_sitetype_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, type, Species_top5_macroinvert, common_top5_macroinvert)]
dat_macroinvertdensity_averages_sitetype.u <- unique(dat_macroinvert_averages_sitetype[,.(Species_top5_macroinvert, type, DepthZone, summed_mean_depthzone_sitetype_density_m2, common_top5_macroinvert)])
dat_macroinvertdensity_averages_sitetype.u[,full_label := ifelse(Species_top5_macroinvert == "Other","Other",paste0(Species_top5_macroinvert,"\n", common_top5_macroinvert))]

#add colors and manually set factor order for plotting
dat_macroinvertdensity_averages_sitetype.u[, full_label := factor(full_label, levels = c("Anthopleura sola\nstarburst anenome",              "Centrostephanus coronatus\ncrowned urchin" ,       "Haliotis fulgens\ngreen abalone" ,                
                                                                                         "Kelletia kelletii\nKellet's whelk",                "Leptogorgia chilensis\nred gorgonian"  ,           "Megastraea undosa\nwavy turban snail" ,           
                                                                                         "Mesocentrotus franciscanus\nred urchin",           "Muricea californica\ngolden gregorian"   ,         "Muricea fruticosa\nbrown gregorian" ,             
                                                                                         "Patiria miniata\nbat star"                   ,     "Strongylocentrotus purpuratus\npurple sea urchin",
                                                                                         "Styela montereyensis\nstalked tunicate","Other"     ))]


######FISHBIOMASS######
dat_fish_averages_sitetype[, Species_top5_fishbiomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,taxa,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, common_top5_fishbiomass := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,common_name_final,"Other"),.(DepthZone,type)]
dat_fish_averages_sitetype[, summed_mean_depthzone_sitetype_biomass_m2 := sum(mean_depthzone_wt_density_g_m2), .(DepthZone, type, Species_top5_fishbiomass, common_top5_fishbiomass)]
dat_fishbiomass_averages_sitetype.u <- unique(dat_fish_averages_sitetype[,.(Species_top5_fishbiomass, type, DepthZone, summed_mean_depthzone_sitetype_biomass_m2, common_top5_fishbiomass)])
dat_fishbiomass_averages_sitetype.u[,full_label := ifelse(Species_top5_fishbiomass == "Other","Other",paste0(Species_top5_fishbiomass,"\n", common_top5_fishbiomass))]

#manually set factor order for plotting
dat_fishbiomass_averages_sitetype.u[, full_label := factor(full_label, levels = c(
  "Anisotremus davidsonii\nxantic sargo"  ,     "Bodianus pulcher\nCalifornia sheephead"  ,   "Chromis punctipinnis\nblacksmith"   ,        
"Girella nigricans\nopaleye"              ,   "Hermosilla azurea\nNA"                    ,  "Hypsypops rubicundus\nGaribaldi damselfish",
"Medialuna californiensis\nhalfmoon"      ,  "Paralabrax clathratus\nkelp bass" ,         
"Paralabrax nebulifer\nbarred sand bass"  ,   "Sebastes serranoides\nolive rockfish"     ,  "Stereolepis gigas\ngiant sea bass","Other"),
labels = c(
  "Anisotremus davidsonii\nxantic sargo"  ,     "Bodianus pulcher\nCalifornia sheephead"  ,   "Chromis punctipinnis\nblacksmith"   ,        
  "Girella nigricans\nopaleye"              ,   "Hermosilla azurea\nzebra perch"                    ,  "Hypsypops rubicundus\nGaribaldi damselfish",
  "Medialuna californiensis\nhalfmoon"      ,  "Paralabrax clathratus\nkelp bass" ,         
  "Paralabrax nebulifer\nbarred sand bass"  ,   "Sebastes serranoides\nolive rockfish"     ,  "Stereolepis gigas\ngiant sea bass","Other"))]


############################################################
#Stacked barplots for average biomass and density by site type and depth, only top 5 species
############################################################
color_palette <- c("#8DD3C7", "#CCCC8F", "#BEBADA" ,"#FB8072", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#96A1FF","#FCCDE5" ,"#BC80BD" ,
                   #"#CCEBC5" ,#barred sand bass, excluded when ARM excluded
                   "#6FB6FF", "#D9D9D9")

color_palette_fishbiomass <- c("#E16A86","#8DD3C7", "#BEBADA" ,"#FB8072","#FCCDE5","#FDB462","#96A1FF","#BC80BD" ,"#CCEBC5" ,"darksalmon", "#00ABB4", "#D9D9D9")

color_palette_kelp <- c("#8DD3C7","#6FB6FF", "#CCCC8F", "#BEBADA" ,"#FB8072","#BC80BD", "#80B1D3" ,"#FDB462" ,"#B3DE69" ,"#E16A86","#FCCDE5", "#D9D9D9")

color_palette_macroinvert <- c("#8DD3C7", "#CCCC8F","#CCEBC5","#E16A86", "#BEBADA" ,"#FB8072", "#80B1D3" ,"#FDB462" ,"#6FB6FF","#B3DE69" ,"#FCCDE5" , "#D9D9D9")
   

#fish density
fish_density_top5_sitetype_stacked <- ggplot(dat_fishdensity_averages_sitetype.u[type != "ARM"], aes(fill=full_label, y=summed_mean_depthzone_sitetype_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Average density per 100m" ^2)), fill = "Fish\nspecies") +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~type, scales = "free_x", 
             #space = "free", #if using artificial reefs, keep this in and use facet_wrap instead
             ncol = 1) +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))

ggsave(fish_density_top5_sitetype_stacked, path = file.path("figures"), filename = "fish_density_top5_sitetype_stacked.jpg", height = 10, width = 7.5, units = "in")

#fish biomass
fish_biomass_top5_sitetype_stacked <- ggplot(dat_fishbiomass_averages_sitetype.u[type != "ARM"], aes(fill=full_label, y=summed_mean_depthzone_sitetype_biomass_m2*100/1000, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = bquote("Average biomass (kg per 100m"^2*")"), fill = "Fish\nspecies") +
  scale_fill_manual(values = color_palette_fishbiomass) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~type, scales = "free_x", 
             #space = "free", #if using artificial reefs, keep this in and use facet_wrap instead
             ncol = 1) +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))

ggsave(fish_biomass_top5_sitetype_stacked, path = file.path("figures"), filename = "fish_biomass_top5_sitetype_stacked.jpg", height = 10, width = 7.5, units = "in")


#macroinvert density
macroinvert_density_top5_sitetype_stacked <- ggplot(dat_macroinvertdensity_averages_sitetype.u[type != "ARM"], aes(fill=full_label, y=summed_mean_depthzone_sitetype_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Average density per 100m" ^2)), fill = "Macroinvertebrate\nspecies") +
  scale_fill_manual(values = color_palette_macroinvert) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~type, scales = "free_x", 
             #space = "free", #if using artificial reefs, keep this in and use facet_wrap instead
             ncol = 1) +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))


ggsave(macroinvert_density_top5_sitetype_stacked, path = file.path("figures"), filename = "macroinvert_top5_sitetype_stacked.jpg", height = 10, width = 7.5, units = "in")

#kelp density
kelp_density_top5_sitetype_stacked <- ggplot(dat_kelpdensity_averages_sitetype.u[type != "ARM"], aes(fill=full_label, y=summed_mean_depthzone_sitetype_density_m2*100, x=DepthZone)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Depth zone", y = expression(paste("Average density per 100m" ^2)), fill = "Macroalgae\nspecies") +
  scale_fill_manual(values = color_palette_kelp) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~type, scales = "free_x", 
             #space = "free", #if using artificial reefs, keep this in and use facet_wrap instead
             ncol = 1) +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"))

ggsave(kelp_density_top5_sitetype_stacked, path = file.path("figures"), filename = "kelp_density_top5_sitetype_stacked.jpg", height = 10, width = 7.5, units = "in")

############################################################
#Now, average biomass and density by Region
############################################################


#average by type and depth zone
dat_fish_averages_region <- dat_fish_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2),
                                                        mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
                                                     .(Species, DepthZone, Region)] 
dat_macroinvert_averages_region <- dat_macroinvert_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                                   .(BenthicReefSpecies, DepthZone, Region)]  
dat_kelp_averages_region <- dat_kelp_site_averages[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
                                                     .(BenthicReefSpecies, DepthZone, Region)]  

#Identify top 5 species, sum other into 'other' category
###FISH DENSITY#####

dat_fish_averages_region[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,Species,"Other"),.(DepthZone,Region)]
dat_fish_averages_region[, summed_mean_depthzone_region_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Region, Species_top5)]
dat_fishdensity_averages_region.u <- unique(dat_fish_averages_region[,.(Species_top5, Region, DepthZone, summed_mean_depthzone_region_density_m2)])

#add colors and manually set factor order for plotting
dat_fishdensity_averages_region.u <- spp_fish_key[dat_fishdensity_averages_region.u, on = c("Species"="Species_top5")]

####KELP######
dat_kelp_averages_region[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,BenthicReefSpecies,"Other"),.(DepthZone,Region)]
dat_kelp_averages_region[, summed_mean_depthzone_region_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Region, Species_top5)]
dat_kelpdensity_averages_region.u <- unique(dat_kelp_averages_region[,.(Species_top5, Region, DepthZone, summed_mean_depthzone_region_density_m2)])

#add colors and manually set factor order for plotting
dat_kelpdensity_averages_region.u <- spp_kelp_key[dat_kelpdensity_averages_region.u, on = c("Species"="Species_top5")]

#####MACRO######
dat_macroinvert_averages_region[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,BenthicReefSpecies,"Other"),.(DepthZone,Region)]
dat_macroinvert_averages_region[, summed_mean_depthzone_region_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Region, Species_top5)]
dat_macroinvertdensity_averages_region.u <- unique(dat_macroinvert_averages_region[,.(Species_top5, Region, DepthZone, summed_mean_depthzone_region_density_m2)])

#add colors and manually set factor order for plotting
dat_macroinvertdensity_averages_region.u <- spp_macro_key[dat_macroinvertdensity_averages_region.u, on = c("Species"="Species_top5")]


######FISHBIOMASS######
dat_fish_averages_region[, Species_top5 := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,Species,"Other"),.(DepthZone,Region)]
dat_fish_averages_region[, summed_mean_depthzone_region_wt_density_g_m2 := sum(mean_depthzone_wt_density_g_m2), .(DepthZone, Region, Species_top5)]
dat_fishbiomass_averages_region.u <- unique(dat_fish_averages_region[,.(Species_top5, Region, DepthZone, summed_mean_depthzone_region_wt_density_g_m2)])

#add colors and manually set factor order for plotting
dat_fishbiomass_averages_region.u <- spp_fish_key[dat_fishbiomass_averages_region.u, on = c("Species"="Species_top5")]


#Below is regional illustrations, and zooming in for OSM poster, ignore for now unless we want to revive regional analyses
# ############################################################
# #Stacked barplots for average biomass and density by region and depth, only top 5 species
# ############################################################
# 
# #Add site types to help with facet plot
# #New Column Identifying ARM vs Island vs Natural Coast
# dat_fishdensity_averages_region.u[,type := ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")]
# dat_fishbiomass_averages_region.u[,type := ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")]
# dat_macroinvertdensity_averages_region.u[,type := ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")]
# dat_macroinvertdensity_averages_region.u[,type := ifelse(Region %in% c("Santa Catalina Island","Santa Barbara Island","San Clemente Island"),"Island","Mainland")]
# 
# #fish density
# fish_density_top5_region_stacked <- ggplot(dat_fishdensity_averages_region.u, aes(fill=Species_label, y=summed_mean_depthzone_region_density_m2*100, x=DepthZone)) + 
#   geom_bar(position="stack", stat="identity") +
#   labs(x = "Depth zone", y = "Average density (per 100m^2)", fill = "Species") +
#   scale_fill_manual(values = unique(dat_fishdensity_averages_region.u[order(Species_label)]$Species_color)) +
#   facet_grid(~Region, scales = "free_x", space = "free") +
#   theme_classic()
# 
# ggsave(fish_density_top5_region_stacked, path = file.path("figures"), filename = "fish_density_top5_region_stacked.jpg", height = 5, width = 17, units = "in")
# 
# #fish biomass
# fish_biomass_top5_region_stacked <- ggplot(dat_fishbiomass_averages_region.u, aes(fill=Species_label, y=summed_mean_depthzone_region_wt_density_g_m2*100/1000, x=DepthZone)) + 
#   geom_bar(position="stack", stat="identity") +
#   labs(x = "Depth zone", y = "Average biomass (kg per 100m^2)", fill = "Species") +
#   scale_fill_manual(values = unique(dat_fishbiomass_averages_region.u[order(Species_label)]$Species_color)) +
#   facet_grid(~Region, scales = "free_x", space = "free") +
#   theme_classic()
# 
# ggsave(fish_biomass_top5_region_stacked, path = file.path("figures"), filename = "fish_biomass_top5_region_stacked.jpg", height = 5, width = 17, units = "in")
# 
# #macroinvert density
# macroinvert_density_top5_region_stacked <- ggplot(dat_macroinvertdensity_averages_region.u, aes(fill=Species_label, y=summed_mean_depthzone_region_density_m2*100, x=DepthZone)) + 
#   geom_bar(position="stack", stat="identity") +
#   labs(x = "Depth zone", y = "Average density (per 100m^2)", fill = "Species") +
#   scale_fill_manual(values = unique(dat_macroinvertdensity_averages_region.u[order(Species_label)]$Species_color)) +
#   facet_grid(~Region, scales = "free_x", space = "free") +
#   theme_classic()
# 
# 
# ggsave(macroinvert_density_top5_region_stacked, path = file.path("figures"), filename = "macroinvert_top5_region_stacked.jpg", height = 5, width = 17, units = "in")
# 
# #kelp density
# kelp_density_top5_region_stacked <- ggplot(dat_kelpdensity_averages_region.u, aes(fill=Species_label, y=summed_mean_depthzone_region_density_m2*100, x=DepthZone)) + 
#   geom_bar(position="stack", stat="identity") +
#   labs(x = "Depth zone", y = "Average density (per 100m^2)", fill = "Species") +
#   scale_fill_manual(values = unique(dat_kelpdensity_averages_region.u[order(Species_label)]$Species_color)) +
#   facet_grid(~Region, scales = "free_x", space = "free") +
#   theme_classic()
# 
# ggsave(kelp_density_top5_region_stacked, path = file.path("figures"), filename = "kelp_density_top5_region_stacked.jpg", height = 5, width = 17, units = "in")
# 
# ##################################################################################
# #Visual summaries for OSM poster
# ##################################################################################
#           
#           ########################
#           ##Load data
#           ########################
#           dat_event_OSM.r <- readRDS(file.path("data","processed_crane", "dat_event_OSM.r.rds"))
#           dat_fish_site_averages_OSM <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages_OSM.rds"))
#           dat_macroinvert_site_averages_OSM <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages_OSM.rds"))
#           dat_kelp_site_averages_OSM <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages_OSM.rds"))
#           
#           ########################
#           ##Deep only
#           ########################
#           dat_event_OSM.r.deep_ar <- dat_event_OSM.r[DepthZone %in% c("ARM","Deep")]
#           dat_fish_site_averages_OSM.deep_ar <- dat_fish_site_averages_OSM[DepthZone %in% c("ARM","Deep")]
#           dat_macroinvert_site_averages_OSM.deep_ar <- dat_macroinvert_site_averages_OSM[DepthZone %in% c("ARM","Deep")]
#           dat_kelp_site_averages_OSM.deep_ar <- dat_kelp_site_averages_OSM[DepthZone %in% c("ARM","Deep")]
#           
#           #######################
#           ##Counts of each site type
#           ######################
#           count_natural <- length(unique(dat_event_OSM.r.deep_ar[DepthZone == "Deep",Site])) #27
#           count_artificial <- length(unique(dat_event_OSM.r.deep_ar[DepthZone == "ARM",Site])) #28
#           
#           ########################
#           ##Averaged across all sites, top species per depth zone
#           ########################
#           
#           dat_fish_averages_deep_ar <- dat_fish_site_averages_OSM.deep_ar[,.(mean_depthzone_density_m2 = mean(mean_density_m2),
#                                                          mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
#                                                       .(Species, DepthZone)] 
#           dat_macroinvert_averages_deep_ar <- dat_macroinvert_site_averages_OSM.deep_ar[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
#                                                                     .(BenthicReefSpecies, DepthZone)]  
#           dat_kelp_averages_deep_ar <- dat_kelp_site_averages_OSM.deep_ar[,.(mean_depthzone_density_m2 = mean(mean_density_m2)),
#                                                       .(BenthicReefSpecies, DepthZone)]  
#           
#             ###################
#             #Same, but for SoS only
#             ###################
#           dat_fish_averages_SoS <- dat_fish_site_averages_OSM.deep_ar[Site == "Star of Scotland",.(mean_depthzone_density_m2 = mean(mean_density_m2),
#                                                                              mean_depthzone_wt_density_g_m2=mean(mean_wt_density_g_m2)),
#                                                                           .(Species, DepthZone)][,DepthZone := "Star of Scotland"] 
#           dat_macroinvert_averages_SoS <- dat_macroinvert_site_averages_OSM.deep_ar[Site == "Star of Scotland",.(mean_depthzone_density_m2 = mean(mean_density_m2)),
#                                                                                         .(BenthicReefSpecies, DepthZone)][,DepthZone := "Star of Scotland"]  
#           dat_kelp_averages_SoS <- dat_kelp_site_averages_OSM.deep_ar[Site == "Star of Scotland",.(mean_depthzone_density_m2 = mean(mean_density_m2)),
#                                                                           .(BenthicReefSpecies, DepthZone)][,DepthZone := "Star of Scotland"]  
#           
#           ######################
#           #Row bind to add Star of Scotland in
#           ######################
#           dat_fish_averages_deep_ar <- rbind(dat_fish_averages_deep_ar, dat_fish_averages_SoS)
#           dat_macroinvert_averages_deep_ar <- rbind(dat_macroinvert_averages_deep_ar, dat_macroinvert_averages_SoS)
#           dat_kelp_averages_deep_ar <- rbind(dat_kelp_averages_deep_ar, dat_kelp_averages_SoS)
#           
#           #Add column for whether species is in top 5 species (by density) for that zone
#           #fish density
#           dat_fish_averages_deep_ar[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,Species,"Other"), .(DepthZone)]
#           dat_fish_averages_deep_ar[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
#           dat_fishdensity_averages_deep_ar.u <- unique(dat_fish_averages_deep_ar[,.(Species_top5, DepthZone, summed_mean_depthzone_density_m2)])
#           
#           #manually set factor order for plotting
#           dat_fishdensity_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Brachyistius frenatus", "Chromis punctipinnis", "Girella nigricans", "Halichoeres semicinctus", "Hypsypops rubicundus",
#                                                                                           "Lythrypnus dalli", "Oxyjulis californica", "Paralabrax clathratus", "Paralabrax nebulifer", "Semicossyphus pulcher",
#                                                                                           "Other"),
#                                                                  labels = c("Brachyistius frenatus\n(kelp perch)", "Chromis punctipinnis\n(blacksmith damselfish)", "Girella nigricans\n(opaleye)", "Halichoeres semicinctus\n(rock wrasse)", "Hypsypops rubicundus\n(garibaldi)",
#                                                                             "Lythrypnus dalli\n(blue-banded goby)", "Oxyjulis californica\n(señorita wrasse)", "Paralabrax clathratus\n(kelp bass)", "Paralabrax nebulifer\n(barred sand bass)", "Semicossyphus pulcher\n(sheephead)",
#                                                                             "Other"))]
#           
#           #site type column
#           dat_fishdensity_averages_deep_ar.u[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]
#           
#           #reorder these factors
#           dat_fishdensity_averages_deep_ar.u[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
#                                                                     labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]
#           
#           #macroinvert density
#           dat_macroinvert_averages_deep_ar[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,BenthicReefSpecies,"Other"), .(DepthZone)]
#           dat_macroinvert_averages_deep_ar[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
#           dat_macroinvert_averages_deep_ar.u <- unique(dat_macroinvert_averages_deep_ar[,.(Species_top5, DepthZone, summed_mean_depthzone_density_m2)])
#           
#           dat_macroinvert_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Anthopleura sola", "Apostichopus parvimensis", "Centrostephanus coronatus","Crassadoma gigantea","Kelletia kelletii", "Leptogorgia chilensis",
#                                                                                           "Megastraea undosa", "Mesocentrotus franciscanus", "Muricea californica","Muricea fruticosa","Neobernaya spadicea", "Pachycerianthus fimbriatus",
#                                                                                           "Patiria miniata","Strongylocentrotus purpuratus", "Other"),
#                                                                  labels = c("Anthopleura sola\n(starburst anenome)", "Apostichopus parvimensis\n(warty sea cucumber)", "Centrostephanus coronatus\n(crowned urchin)","Crassadoma gigantea\n(giant rock scallop)","Kelletia kelletii\n(Kellet's whelk)", "Leptogorgia chilensis\n(red gorgonian)",
#                                                                             "Megastraea undosa\n(wavy turban snail)", "Mesocentrotus franciscanus\n(red urchin)", "Muricea californica\n(golden gregorian)","Muricea fruticosa\n(brown gregorian)","Neobernaya spadicea\n(chesnut cowrie)", "Pachycerianthus fimbriatus\n(tube dwelling anenome)",
#                                                                             "Patiria miniata\n(bat star)","Strongylocentrotus purpuratus\n(purple urchin)", "Other"))]
#           #site type column
#           dat_macroinvert_averages_deep_ar.u[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]
#           
#           #reorder these factors
#           dat_macroinvert_averages_deep_ar.u[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
#                                                                     labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]
#           
#           
#           #Kelp density
#           dat_kelp_averages_deep_ar[, Species_top5 := ifelse(frank(-mean_depthzone_density_m2)<=5,BenthicReefSpecies,"Other"), .(DepthZone)]
#           dat_kelp_averages_deep_ar[, summed_mean_depthzone_density_m2 := sum(mean_depthzone_density_m2), .(DepthZone, Species_top5)]
#           dat_kelp_averages_deep_ar.u <- unique(dat_kelp_averages_deep_ar[,.(Species_top5, DepthZone, summed_mean_depthzone_density_m2)])
#           
#           dat_kelp_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Agarum fimbriatum", "Eisenia arborea","Egregia menziesii" , "Laminaria farlowii", "Macrocystis pyrifera","Pelagophycus porra",
#                                                                                    "Pterygophora californica", "Sargassum horneri","Sargassum palmeri", "Stephanocystis spp.",  "Other"))]
#           
#           #site type column
#           dat_kelp_averages_deep_ar.u[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]
#           
#           #reorder these factors
#           dat_kelp_averages_deep_ar.u[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
#                                                                     labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]
#           
#           #Add column for whether species is in top 5 species (by biomass) for that zone, and then sum biomass for all others
#           #fish biomass
#           dat_fish_averages_deep_ar[, Species_top5 := ifelse(frank(-mean_depthzone_wt_density_g_m2)<=5,Species,"Other"), .(DepthZone)]
#           dat_fish_averages_deep_ar[, summed_mean_depthzone_wt_density_g_m2 := sum(mean_depthzone_wt_density_g_m2), .(DepthZone, Species_top5)]
#           dat_fishbiomass_averages_deep_ar.u <- unique(dat_fish_averages_deep_ar[,.(Species_top5, DepthZone, summed_mean_depthzone_wt_density_g_m2)])
#           
#           dat_fishbiomass_averages_deep_ar.u[,Species_forlabel := factor(Species_top5, levels = c("Anisotremus davidsonii","Chromis punctipinnis", "Girella nigricans", "Hypsypops rubicundus",
#                                                                                           "Paralabrax clathratus", "Paralabrax nebulifer","Sebastes serranoides", "Semicossyphus pulcher","Stereolepis gigas",
#                                                                                           "Other"),
#                                                                  labels = c("Anisotremus davidsonii\n(sargo)", "Chromis punctipinnis\n(blacksmith damselfish)", "Girella nigricans\n(opaleye)", "Hypsypops rubicundus\n(garibaldi)",
#                                                                             "Paralabrax clathratus\n(kelp bass)", "Paralabrax nebulifer\n(barred sand bass)","Sebastes serranoides\n(olive rockfish)", "Semicossyphus pulcher\n(sheephead)","Stereolepis gigas\n(giant sea bass)",
#                                                                             "Other"))]
#           
#           #site type column
#           dat_fishbiomass_averages_deep_ar.u[,`Site type` := ifelse(DepthZone=="ARM","Artificial",ifelse(DepthZone=="Deep","Natural","Star of Scotland"))]
#           
#           #reorder these factors
#           dat_fishbiomass_averages_deep_ar.u[,`Site type` := factor(`Site type`, levels = c("Natural","Artificial","Star of Scotland"),
#                                                                     labels = c(paste0("Natural\nn = ",count_natural),paste0("Artificial\nn = ",count_artificial),paste0("Star of Scotland\nn = ","1")))]
#           
#           
#           ############################################################
#           #Stacked barplots for all sites together, if not in top 5 species, summed into OTHER category
#           ############################################################
#           
#           #fish density palette with 9 colors
#           pal_fishdens9 <- c("#CCEBC5",  #"Brachyistius frenatus"
#                               "#80B1D3",  #"Chromis punctipinnis"
#                               "#8ea489",  #"Girella nigricans"
#                            #   "#FCCDE5",  #"Halichoeres semicinctus"
#                            #   "#FDB462",  #"Hypsypops rubicundus"
#                               "lightskyblue",  #"Lythrypnus dalli"
#                               "#FFFFB3",  #"Oxyjulis californica"
#                               "#8DD3C7",  #"Paralabrax clathratus"
#                               "#FFED6F",  #"Paralabrax nebulifer"
#                               "darksalmon",  #"Semicossyphus pulcher"
#                               "#D9D9D9")  #"Other"
#           
#           #fish biomass palette with 7 colors
#           pal_fishbio7 <- c(
#                              #"#CEAD64",#"Anisotremus davidsonii" 
#                              "#80B1D3",#"Chromis punctipinnis"
#                              "#8ea489",#"Girella nigricans" 
#                             # "#FDB462",#"Hypsypops rubicundus"
#                              "#8DD3C7",#"Paralabrax clathratus"
#                              "#FFED6F", #"Paralabrax nebulifer"
#                              #"#9D9E39", #"Sebastes serranoides"
#                              "darksalmon",#"Semicossyphus pulcher"
#                              "#A38389",#"Stereolepis gigas"
#                              "#D9D9D9")
#           
#           #macro density palette with 11 colors
#           pal_macro11 <- c(
#                             #"#97D4BA",#"Anthopleura sola"
#                            #"#B67436",#"Apostichopus parvimensis"
#                            "#3F4965",#"Centrostephanus coronatus"
#                             "#EDC75C", #Crassadoma gigantea ADD 
#                            "#814A23", # "Kelletia kelletii"
#                            "#DA7E80",#"Leptogorgia chilensis"
#                            "#9C8074", #"Megastraea undosa"
#                            "#BA4C61", #"Mesocentrotus franciscanus"
#                            "#BF9D5D",#"Muricea californica"
#                            "#C9664B",#"Muricea fruticosa" Brown grogorian
#                            "#411009", #"Neobernaya spadicea" Chesnut cowrie
#                           # "#D7CDAA",# "Pachycerianthus fimbriatus"
#                           # "#E08454",# "Patiria miniata"
#                            "#BC80BD",#"Strongylocentrotus purpuratus"
#                            "#D9D9D9")#"Other"
# 
#           
#           #kelp density palette with 9 colors
#           pal_kelp8 <- c("#8DD3C7",#"Agarum fimbriatum"
#                         # "#E9D677",#"Eisenia arborea"
#                          "#FCEFDB", # "Egregia menziesii" 
#                          "#BEBADA",#"Laminaria farlowii"
#                          "#FB8072",#"Macrocystis pyrifera"
#                         "#38597A", # "Pelagophycus porra"  
#                          "#3E8AAF",#"Pterygophora californica"
#                         # "turquoise",# "Sargassum horneri"
#                         # "#B3DE69",#"Sargassum palmeri"
#                          "#CDA3AB",#"Stephanocystis spp."
#                          "#D9D9D9") #other
#   
#           
#           #fish density
#           fish_density_top5deep_ar_stacked <- ggplot(dat_fishdensity_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_density_m2*100, x=`Site type`)) + 
#             geom_bar(position="stack", stat="identity") +
#             labs(x = "Site type", y = expression(paste("Average density per 100m" ^2)), fill = "Fish species") +
#             scale_fill_manual(values = pal_fishdens9) +
#             scale_y_continuous(expand = c(0,0)) +
#             theme_classic() +
#             theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))
#           
#           ggsave(fish_density_top5deep_ar_stacked, path = file.path("figures"), filename = "fish_density_top5deep_ar_stacked.jpg", height = 3.5, width = 5, units = "in")
#           
#           #fish biomass
#           fish_biomass_top5deep_ar_stacked <- ggplot(dat_fishbiomass_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_wt_density_g_m2*100/1000, x=`Site type`)) + 
#             geom_bar(position="stack", stat="identity") +
#             labs(x = "Site type", y = expression(paste("Average biomass in kg per 100m" ^2)), fill = "Fish species") +
#             scale_fill_manual(values = pal_fishbio7) +
#             scale_y_continuous(expand = c(0,0)) +
#             theme_classic() +
#             theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))
#           
#           ggsave(fish_biomass_top5deep_ar_stacked, path = file.path("figures"), filename = "fish_biomass_top5deep_ar_stacked.jpg", height = 3, width = 5, units = "in")
#           
#           #macroinvert density
#           macroinvert_density_top5deep_ar_stacked <- ggplot(dat_macroinvert_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_density_m2*100, x=`Site type`)) + 
#             geom_bar(position="stack", stat="identity") +
#             labs(x = "Site type", y = expression(paste("Average density per 100m" ^2)), fill = "Macroinvertebrate species") +
#             scale_fill_manual(values = pal_macro10) +
#             scale_y_continuous(expand = c(0,0)) +
#             theme_classic() +
#             theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))
#           
#           ggsave(macroinvert_density_top5deep_ar_stacked, path = file.path("figures"), filename = "macroinvert_density_top5deep_ar_stacked.jpg", height = 3.5, width = 5, units = "in")
#           
#           #kelp density
#           kelp_density_top5deep_ar_stacked <- ggplot(dat_kelp_averages_deep_ar.u, aes(fill=Species_forlabel, y=summed_mean_depthzone_density_m2*100, x=`Site type`)) + 
#             geom_bar(position="stack", stat="identity") +
#             labs(x = "Site type", y = expression(paste("Average density per 100m" ^2)), fill = "Kelp species") +
#             scale_fill_manual(values = pal_kelp8) +
#             scale_y_continuous(expand = c(0,0)) +
#             theme_classic() +
#             theme(legend.text = element_text(size = 5), legend.title = element_text(size = 6, face = "bold"))
#           
#           ggsave(kelp_density_top5deep_ar_stacked, path = file.path("figures"), filename = "kelp_density_top5deep_ar_stacked.jpg", height = 3.5, width = 5, units = "in")
#           
#           #Merge figures
#       CRANE_deep_artificial_SoS_2022_23 <- plot_grid(fish_density_top5deep_ar_stacked + theme(axis.title.x = element_blank()),
#                     fish_biomass_top5deep_ar_stacked + theme(axis.title.x = element_blank()),
#                     macroinvert_density_top5deep_ar_stacked,
#                     kelp_density_top5deep_ar_stacked + labs(y=" "))
#       
#       ggsave(CRANE_deep_artificial_SoS_2022_23, path = file.path("figures"), filename = "CRANE_deep_artificial_SoS_2022_23.jpg",
#              height = 8, width = 11, units = "in")
#       ggsave(CRANE_deep_artificial_SoS_2022_23, path = file.path("figures"), filename = "CRANE_deep_artificial_SoS_2022_23.pdf",
#              height = 8, width = 11, units = "in")
#       
# 
# 
# ##################################################################################
# #TO DO: average biomass and density for natural reefs only, MPA vs non-MPA (not sure how to determine this right now)
# ##################################################################################
# 
# #PARKING LOT
#       
#       #old coloring and naming
#       
#       
#       ########################
#       ##Species, label, color key
#       ########################
#       spp_fish_key <- data.table(Species = factor(c("Anisotremus davidsonii","Brachyistius frenatus", "Chromis punctipinnis","Embiotoca jacksoni", "Girella nigricans", "Halichoeres semicinctus","Hermosilla azurea", "Hypsypops rubicundus",
#                                                     "Lythrypnus dalli","Medialuna californiensis", "Oxyjulis californica", "Paralabrax clathratus", "Paralabrax nebulifer","Rhacochilus toxotes", "Rhinogobiops nicholsii","Sebastes serranoides", "Semicossyphus pulcher",
#                                                     "Stereolepis gigas", "Other"), 
#                                                   levels = c("Anisotremus davidsonii","Brachyistius frenatus", "Chromis punctipinnis","Embiotoca jacksoni", "Girella nigricans", "Halichoeres semicinctus","Hermosilla azurea", "Hypsypops rubicundus",
#                                                              "Lythrypnus dalli","Medialuna californiensis", "Oxyjulis californica", "Paralabrax clathratus", "Paralabrax nebulifer","Rhacochilus toxotes", "Rhinogobiops nicholsii","Sebastes serranoides", "Semicossyphus pulcher",
#                                                              "Stereolepis gigas", "Other")),
#                                  Species_label = factor(c("Anisotremus davidsonii\n(sargo)","Brachyistius frenatus\n(kelp perch)", "Chromis punctipinnis\n(blacksmith damselfish)","Embiotoca jacksoni\n(black surfperch)", "Girella nigricans\n(opaleye)", "Halichoeres semicinctus\n(rock wrasse)","Hermosilla azurea\n(zebra perch)", "Hypsypops rubicundus\n(garibaldi)",
#                                                           "Lythrypnus dalli\n(blue-banded goby)","Medialuna californiensis\n(halfmoon)", "Oxyjulis californica\n(señorita wrasse)", "Paralabrax clathratus\n(kelp bass)", "Paralabrax nebulifer\n(barred sand bass)","Rhacochilus toxotes\n(rubberlip seaperch)", "Rhinogobiops nicholsii\n(blackeye goby)","Sebastes serranoides\n(olive rockfish)", "Semicossyphus pulcher\n(sheephead)",
#                                                           "Stereolepis gigas\n(giant sea bass)", "Other"), 
#                                                         levels = c("Anisotremus davidsonii\n(sargo)","Brachyistius frenatus\n(kelp perch)", "Chromis punctipinnis\n(blacksmith damselfish)","Embiotoca jacksoni\n(black surfperch)", "Girella nigricans\n(opaleye)", "Halichoeres semicinctus\n(rock wrasse)","Hermosilla azurea\n(zebra perch)", "Hypsypops rubicundus\n(garibaldi)",
#                                                                    "Lythrypnus dalli\n(blue-banded goby)","Medialuna californiensis\n(halfmoon)", "Oxyjulis californica\n(señorita wrasse)", "Paralabrax clathratus\n(kelp bass)", "Paralabrax nebulifer\n(barred sand bass)","Rhacochilus toxotes\n(rubberlip seaperch)", "Rhinogobiops nicholsii\n(blackeye goby)","Sebastes serranoides\n(olive rockfish)", "Semicossyphus pulcher\n(sheephead)",
#                                                                    "Stereolepis gigas\n(giant sea bass)", "Other")),
#                                  Species_common = factor(c("sargo","kelp perch", "blacksmith damselfish","black surfperch", "opaleye", "rock wrasse","zebra perch", "garibaldi",
#                                                            "blue-banded goby","halfmoon", "señorita wrasse", "kelp bass", "barred sand bass","rubberlip seaperch", "blackeye goby","olive rockfish", "sheephead",
#                                                            "giant sea bass", "Other"), 
#                                                          levels = c("sargo","kelp perch", "blacksmith damselfish","black surfperch", "opaleye", "rock wrasse","zebra perch", "garibaldi",
#                                                                     "blue-banded goby","halfmoon", "señorita wrasse", "kelp bass", "barred sand bass","rubberlip seaperch", "blackeye goby","olive rockfish", "sheephead",
#                                                                     "giant sea bass", "Other")),
#                                  Species_color = c("#CEAD64", #"Anisotremus davidsonii" 
#                                                    "#CCEBC5",  #"Brachyistius frenatus"
#                                                    "#80B1D3",  #"Chromis punctipinnis"
#                                                    "purple",#black surfperch
#                                                    "#8ea489",  #"Girella nigricans"
#                                                    "#FCCDE5",  #"Halichoeres semicinctus"
#                                                    "#2C3D42", #Hermosilla azurea
#                                                    "#FDB462",  #"Hypsypops rubicundus"
#                                                    "lightskyblue",  #"Lythrypnus dalli"
#                                                    "#DAB776",#Medialuna californiensis
#                                                    "#FFFFB3",  #"Oxyjulis californica"
#                                                    "#8DD3C7",  #"Paralabrax clathratus"
#                                                    "#FFED6F",  #"Paralabrax nebulifer"
#                                                    "#CDC480", #Rhacochilus toxotes
#                                                    "#DEBB83", #Rhinogobiops nicholsii
#                                                    "#5F5C25", #olive rockfish
#                                                    "darksalmon",  #"Semicossyphus pulcher",
#                                                    "#A38389",#"Stereolepis gigas"
#                                                    "#D9D9D9")  #"Other"
#       )
#       
#       spp_macro_key <- data.table(Species = factor(c("Anthopleura sola", "Apostichopus parvimensis", "Centrostephanus coronatus","Haliotis fulgens", "Kelletia kelletii", "Leptogorgia chilensis",
#                                                      "Megastraea undosa", "Mesocentrotus franciscanus", "Muricea californica","Muricea fruticosa", "Pachycerianthus fimbriatus","Panulirus interruptus",
#                                                      "Patiria miniata","Strongylocentrotus purpuratus", "Other"),
#                                                    levels = c("Anthopleura sola", "Apostichopus parvimensis", "Centrostephanus coronatus","Haliotis fulgens", "Kelletia kelletii", "Leptogorgia chilensis",
#                                                               "Megastraea undosa", "Mesocentrotus franciscanus", "Muricea californica","Muricea fruticosa", "Pachycerianthus fimbriatus","Panulirus interruptus",
#                                                               "Patiria miniata","Strongylocentrotus purpuratus", "Other")),
#                                   Species_label = factor(c("Anthopleura sola\n(starburst anenome)", "Apostichopus parvimensis\n(warty sea cucumber)", "Centrostephanus coronatus\n(crowned urchin)","Haliotis fulgens\nGreen abalone", "Kelletia kelletii\n(Kellet's whelk)", "Leptogorgia chilensis\n(red gorgonian)",
#                                                            "Megastraea undosa\n(wavy turban snail)", "Mesocentrotus franciscanus\n(red urchin)", "Muricea californica\n(golden gregorian)","Muricea fruticosa\n(brown gregorian)", "Pachycerianthus fimbriatus\n(tube dwelling anenome)",
#                                                            "Panulirus interruptus\n(CA spiny lobster)", "Patiria miniata\n(bat star)","Strongylocentrotus purpuratus\n(purple urchin)", "Other"),
#                                                          levels = c("Anthopleura sola\n(starburst anenome)", "Apostichopus parvimensis\n(warty sea cucumber)", "Centrostephanus coronatus\n(crowned urchin)","Haliotis fulgens\nGreen abalone", "Kelletia kelletii\n(Kellet's whelk)", "Leptogorgia chilensis\n(red gorgonian)",
#                                                                     "Megastraea undosa\n(wavy turban snail)", "Mesocentrotus franciscanus\n(red urchin)", "Muricea californica\n(golden gregorian)","Muricea fruticosa\n(brown gregorian)", "Pachycerianthus fimbriatus\n(tube dwelling anenome)",
#                                                                     "Panulirus interruptus\n(CA spiny lobster)",  "Patiria miniata\n(bat star)","Strongylocentrotus purpuratus\n(purple urchin)", "Other")),
#                                   Species_common =factor(c("starburst anenome", "warty sea cucumber", "crowned urchin","Haliotis fulgens", "Kellet's whelk", "red gorgonian",
#                                                            "wavy turban snail", "red urchin", "golden gregorian","brown gregorian", "tube dwelling anenome",
#                                                            "CA spiny lobster", "bat star","purple urchin", "Other"),levels = c("starburst anenome", "warty sea cucumber", "crowned urchin","Haliotis fulgens", "Kellet's whelk", "red gorgonian",
#                                                                                                                                "wavy turban snail", "red urchin", "golden gregorian","brown gregorian", "tube dwelling anenome",
#                                                                                                                                "CA spiny lobster", "bat star","purple urchin", "Other")),
#                                   Species_color = c("#97D4BA", #starburst anenome
#                                                     "#B67436",#warty sea cucumber
#                                                     "#3F4965",#crowned urchin
#                                                     "#C6E4C9",#green abalone
#                                                     "#814A23", #Kellet's whelk
#                                                     "#DA7E80",#red gorgonian
#                                                     "#9C8074",#wavy turban snail
#                                                     "#BA4C61",#red urchin
#                                                     "#BF9D5D",#golden gregorian
#                                                     "#D7CDAA",#brown gregorian
#                                                     "#BD2D26",#tube dwelling anenome
#                                                     "#592C2F",#spiny lobster
#                                                     "#E08454",#bat star
#                                                     "#BC80BD",#purple urchin
#                                                     "#D9D9D9")  #"Other"
#       )
#       
#       spp_kelp_key <- data.table(Species = factor(c("Agarum fimbriatum", "Egregia menziesii","Eisenia arborea", "Laminaria farlowii", "Macrocystis pyrifera",
#                                                     "Pterygophora californica", "Sargassum horneri","Sargassum muticum","Sargassum palmeri","Sargassum sp", "Stephanocystis spp.",  "Other"),
#                                                   levels = c("Agarum fimbriatum","Egregia menziesii", "Eisenia arborea", "Laminaria farlowii", "Macrocystis pyrifera",
#                                                              "Pterygophora californica", "Sargassum horneri","Sargassum muticum","Sargassum palmeri","Sargassum sp", "Stephanocystis spp.",  "Other")),
#                                  Species_label = factor(c("Agarum fimbriatum\nfringed sieve kelp", "Egregia menziesii\nfeather boa kelp","Eisenia arborea\nsouthern sea palm", "Laminaria farlowii\ngolden kombu", "Macrocystis pyrifera\ngiant kelp",
#                                                           "Pterygophora californica\nstalked kelp", "Sargassum horneri","Sargassum muticum","Sargassum palmeri","Sargassum sp", "Stephanocystis spp.\nchainbladder kelp",  "Other"),
#                                                         levels = c("Agarum fimbriatum\nfringed sieve kelp","Egregia menziesii\nfeather boa kelp", "Eisenia arborea\nsouthern sea palm", "Laminaria farlowii\ngolden kombu", "Macrocystis pyrifera\ngiant kelp",
#                                                                    "Pterygophora californica\nstalked kelp", "Sargassum horneri","Sargassum muticum","Sargassum palmeri","Sargassum sp", "Stephanocystis spp.\nchainbladder kelp",  "Other")),
#                                  Species_common = factor(c("fringed sieve kelp","feather boa kelp", "southern sea palm", "golden kombu", "giant kelp",
#                                                            "stalked kelp", "S. horneri","S. muticum","S. palmeri","Sargassum sp", "chainbladder kelp",  "Other"),
#                                                          levels = c("fringed sieve kelp","feather boa kelp", "southern sea palm", "golden kombu", "giant kelp",
#                                                                     "stalked kelp", "S. horneri","S. muticum","S. palmeri","Sargassum sp", "chainbladder kelp",  "Other")),
#                                  Species_color = c("#8DD3C7", "#E9D677", "#BEBADA","#BC80BD", "#FB8072", "#3E8AAF",
#                                                    "turquoise", "#B3DE69", "#CDA3AB","coral","#BA4C61", "#D9D9D9"
#                                  ))
#       
#       #top_spp_key <- rbind(spp_fish_key, spp_kelp_key, spp_macro_key)
#       
#       #saveRDS(top_spp_key, file = file.path("keys","top_spp_key.Rds"))
#       
