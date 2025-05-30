# CREATION DATE 14 Apr 2023

# AUTHOR: kitchel@oxy.edu

# PURPOSE: Links taxonomy, colors, and common names with species

#############################
##Setup
#############################
library(data.table)
library(fishualize)

source(file.path("functions","return_spptaxonomy_function.R"))

########################
##Load data
########################
dat_fish_site_averages <- readRDS(file.path("data","processed_crane", "dat_fish_site_averages.rds"))
dat_macroinvert_site_averages <- readRDS(file.path("data","processed_crane", "dat_macroinvert_site_averages.rds"))
dat_kelp_site_averages <- readRDS(file.path("data","processed_crane", "dat_kelp_site_averages.rds"))

#I will use function, but not all fish named from function, so in these cases tie in VRG key from All Fish file
all_fish_key_VRG <- fread(file.path("keys", "all_fish_key_VRG.csv"))

all_fish_key_VRG <- all_fish_key_VRG[,1:2] #delete family column

#################
#Develop species list
#################
fish <- unique(dat_fish_site_averages$Species)
kelp <- unique(dat_kelp_site_averages[BenthicReefSpecies != "Macrocystis pyrifera stipes",]$BenthicReefSpecies)
macroinvert <- unique(dat_macroinvert_site_averages$BenthicReefSpecies)

all_species <- c(fish, kelp, macroinvert)

#Prep data table
species_key <- data.table(Species = all_species)

#apply function
species_info <- get_taxa(all_species)

#merge two data tables
species_key <- species_info[species_key, on = c("query" = "Species")]

#merge with VRG common names
species_key <- all_fish_key_VRG[species_key, on = c("Taxa" = "query")]

#merge with manual kelp/invert names

spp_top_macro_kelp_key <- data.table(Species = c("Anthopleura sola", "Apostichopus parvimensis", "Centrostephanus coronatus","Haliotis fulgens", "Kelletia kelletii", "Leptogorgia chilensis",
                                                 "Megastraea undosa", "Mesocentrotus franciscanus", "Muricea californica","Muricea fruticosa", "Pachycerianthus fimbriatus","Panulirus interruptus",
                                                 "Patiria miniata","Strongylocentrotus purpuratus",
                                                 "Agarum fimbriatum", "Egregia menziesii","Eisenia arborea", "Laminaria farlowii", "Macrocystis pyrifera",
                                                 "Pterygophora californica", "Sargassum horneri","Sargassum muticum","Sargassum palmeri","Sargassum sp", "Stephanocystis spp.", "Styela montereyensis"),
                                     Species_common =c("starburst anenome", "warty sea cucumber", "crowned urchin","Haliotis fulgens", "Kellet's whelk", "red gorgonian",
                                                       "wavy turban snail", "red urchin", "golden gorgonian","brown gorgonian", "tube dwelling anenome",
                                                       "CA spiny lobster", "bat star","purple urchin",
                                                       "fringed sieve kelp","feather boa kelp", "southern sea palm", "golden kombu", "giant kelp",
                                                       "stalked kelp", "sargassum (horneri)","sargassum (muticum)","sargassum (palmeri)","Sargassum sp", "chainbladder kelp", "stalked tunicate"))

species_key <- spp_top_macro_kelp_key[species_key, on = c("Species" = "Taxa")]


#rename some columns
species_key[,common_name_final := ifelse(!is.na(common_name), common_name, #from taxize function
                                         ifelse(!is.na(CommonName), tolower(CommonName), #from VRG All fish
                                                ifelse(!is.na(Species_common),Species_common, #from manual additions of most common macro and kelp species
                                                NA)))]

species_key <- species_key[,.(worms_id,taxa,common_name_final,kingdom, phylum, class, order, family, genus, rank)]

#In case sargassum fix didn't work
species_key[, common_name_final := gsub("Sargassum horneri", "Sargassum (horneri)", common_name_final)]
species_key[, common_name_final := gsub("Sargassum muticum", "Sargassum (muticum)", common_name_final)]
species_key[, common_name_final := gsub("Sargassum palmeri", "Sargassum (palmeri)", common_name_final)]

fwrite(species_key, file.path("keys","species_key.csv"))

