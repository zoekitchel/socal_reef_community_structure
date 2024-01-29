#notes from talking with Jeremy
#1. Use all VRG data w/ typical DepthZone (Islands: SB, SCL, SCLI, San Nic?, Begg Rock, mainland: Malibu to SD regions)
#2. Assemblages

    # fish - density w/ YOY removed
# 
    # fish - biomass
# 
    # macroinvert (see Zahn 2016) (kelps?)
# 
    # UPC? (or just use as explanatory variables)
# 
    # Rocky reef/kelp forest combined assemblage (fish density + macroinvert + kelp density)
# 
# 3. Explanatory variables
# 
    # Natural vs. AR
# 
    # MPA vs outside MPA
# 
    # DepthZone specific versions of fish biomass paper (+ others) explanatory variables 
                #(remote sensed and in-situ versions of relief and substrate etc)

#Jeremy: AR_NR_Aseemblage_Comparison.R from SSINP_ARs
# Read in reference tables -----------------------------------------
source("~/Dropbox/VRG Files/R Code/DataFiles/READ_Sites_Fish_BRS.r")

# Load libraries------------------------------------------------------------
# done after reference tables since that loads package data.table which causes conflicts with tidyverse
library(tidyverse) 
library(vegan)

#ggvegan is not yet on CRAN, so install from GitHub
#remotes::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(ggrepel)

#Color Palettes Packege
library(paletteer)

# palettes_d_names |> 
#   arrange(desc(length)) |> 
#   print(n = 200)
# #https://github.com/EmilHvitfeldt/r-color-palettes/blob/master/type-sorted-palettes.md

# Read in Data Tables -----------------------------------------
# 
#Processed CRANE Fish 
dat_fish_t <- read_csv("Data Tables/dat_fish_t_03_Jan_2024.csv")

#Processed CRANE Event 
dat_event <- read_csv("Data Tables/dat_event_03_Jan_2024.csv")



# Add Level Orders to Region, Site, DepthZone -----------------------------
source("~/Dropbox/VRG Files/R Code/Integrated Dive General/CRANE_level_orders.R")

dat_fish_t <- CRANE_level_orders(dat_fish_t, AR_Region = T, AR_Complex = T) |> 
  droplevels()

glimpse(dat_fish_t)

## Check if any ARs have AR Complex assigned and Region recoded properly
dat_fish_t |>
  filter(DepthZone == "ARM") |>
  distinct(AR_Complex, Site, Region) |> 
  arrange(AR_Complex, Site) |> 
  print(n = Inf)

# Filtering ---------------------------------------------------------------

# Filter out PVR Pre-Construction & sampling in 2020 right after construction
dat_fish_t <- dat_fish_t |>
  filter(!(str_detect(Site, "PVR") & SampleYear < 2021 )) |> 
  droplevels()

dat_fish_t <- dat_fish_t |>
  filter(!Site %in% c("Leucadia", "Old 18th")) |> #some depth zones were way out there in nMDS (probably very few species?)
  droplevels() 

# Filter out species with very few observations???


#TO DO? Add PVR replicate descriptions (reef vs. sand) - ask Jonathon for these, filter out sand/halo classified transects? 

# Summarize Sampling Effort -----------------------------------------------

# count of fish transects by site, depth zone
t_count_s_dz <- dat_fish_t |> 
  group_by(AR_Complex, Site, DepthZone) |> 
  summarise(n_transects = n_distinct(Transect)) |> 
  pivot_wider(names_from = DepthZone, values_from = n_transects)

# Only ARs, across years by site, depth zone
t_count_ARs_yr <- dat_fish_t |>
  filter(DepthZone == "ARM") |> 
  group_by(SampleYear, AR_Complex, Site) |> #all DepthZone == ARM, so don't need to include
  summarise(n_transects = n_distinct(Transect)) |> 
  pivot_wider(names_from = SampleYear, values_from = n_transects)

# count of fish transects by AR complex, yr
t_count_ARC_yr <- dat_fish_t |>
  filter(DepthZone == "ARM") |> 
  group_by(AR_Complex, SampleYear) |> 
  summarise(n_transects = n_distinct(Transect)) |> 
  pivot_wider(names_from = SampleYear, values_from = n_transects)

# Calculate means across transects ----------------------------------------

#Calculate means for each year (across transects), then over years by Site + DepthZone, for density and biomass
# + multiply densities by 100 to get per 100m2
dat_fish_s_dz_sp <- dat_fish_t |> 
  #across transects for each year
  group_by(Species, SampleYear, Region, AR_Complex, Site, DepthZone) |> 
  summarise(mean_density_100m2 = mean(density_m2*100, na.rm = TRUE),
            mean_biomass_g_100m2 = mean(wt_density_g_m2*100, na.rm = TRUE)) |> 
  ungroup() |>
  #across years for each site by depthzone
  group_by(Species, Region, AR_Complex, Site, DepthZone) |> 
  summarise(mean_density_100m2 = mean(mean_density_100m2, na.rm = TRUE),
            mean_biomass_g_100m2 = mean(mean_biomass_g_100m2, na.rm = TRUE)) |> 
  ungroup() 


# count of sites, by AR Complex and depth zone
s_count_ARC_dz <- dat_fish_s_dz_sp |> 
  group_by(AR_Complex, DepthZone) |> 
  summarise(n_sites = n_distinct(Site)) |> 
  pivot_wider(names_from = DepthZone, values_from = n_sites)

# Stacked Bar Plots -------------------------------------------------------

# Find overall proportions of species (across all DepthZones (Sites) in data set)
dat_fish_s_dz_sp.order <- dat_fish_s_dz_sp |> 
  group_by(Species) |> 
  summarize(mean_density_100m2 = mean(mean_density_100m2), 
            mean_biomass_g_100m2 = mean(mean_biomass_g_100m2)) |> 
  mutate(sum_density_100m2 = sum(mean_density_100m2),
         sum_biomass_g_100m2 = sum(mean_biomass_g_100m2),
         prop_density_100m2 = mean_density_100m2/sum_density_100m2,
         prop_biomass_g_100m2 = mean_biomass_g_100m2/sum_biomass_g_100m2) |> 
  ungroup() |>
  arrange(desc(mean_biomass_g_100m2)) #can sort by numeric density or biomass density

dat_fish_s_dz_sp.order

# TO-DO: filter out rare species? Or just for plotting?
# w/ square-root (as compared to 4th-root) transformation, will make rare species less impactful
#At least 1% of biomass averaged across whole data set
dat_fish_s_dz_sp.order |> 
  filter(prop_biomass_g_100m2 > 0.01) |> 
  select(Species, prop_biomass_g_100m2) |>
  arrange(desc(prop_biomass_g_100m2))

#At least 1% of density averaged across whole data set
dat_fish_s_dz_sp.order |> 
  filter(prop_density_100m2 > 0.01) |> 
  select(Species, prop_density_100m2) |>
  arrange(desc(prop_density_100m2))

# define order for plotting/display of Species
dat_fish_s_dz_sp$Species <- factor(dat_fish_s_dz_sp$Species, 
                                   levels = dat_fish_s_dz_sp.order$Species)

# Species density and biomass by Region, AR_Complex and DepthZone
dat_fish_ARC_dz_sp <- dat_fish_s_dz_sp |> 
  group_by(Species, Region, AR_Complex, DepthZone) |> 
  summarize(mean_density_100m2 = mean(mean_density_100m2),
            mean_biomass_g_100m2 = mean(mean_biomass_g_100m2)) |> 
  ungroup()

# Convert mean Species density to proportions by AR_Complex and DepthZone
dat_fish_ARC_dz_sp <- dat_fish_ARC_dz_sp |> 
  group_by(Region, AR_Complex, DepthZone) |> 
  mutate(sum_density_100m2 = sum(mean_density_100m2),
         sum_biomass_g_100m2 = sum(mean_biomass_g_100m2),
         prop_density_100m2 = mean_density_100m2/sum_density_100m2,
         prop_biomass_g_100m2 = mean_biomass_g_100m2/sum_biomass_g_100m2) |> 
  ungroup() |> 
  mutate(Variables_Label = paste(AR_Complex, DepthZone, sep = " - ")) #create new column for plotting

#At least 1% of biomass in any AR_Complex + DepthZone
biomass_1perc_ARC_dz_sp <- dat_fish_ARC_dz_sp |> 
  filter(prop_biomass_g_100m2 > 0.01) |> 
  select(Species, prop_biomass_g_100m2) |>
  arrange(desc(prop_biomass_g_100m2)) |> 
  distinct(Species) |> 
  droplevels()

biomass_1perc_ARC_dz_sp

#SPECIES ORDER w/ "Other" at least 1% of biomass in any AR_Complex + DepthZone
biomass_1perc_ARC_dz_sp_ORDER <- dat_fish_ARC_dz_sp |> 
  filter(prop_biomass_g_100m2 > 0.01) |>
  group_by(Species) |>
  summarise(prop_biomass_g_100m2 = mean(prop_biomass_g_100m2)) |>
  arrange(desc(prop_biomass_g_100m2)) |>
  ungroup() |> 
  select(Species) |>
  droplevels() |> 
  add_row(Species = "Other")

biomass_1perc_ARC_dz_sp_ORDER

#At least 1% of density in any AR_Complex + DepthZone
density_1perc_ARC_dz_sp <- dat_fish_ARC_dz_sp |> 
  filter(prop_density_100m2 > 0.01) |> 
  select(Species, prop_density_100m2) |>
  arrange(desc(prop_density_100m2)) |> 
  distinct(Species) |> 
  print(n = Inf)|> 
  droplevels()

density_1perc_ARC_dz_sp

#SPECIES ORDER w/ "Other" at least 1% of density in any AR_Complex + DepthZone
density_1perc_ARC_dz_sp_ORDER <- dat_fish_ARC_dz_sp |> 
  filter(prop_density_100m2 > 0.01) |>
  group_by(Species) |>
  summarise(prop_density_100m2 = mean(prop_density_100m2)) |>
  arrange(desc(prop_density_100m2)) |>
  ungroup() |> 
  select(Species) |>
  droplevels() |> 
  add_row(Species = "Other")

density_1perc_ARC_dz_sp_ORDER


#Barplot of proportions of Species Biomass by AR_Complex and DepthZone
# TO DO - revise plot to be more focused on depthzone patterns?

barplot_biomass_ARC_dz_sp <- dat_fish_ARC_dz_sp |>
  #Any species <1% biomass in any AR_Complex + DepthZone is grouped as "Other"
  mutate(Species = case_when(Species %in% biomass_1perc_ARC_dz_sp$Species ~ Species, 
                             TRUE ~ "Other")) |>
  group_by(AR_Complex, DepthZone, Species, Variables_Label) |>
  summarize(prop_biomass_g_100m2 = sum(prop_biomass_g_100m2)) |>
  ungroup() |>
  droplevels() |> 
  mutate(Species = factor(Species, levels = biomass_1perc_ARC_dz_sp_ORDER$Species)) |>
  ggplot( aes(x = Variables_Label, 
              y = prop_biomass_g_100m2, 
              fill = Species)) +
  #geom_bar_pattern replaces geom_bar - part of the ggpattern package needed to add diagonal lines to algae categories
  geom_bar(stat="identity", 
           position = "fill") +
  #https://github.com/EmilHvitfeldt/r-color-palettes/blob/master/type-sorted-palettes.md
  scale_fill_paletteer_d("Polychrome::palette36") +
  facet_grid(AR_Complex ~., scales = "free", space = "free") +
  coord_flip() +
  labs(x = " ", 
       y = "Proportion Biomass") +
  theme_bw() +
  theme(axis.text = element_text(size = 13),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 15)) +
  theme(legend.title = element_blank(),
        legend.position = "top") +
  theme(strip.text = element_text(size = 11),
        strip.background = element_rect(color = "black", fill = "white"))

barplot_biomass_ARC_dz_sp

ggsave("barplot_biomass_ARC_dz_sp.png", barplot_biomass_ARC_dz_sp,
       width = 10,
       height = 9, 
       dpi = 600)

# Prep data for multivariate analysis (community table format) ----------------

# TO DO - (both Biomass + Density) average to AR_Complex (Regions) + DepthZone 
# and use those as points for overall nMDS & stacked bar plots
# 
# Then for AR_Complex + most similar regions/depthzones repeat analyses at the 
# site/module level.

## Site + DepthZone level data
# square root transform numerical density and biomass
dat_fish_s_dz_sp <- dat_fish_s_dz_sp |> 
  mutate(sqrt_mean_density_100m2 = mean_density_100m2^0.5, #sqr-root transformed  
         sqrt_mean_biomass_g_100m2 = mean_biomass_g_100m2^0.5) #sqr-root transformed  

#create community data using pivot_wider (Region, Reef Type)
wide_dat_fish_s_dz_sp <- dat_fish_s_dz_sp |>
  mutate(plot_label = paste(AR_Complex, Site, DepthZone)) |>
  pivot_wider(id_cols = c(plot_label, Region, AR_Complex, Site, DepthZone),
              names_from = Species, 
              values_from = sqrt_mean_biomass_g_100m2,
              values_fill = 0)

## make "community data" table format to go into vegan package functions
comm_wide_dat_fish_s_dz_sp <- wide_dat_fish_s_dz_sp |> 
  column_to_rownames(var = "plot_label") |>
  select(-Region, -AR_Complex, -Site, -DepthZone) #remove the columns to get into comm format 



# nMDS Site DepthZone Analysis ---------------------------------------------------

#create nMDS
s_dz_NMDS <- metaMDS(comm_wide_dat_fish_s_dz_sp, 
                     distance = "bray", #using bray-curtis dissimilarity 
                     autotransform = FALSE, trymax=2000)

#extract nMDS stress for plotting
s_dz_NMDS_stress <-round(s_dz_NMDS$stress, 2)
# <0.2 good, <0.3 ok according to some soures, >0.3 poor (and/but close to 0 can also be an issue)
s_dz_NMDS_stress

# put output into tibble format & add back columns for plotting
tibble_s_dz_NMDS <- fortify(s_dz_NMDS) |>
  filter(score == "sites") |>
  rename(plot_label = "label") |>
  left_join(wide_dat_fish_s_dz_sp)

# extract species scores
tibble_s_dz_NMDS_species <- fortify(s_dz_NMDS) |>
  filter(score == "species")


# Plot nMDS ----------------------------------------------------------------

#To Do - separate out PVR - DepthZone ARM and/or Region Artificial reef
# Show both depth zone and region (filled shapes and/or labels reduced region?)

## ARs filled circles & NRs filled shapes?
## certain colors by AR_Complex & DepthZone
## AR complex colors differ between NR Regions and AR complexes (by latitude warm/cold?)
## How to distinguish between AR_Complex and DepthZone?

## ADD MCPs or ellipses around AR_Complex + DepthZone?


plot_s_dz_NMDS <- ggplot(tibble_s_dz_NMDS, 
                         aes(x = NMDS1, y = NMDS2,
                             color = DepthZone
                         )) +
  geom_point() +
  # stat_ellipse() +
  scale_shape_manual(values=c(16, 18, 15, 3, 0)) +
  geom_text_repel(aes(label = Site), 
                  hjust = 0.3, 
                  nudge_x = 0.03, 
                  min.segment.length = 1,
                  size = 2) + #text size of location labels
  theme_classic() +
  theme(panel.background = element_rect(fill = NA, colour = "black", linewidth = 1, linetype = "solid")) +
  #nMDS should not have axes
  xlab("") +
  ylab("") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_text(aes(x = Inf, y = Inf, label = paste("2D Stress:", s_dz_NMDS_stress)),
            vjust = "inward", 
            hjust = "inward",
            color = "black")


print(plot_s_dz_NMDS)

#"Species" Scores
comm_wide_dat_fish_s_dz_sp.spp.fit <- envfit(s_dz_NMDS, 
                                             comm_wide_dat_fish_s_dz_sp, 
                                             permutations = 999) # this fits species vectors

Spp_comm_wide_dat_fish_s_dz_sp.scrs <- as_tibble(scores(comm_wide_dat_fish_s_dz_sp.spp.fit, 
                                                        display = "vectors"), 
                                                 rownames = "Species") |> #save species intrinsic values into dataframe
  mutate(pval = comm_wide_dat_fish_s_dz_sp.spp.fit$vectors$pvals) |> #add pvalues to dataframe so you can select species which are significant
  filter(pval <= 0.05) |>
  # can make species arrow vectors shorter by reducing the multiplier to <1
  mutate(NMDS1 = NMDS1 * 1,
         NMDS2 = NMDS2 * 1)

#Plot "species" vectors on top of nMDS
plot_s_dz_NMDS_final <- plot_s_dz_NMDS +
  geom_segment(data = Spp_comm_wide_dat_fish_s_dz_sp.scrs, 
               aes(x = 0, 
                   xend = NMDS1, 
                   y = 0, 
                   yend = (NMDS2*-1), 
                   shape = NULL, 
                   fill =  NULL), # need fill Null because specified in mds above
               arrow = arrow(length = unit(0.2, "cm")), #size of the arrow head
               colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  geom_text_repel(data = Spp_comm_wide_dat_fish_s_dz_sp.scrs, 
                  aes(x=NMDS1, 
                      y=(NMDS2*-1), 
                      label = Species, 
                      shape = NULL, 
                      color = NULL, 
                      fill =  NULL), # need color Null because specified in mds above
                  cex = 4, #text size of the arrow labels
                  direction = "y", 
                  segment.size = 0.25) #add labels for species, use ggrepel

print(plot_s_dz_NMDS_final)

#save plot as png file
ggsave("plot_s_dz_NMDS_final.png", plot_s_dz_NMDS_final,
       width = 9,
       height = 7, 
       dpi = 600)

# Cluster analysis-------------------------------------------------------------
#Create a distance matrix based on the community assemblages
dis.comm_wide_dat_fish_s_dz_sp <- vegdist(comm_wide_dat_fish_s_dz_sp, method = "bray")

#Create a cluster dendrogram
clust.comm_wide_dat_fish_s_dz_sp <- hclust(dis.comm_wide_dat_fish_s_dz_sp, "average")

## create dendogram object
dend.comm_wide_dat_fish_s_dz_sp <- as.dendrogram(clust.comm_wide_dat_fish_s_dz_sp)

# Plot the dendrogram
library(dendextend)
plot(dend.comm_wide_dat_fish_s_dz_sp, # 
     horiz = T,
     xlab = "Dissimilarity")

# from Brenda's code (I think the c() order allows you to rotate various branches to make it read better)
# par(mar=c(5,1,1,12))
# plot(rotate(dend.comm_wide_dat_fish_s_dz_sp, c(1:15, 20, 24, 26:27, 25, 29:30, 28, 23:21, 16:17, 18:19)), # 
#      horiz = T,
#      xlab = "Dissimilarity",
#      xlim = c(0.4, 0))

# Color the branches based on height
dend_color.comm_wide_dat_fish_s_dz_sp <- dend.comm_wide_dat_fish_s_dz_sp %>% color_branches(k = 20)  # You can adjust the number of colors (k) as needed

plot(dend_color.comm_wide_dat_fish_s_dz_sp, # 
     horiz = T,
     xlab = "Dissimilarity")

# pairwise adonis (PERMANOVA) -------------------------------------------------

#######################
## Pairwise Adonis tests - ADAPTED FROM JTC's MPA BASELINE CODE (and used in platform project)


## test if there is significant difference in dispersion among groups  # if there is, this doesn't disqualify results, but there is potential that it could affect them. 
betadis.ARC_dz<-betadisper(dis.comm_wide_dat_fish_s_dz_sp,
                           paste(wide_dat_fish_s_dz_sp$AR_Complex , wide_dat_fish_s_dz_sp$DepthZone))
anova(betadis.ARC_dz)
permutest(betadis.ARC_dz)
TukeyHSD(betadis.ARC_dz)
THSD.betadis.ARC_dz<-TukeyHSD(betadis.ARC_dz)
# output just significant pairwise beta disper comparisons
bd_sig_pairs<-as.matrix(round(THSD.betadis.ARC_dz$group [which(THSD.betadis.ARC_dz$group[,4] < 0.05),4],3))
data.frame(Regions=rownames(bd_sig_pairs),p_adj = bd_sig_pairs[,1])
write.table(data.frame(ARC_dz=rownames(bd_sig_pairs),p_adj = bd_sig_pairs[,1]), file="betadis.ARC_dz.pairwise.sig_table.csv", sep=",", quote = FALSE, row.names=F)

levels(dat_fish_t$AR_Complex)

#Pairwise comparisons Using adonis (instead of ANOSIM)
ARC_DZ_list<-as.character(sort(unique(paste(wide_dat_fish_s_dz_sp$AR_Complex , wide_dat_fish_s_dz_sp$DepthZone))))  # create vector of tt_dz_list
reg_pairs<-combn(ARC_DZ_list, 2)    #create list of tt_dz pairs
RESULTS<-data.frame(reg_1=factor(reg_pairs[1,]),
                    reg_2=factor(reg_pairs[2,]),
                    R2=NA,
                    PrF=NA )  # create results table

for(i in 1:(length(reg_pairs)/2)){
  #i=11 # for testing
  #reg_pairs[,i]
  # select ARC_DZ pair vector
  ARC_DZ <- paste(wide_dat_fish_s_dz_sp$AR_Complex , wide_dat_fish_s_dz_sp$DepthZone)
  regs_pair <- ARC_DZ[ARC_DZ %in%  reg_pairs[,i]]
  # select data just for sites in the ARC_DZ pair
  comm_pair<-comm_wide_dat_fish_s_dz_sp[ARC_DZ %in%  reg_pairs[,i],]        
  ado_pair<-adonis2(comm_pair ~ regs_pair, comm_pair)
  RESULTS[i,3]<- ado_pair$R2[1]
  RESULTS[i,4]<- ado_pair[1,"Pr(>F)"]
}

RESULTS
RESULTS$R2<-round(RESULTS$R2,2)
#dcast(RESULTS, reg_2 ~ reg_1,  value.var="R2")
RESULTS$p.symb[RESULTS$PrF>=0.05]  <- ""
RESULTS$p.symb[RESULTS$PrF<0.005]  <- "**"
RESULTS$p.symb[RESULTS$PrF<0.05 & RESULTS$PrF>=0.005]  <- "*"
RESULTS$R_p.symb<-paste(RESULTS$R2,RESULTS$p.symb, sep="")
RESULTS <- RESULTS %>% arrange(reg_1, reg_2)
# select those with 3 or more sites (the rest will be removed from results table)
site_count <- wide_dat_fish_s_dz_sp |>  
  group_by(AR_Complex, DepthZone) |>  
  summarize(site_count= n()) |> 
  filter(site_count>3)
min_site_count <- paste(site_count$AR_Complex,
                        site_count$DepthZone)
# replace selected results from table
RESULTS$R_p.symb[!(RESULTS$reg_1 %in% min_site_count)] <- "-" # select those not in list of 3 or more sites
RESULTS$R_p.symb[!(RESULTS$reg_2 %in% min_site_count)] <- "-" # select those not in list of 3 or more sites
display_RESULTS<-dcast(RESULTS, reg_2 ~ reg_1,  value.var="R_p.symb")
display_RESULTS
#write RESULTS to a comma separated text file
write.table(display_RESULTS, file="ADONIS_pairwise_ARC_DZ.csv", sep=",", quote = FALSE, row.names=F)


# using function Jonathan found somewhere
source("~/Dropbox/VRG Files/R Code/custom_functions/parwise.adonis.r")
pairwise_dz_ARC <- pairwise.adonis(
  wide_dat_fish_s_dz_sp[, 6:ncol(wide_dat_fish_s_dz_sp)],
  paste(wide_dat_fish_s_dz_sp$AR_Complex , wide_dat_fish_s_dz_sp$DepthZone),
  sim.function = 'vegdist',
  sim.method = 'bray',
  p.adjust.m = 'none' #permutation tests shouldn't need multiple comparison test p-value corrections
)
pairwise_dz_ARC <- pairwise_dz_ARC |> 
  arrange(desc(R2)) #sort by R2 value



# nMDS AR_Complex DepthZone Analysis ---------------------------------------

## AR_Complex + DepthZone level communuty data table

# square root transform numerical density and biomass
dat_fish_ARC_dz_sp <- dat_fish_ARC_dz_sp |> 
  mutate(sqrt_mean_density_100m2 = mean_density_100m2^0.5, #sqr-root transformed  
         sqrt_mean_biomass_g_100m2 = mean_biomass_g_100m2^0.5) #sqr-root transformed  

#create community data using pivot_wider (Region, Reef Type)
wide_dat_fish_ARC_dz_sp <- dat_fish_ARC_dz_sp |>
  mutate(plot_label = paste(AR_Complex, DepthZone)) |>
  pivot_wider(id_cols = c(plot_label, Region, AR_Complex, DepthZone),
              names_from = Species,
              #PICK EITHER DENSITY OR BIOMASS (sqrt transformed)
              values_from = sqrt_mean_biomass_g_100m2,
              values_fill = 0)

## make "community data" table format to go into vegan package functions
comm_wide_dat_fish_ARC_dz_sp <- wide_dat_fish_ARC_dz_sp |> 
  column_to_rownames(var = "plot_label") |>
  select(-Region, -AR_Complex, -DepthZone) #remove the columns to get into comm format 


#create nMDS
ARC_dz_NMDS <- metaMDS(comm_wide_dat_fish_ARC_dz_sp, 
                       distance = "bray", #using bray-curtis dissimilarity 
                       autotransform = FALSE, trymax=2000)

#extract nMDS stress for plotting
ARC_dz_NMDS_stress <-round(ARC_dz_NMDS$stress, 2)
# <0.2 good, <0.3 ok according to some soures, >0.3 poor (and/but close to 0 can also be an issue)
ARC_dz_NMDS_stress

# put output into tibble format & add back columns for plotting
tibble_ARC_dz_NMDS <- fortify(ARC_dz_NMDS) |>
  filter(score == "sites") |>
  rename(plot_label = "label") |> 
  left_join(wide_dat_fish_ARC_dz_sp)

# extract species scores
tibble_ARC_dz_NMDS_species <- fortify(ARC_dz_NMDS) |>
  filter(score == "species")


# Plot nMDS ----------------------------------------------------------------

#To Do - separate out PVR - DepthZone ARM and/or Region Artificial reef
# Show both depth zone and region (filled shapes and/or labels reduced region?)

## ARs filled circles & NRs filled shapes?
## certain colors by AR_Complex & DepthZone
## AR complex colors differ between NR Regions and AR complexes (by latitude warm/cold?)
## How to distinguish between AR_Complex and DepthZone?

## MCPs or ellipses around AR_Complex + DepthZone?


plot_ARC_dz_NMDS <- ggplot(tibble_ARC_dz_NMDS, 
                           aes(x = NMDS1, y = NMDS2,
                               color = AR_Complex,
                               shape = DepthZone)) + #sizing points by temp
  geom_point() +
  # stat_ellipse() +
  scale_shape_manual(values=c(16, 18, 15, 3, 0)) +
  geom_text_repel(aes(label = plot_label), 
                  hjust = 0.3, 
                  nudge_x = 0.03, 
                  min.segment.length = 1,
                  size = 2) + #text size of location labels
  theme_classic() +
  theme(panel.background = element_rect(fill = NA, colour = "black", linewidth = 1, linetype = "solid")) +
  #nMDS should not have axes
  xlab("") +
  ylab("") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  geom_text(aes(x = Inf, y = Inf, label = paste("2D Stress:", ARC_dz_NMDS_stress)),
            vjust = "inward", 
            hjust = "inward",
            color = "black")


print(plot_ARC_dz_NMDS)

#"Species" Scores
comm_wide_dat_fish_ARC_dz_sp.spp.fit <- envfit(ARC_dz_NMDS, 
                                               comm_wide_dat_fish_ARC_dz_sp, 
                                               permutations = 999) # this fits species vectors

Spp_comm_wide_dat_fish_ARC_dz_sp.scrs <- as_tibble(scores(comm_wide_dat_fish_ARC_dz_sp.spp.fit, 
                                                          display = "vectors"), 
                                                   rownames = "Species") |> #save species intrinsic values into dataframe
  mutate(pval = comm_wide_dat_fish_ARC_dz_sp.spp.fit$vectors$pvals) |> #add pvalues to dataframe so you can select species which are significant
  filter(pval <= 0.05) |>
  # can make species arrow vectors shorter by reducing the multiplier to <1
  mutate(NMDS1 = NMDS1 * 1,
         NMDS2 = NMDS2 * 1)

#Plot "species" vectors on top of nMDS
plot_ARC_dz_NMDS_final <- plot_ARC_dz_NMDS +
  geom_segment(data = Spp_comm_wide_dat_fish_ARC_dz_sp.scrs, 
               aes(x = 0, 
                   xend = NMDS1, 
                   y = 0, 
                   yend = (NMDS2*-1), 
                   shape = NULL, 
                   fill =  NULL), # need fill Null because specified in mds above
               arrow = arrow(length = unit(0.2, "cm")), #size of the arrow head
               colour = "grey10", lwd=0.3) + #add vector arrows of significant species
  geom_text_repel(data = Spp_comm_wide_dat_fish_ARC_dz_sp.scrs, 
                  aes(x=NMDS1, 
                      y=(NMDS2*-1), 
                      label = Species, 
                      shape = NULL, 
                      color = NULL, 
                      fill =  NULL), # need color Null because specified in mds above
                  cex = 4, #text size of the arrow labels
                  direction = "y", 
                  segment.size = 0.25) #add labels for species, use ggrepel

print(plot_ARC_dz_NMDS_final)

#save plot as png file
ggsave("plot_ARC_dz_NMDS_final.png", plot_s_dz_NMDS_final,
       width = 9,
       height = 7, 
       dpi = 600)


# Cluster analysis-------------------------------------------------------------
#Create a distance matrix based on the community assemblages
dis.comm_wide_dat_fish_ARC_dz_sp <- vegdist(comm_wide_dat_fish_ARC_dz_sp, method = "bray")

#Create a cluster dendrogram
clust.comm_wide_dat_fish_ARC_dz_sp <- hclust(dis.comm_wide_dat_fish_ARC_dz_sp, "average")

## create dendogram object
dend.comm_wide_dat_fish_ARC_dz_sp <- as.dendrogram(clust.comm_wide_dat_fish_ARC_dz_sp)

# Plot the dendrogram
library(dendextend)
plot(dend.comm_wide_dat_fish_ARC_dz_sp, # 
     horiz = T,
     xlab = "Dissimilarity")

# from Brenda's code (I think the c() order allows you to rotate various branches to make it read better)
# par(mar=c(5,1,1,12))
# plot(rotate(dend.comm_wide_dat_fish_ARC_dz_sp, c(1:15, 20, 24, 26:27, 25, 29:30, 28, 23:21, 16:17, 18:19)), # 
#      horiz = T,
#      xlab = "Dissimilarity",
#      xlim = c(0.4, 0))

# Color the branches based on height
dend_color.comm_wide_dat_fish_ARC_dz_sp <- dend.comm_wide_dat_fish_ARC_dz_sp %>% color_branches(k = 7)  # You can adjust the number of colors (k) as needed

plot(dend_color.comm_wide_dat_fish_ARC_dz_sp, # 
     horiz = T,
     xlab = "Dissimilarity")