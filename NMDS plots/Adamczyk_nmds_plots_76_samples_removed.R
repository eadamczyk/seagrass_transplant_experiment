# Project: Seagrass (Zostera marina) transplant experiment reveals core microbiome and resistance to environmental change 
# Code by: Emily Adamczyk
# Purpose: NMDS plots - 76 samples removed


# Libraries
library(vegan)
library(ggplot2)
library(reshape2)
library(dplyr)

# setting working directory
setwd("C:/Users/Emily/Documents/seagrass_transplant_experiment/data_analysis")

###Read table metadata and abundances
metadata_NMDS <- read.csv("seagrass_experiment_16s_meta_otu_900_RAREFIED_76removed.csv", header=T)

# recode sand spit and wolf to reflect choked region, renaming treatments
metadata_NMDS <- metadata_NMDS %>% 
  dplyr::mutate(treatment=recode(treatment, "control" = "live_eelgrass"),
                treatment=recode(treatment, "sterile" = "sterile_eelgrass"),
                origin_site=recode(origin_site, "sandspit" = "choked_sand_spit"),
                origin_site=recode(origin_site, "wolf" = "choked_wolf"),
                destination=recode(destination, "sandspit" = "choked_sand_spit"),
                destination=recode(destination, "wolf" = "choked_wolf"))


#######################################################################################
#######################################################################################
#######################################################################################

# *** Note: Eelgrass "control" and "live" treatments are synonymous with "natural" (all refer to natural eelgrass)
# *** Note: Eelgrass "asu" (artificial seagrass unit) treatment is the same as artificial eelgrass


#### creating transplant column to make plots ####
names(metadata_NMDS)
metadata_NMDS$transplant_type = NA

######################
### NATURAL SHOOTS ###
######################

# initial
metadata_NMDS$transplant_type[metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "initial_eelgrass_swab"] <- "initial"

# replant
metadata_NMDS$transplant_type[as.character(metadata_NMDS$origin_site) == as.character(metadata_NMDS$destination) & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "replant"

# transplant
metadata_NMDS$transplant_type[as.character(metadata_NMDS$origin_site) != as.character(metadata_NMDS$destination) & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant"

######################
### STERILE SHOOTS ###
######################

# initial
metadata_NMDS$transplant_type[metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "initial_eelgrass_swab"] <- "initial"

# replant
metadata_NMDS$transplant_type[as.character(metadata_NMDS$origin_site) == as.character(metadata_NMDS$destination) & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "replant"

# transplant
metadata_NMDS$transplant_type[as.character(metadata_NMDS$origin_site) != as.character(metadata_NMDS$destination) & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant"

#########################
### ARTIFICIAL SHOOTS ###
#########################

# initial
metadata_NMDS$transplant_type[metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "initial_eelgrass_swab"] <- "initial"

# replant
metadata_NMDS$transplant_type[as.character(metadata_NMDS$origin_site) == as.character(metadata_NMDS$destination) & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "replant"

# transplant
metadata_NMDS$transplant_type[as.character(metadata_NMDS$origin_site) != as.character(metadata_NMDS$destination) & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant"

##################################################

#### creating secondary transplant column ####
names(metadata_NMDS)
metadata_NMDS$transplant_type_2 = NA

####################
# NATURAL SEAGRASS #
####################

# initial
metadata_NMDS$transplant_type_2[metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "initial_eelgrass_swab"] <- "initial"

# replant
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == as.character(metadata_NMDS$destination) & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "replant"

### transplant same region ###
# transplant same region (pruth bay to pruth pocket)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_bay" & as.character(metadata_NMDS$destination) == "pruth_pocket" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

# transplant same region (pruth pocket to pruth bay)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_pocket" & as.character(metadata_NMDS$destination) == "pruth_bay" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

# transplant same region (choked wolf to choked sand spit)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_wolf" & as.character(metadata_NMDS$destination) == "choked_sand_spit" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

# transplant same region (choked sand spit to choked wolf)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_sand_spit" & as.character(metadata_NMDS$destination) == "choked_wolf" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

### transplant different region ###
# transplant different region (pruth bay to choked sand spit)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_bay" & as.character(metadata_NMDS$destination) == "choked_sand_spit" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (pruth bay to choked wolf)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_bay" & as.character(metadata_NMDS$destination) == "choked_wolf" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (pruth pocket to choked sand spit)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_pocket" & as.character(metadata_NMDS$destination) == "choked_sand_spit" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (pruth pocket to choked wolf)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_pocket" & as.character(metadata_NMDS$destination) == "choked_wolf" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked wolf to pruth bay)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_wolf" & as.character(metadata_NMDS$destination) == "pruth_bay" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked wolf to pruth pocket)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_wolf" & as.character(metadata_NMDS$destination) == "pruth_pocket" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked sand spit to pruth bay)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_sand_spit" & as.character(metadata_NMDS$destination) == "pruth_bay" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked sand spit to pruth pocket)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_sand_spit" & as.character(metadata_NMDS$destination) == "pruth_pocket" & metadata_NMDS$treatment == "live_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"


####################
# STERILE SEAGRASS #
####################

# initial
metadata_NMDS$transplant_type_2[metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "initial_eelgrass_swab"] <- "initial"

# replant
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == as.character(metadata_NMDS$destination) & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "replant"

### transplant same region ###
# transplant same region (pruth bay to pruth pocket)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_bay" & as.character(metadata_NMDS$destination) == "pruth_pocket" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

# transplant same region (pruth pocket to pruth bay)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_pocket" & as.character(metadata_NMDS$destination) == "pruth_bay" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

# transplant same region (choked wolf to choked sand spit)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_wolf" & as.character(metadata_NMDS$destination) == "choked_sand_spit" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

# transplant same region (choked sand spit to choked wolf)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_sand_spit" & as.character(metadata_NMDS$destination) == "choked_wolf" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

### transplant different region ###
# transplant different region (pruth bay to choked sand spit)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_bay" & as.character(metadata_NMDS$destination) == "choked_sand_spit" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (pruth bay to choked wolf)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_bay" & as.character(metadata_NMDS$destination) == "choked_wolf" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (pruth pocket to choked sand spit)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_pocket" & as.character(metadata_NMDS$destination) == "choked_sand_spit" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (pruth pocket to choked wolf)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_pocket" & as.character(metadata_NMDS$destination) == "choked_wolf" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked wolf to pruth bay)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_wolf" & as.character(metadata_NMDS$destination) == "pruth_bay" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked wolf to pruth pocket)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_wolf" & as.character(metadata_NMDS$destination) == "pruth_pocket" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked sand spit to pruth bay)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_sand_spit" & as.character(metadata_NMDS$destination) == "pruth_bay" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked sand spit to pruth pocket)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_sand_spit" & as.character(metadata_NMDS$destination) == "pruth_pocket" & metadata_NMDS$treatment == "sterile_eelgrass" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

#######################
# ARTIFICIAL SEAGRASS #
#######################

# initial
metadata_NMDS$transplant_type_2[metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "initial_eelgrass_swab"] <- "initial"

# replant
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == as.character(metadata_NMDS$destination) & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "replant"

### transplant same region ###
# transplant same region (pruth bay to pruth pocket)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_bay" & as.character(metadata_NMDS$destination) == "pruth_pocket" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

# transplant same region (pruth pocket to pruth bay)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_pocket" & as.character(metadata_NMDS$destination) == "pruth_bay" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

# transplant same region (choked wolf to choked sand spit)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_wolf" & as.character(metadata_NMDS$destination) == "choked_sand_spit" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

# transplant same region (choked sand spit to choked wolf)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_sand_spit" & as.character(metadata_NMDS$destination) == "choked_wolf" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_same_region"

### transplant different region ###
# transplant different region (pruth bay to choked sand spit)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_bay" & as.character(metadata_NMDS$destination) == "choked_sand_spit" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (pruth bay to choked wolf)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_bay" & as.character(metadata_NMDS$destination) == "choked_wolf" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (pruth pocket to choked sand spit)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_pocket" & as.character(metadata_NMDS$destination) == "choked_sand_spit" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (pruth pocket to choked wolf)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "pruth_pocket" & as.character(metadata_NMDS$destination) == "choked_wolf" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked wolf to pruth bay)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_wolf" & as.character(metadata_NMDS$destination) == "pruth_bay" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked wolf to pruth pocket)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_wolf" & as.character(metadata_NMDS$destination) == "pruth_pocket" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked sand spit to pruth bay)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_sand_spit" & as.character(metadata_NMDS$destination) == "pruth_bay" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"

# transplant different region (choked sand spit to pruth pocket)
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$origin_site) == "choked_sand_spit" & as.character(metadata_NMDS$destination) == "pruth_pocket" & metadata_NMDS$treatment == "asu" & metadata_NMDS$category == "final_eelgrass_swab"] <- "transplant_different_region"



#######################################################################################
#######################################################################################
#######################################################################################

############################################
####     Environment Settings         ######
# Making custom ggplot theme for all plots #
############################################

theme.nmds <- function(){
  theme_bw()+
    theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
           axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
           axis.text = element_text(size = 16), #font size of numbers in axis
           panel.grid.major = element_blank(), #remove major grid
           panel.grid.minor = element_blank(), #remove minor grid
           axis.line = element_line(colour = "black"), #draw line in the axis
           panel.border = element_blank(), #remove lines outside the graph
           legend.title=element_text(size = 16), #remove legend title
           legend.direction = "vertical", #direction
           legend.justification = c(1, 1), legend.position = "right", #legend is top right
           legend.key.size = unit(2.0, 'lines'), #spacing between legends
           legend.text = element_text(size = 16), #font size of legend
           plot.title = element_text(hjust = 0.5, size = 16)) #center plot title and set font size
}

####################################################
# renaming treatments, sites, and transplant types #
####################################################

metadata_NMDS <- metadata_NMDS %>% 
  dplyr::mutate(treatment=recode(treatment, "live_eelgrass" = "Live seagrass"),
                treatment=recode(treatment, "sterile_eelgrass" = "Sterilized seagrass"),
                treatment=recode(treatment, "asu" = "Artificial seagrass"),
                destination=recode(destination, "choked_sand_spit" = "Choked Sand Spit"),
                destination=recode(destination, "choked_wolf" = "Choked Wolf"),
                destination=recode(destination, "pruth_bay" = "Pruth Bay"),
                destination=recode(destination, "pruth_pocket" = "Pruth Pocket"),
                origin_site=recode(origin_site, "choked_sand_spit" = "Choked Sand Spit"),
                origin_site=recode(origin_site, "choked_wolf" = "Choked Wolf"),
                origin_site=recode(origin_site, "pruth_bay" = "Pruth Bay"),
                origin_site=recode(origin_site, "pruth_pocket" = "Pruth Pocket"),
                transplant_type_2=recode(transplant_type_2, "initial" = "Initial (Day 1)"),
                transplant_type_2=recode(transplant_type_2, "replant" = "Replant (Day 5)"),
                transplant_type_2=recode(transplant_type_2, "transplant_different_region" = "Transplant different region (Day 5)"),
                transplant_type_2=recode(transplant_type_2, "transplant_same_region" = "Transplant same region (Day 5)"))

#########################
# RENAMING WATER LABELS #
#########################

# Pruth Bay
metadata_NMDS$destination[as.character(metadata_NMDS$origin_site) == "Pruth Bay" & as.character(metadata_NMDS$destination) == "na"] <- "Pruth Bay"
# Pruth Pocket
metadata_NMDS$destination[as.character(metadata_NMDS$origin_site) == "Pruth Pocket" & as.character(metadata_NMDS$destination) == "na"] <- "Pruth Pocket"
# Choked Sand Spit
metadata_NMDS$destination[as.character(metadata_NMDS$origin_site) == "Choked Sand Spit" & as.character(metadata_NMDS$destination) == "na"] <- "Choked Sand Spit"
# Choked Wolf
metadata_NMDS$destination[as.character(metadata_NMDS$origin_site) == "Choked Wolf" & as.character(metadata_NMDS$destination) == "na"] <- "Choked Wolf"
# water final
metadata_NMDS$category[as.character(metadata_NMDS$treatment) == "recovery" & as.character(metadata_NMDS$category) == "water_filters_(0.22um)"] <- "water_final"
# water initial
metadata_NMDS$category[as.character(metadata_NMDS$treatment) == "deployment" & as.character(metadata_NMDS$category) == "water_filters_(0.22um)"] <- "water_initial"
# water name
metadata_NMDS <- metadata_NMDS %>% 
  dplyr::mutate(treatment=recode(treatment, "recovery" = "Water"),
                treatment=recode(treatment, "deployment" = "Water"))
# transplant_type_2 labels
# water final
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$category) == "water_final"] <- "water_final"
# water final
metadata_NMDS$transplant_type_2[as.character(metadata_NMDS$category) == "water_initial"] <- "water_initial"



#######################################################################################
#######################################################################################
#######################################################################################

########################################
### NATURAL SEGRASS DAY 1 - FIGURE 3 ###
########################################

###select initial swab samples
initial_eelgrass_swab <- c("initial_eelgrass_swab")
initial_only <- metadata_NMDS %>% 
  dplyr::filter(category %in% initial_eelgrass_swab )

#View(initial_only)

###select control (live) only
control <- c("Live seagrass")
initial_control <- initial_only %>% 
  dplyr::filter(treatment %in% control)

#View(initial_control)

### Creating an object to store abundances only (remove first 6 columns with dplyr)
abundances_NMDS <- initial_control %>% 
  dplyr::select(-c(1:13, 1688:1689))

#Get MDS stats using matrix with only numericals (without sites' column) 
set.seed(2)
NMDS.initial.control <- metaMDS(log(abundances_NMDS+1), distance = "bray", k=2)  
NMDS.initial.control

stressplot(NMDS.initial.control)
plot(NMDS.initial.control)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.initial.control$points[,1]
MDS2 = NMDS.initial.control$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, origin_site = initial_control$origin_site, live_seagrass = initial_control$treatment, sample_id = initial_control$sample_id)
NMDS

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS, old=c("MDS1", "MDS2", "origin_site", "live_seagrass"), new=c("NMDS1","NMDS2", "origin_site", "live_seagrass"))
NMDS

# re-order the factor levels before the plot
NMDS$origin_site <- factor(NMDS$origin_site, levels=c("Pruth Bay", "Pruth Pocket", "Choked Wolf", "Choked Sand Spit"))
NMDS$live_seagrass <- factor(NMDS$live_seagrass, levels=c("Live seagrass"))

# NMDS plot - Fig. 3
p2 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, shape = live_seagrass, colour= origin_site)) +
  geom_point(size = 6, alpha = 0.9) +
  theme.nmds() +
  scale_shape_manual(values = c(16, 15, 17, 9)) +
  scale_colour_manual(values=c("dodgerblue", "blue4", "tomato", "darkred"))

p2


#######################################################################################
#######################################################################################
#######################################################################################

##################################################################################################
### FINAL NATURAL, ARTIFICIAL, AND STERILE DAY 5 (transplant and replant same region) - FIGURE 4 #
##################################################################################################

###select final swab samples
final <- c("final_eelgrass_swab")
final_only <- metadata_NMDS %>% 
  dplyr::filter(category %in% final)

###select final natural, artificial, sterile
all <- c("Live seagrass", "Artificial seagrass", "Sterilized seagrass")
all_final <- final_only %>% 
  dplyr::filter(treatment %in% all)

###select Transplant same region (Day 5), Replant (Day 5), and final water
region_final <- c("Transplant same region (Day 5)", "Replant (Day 5)")
all_region_final <- all_final %>% 
  dplyr::filter(transplant_type_2 %in% region_final)
arf <- all_region_final 


### Creating an object to store abundances only (remove first 6 columns with dplyr)
abundances_NMDS <- arf %>% 
  dplyr::select(-c(1:13, 1688:1689))

#Get MDS stats using matrix with only numericals (without sites' column) 
set.seed(2)
NMDS.arf <- metaMDS(log(abundances_NMDS+1), distance = "bray", k=2)  
NMDS.arf

stressplot(NMDS.arf)
plot(NMDS.arf)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.arf$points[,1]
MDS2 = NMDS.arf$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Destination_site = arf$destination, Treatment = arf$treatment)
NMDS

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS, old=c("MDS1", "MDS2", "Destination_site", "Treatment"), new=c("NMDS1","NMDS2", "Destination_site", "Treatment"))
NMDS

# re-order the factor levels before the plot
NMDS$Destination_site <- factor(NMDS$Destination_site, levels=c("Pruth Bay", "Pruth Pocket", "Choked Wolf", "Choked Sand Spit"))
NMDS$Destination_site
NMDS$Treatment <- factor(NMDS$Treatment, levels=c("Live seagrass", "Sterilized seagrass", "Artificial seagrass"))
NMDS$Treatment

# NMDS plot - Figure 4
p2 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, shape = Treatment, colour=Destination_site)) +
  geom_point(size = 5, alpha = 0.8) +
  theme.nmds() +
  scale_shape_manual(values = c(16, 15, 17)) +
  #labs(title= "Final Water, ASU, live, sterile") +
  scale_colour_manual(values=c("dodgerblue", "blue4", "tomato", "darkred")) #+
#stat_ellipse(linetype = 1)

p2

#######################################################################################
#######################################################################################
#######################################################################################


#######################################################
### NATURAL ORIGIN AND DESTINATION DAY 5 - FIGURE 6A ##
#######################################################

# final NATURAL only
final <- ("final_eelgrass_swab")
live <- ("Live seagrass")
u_f <- metadata_NMDS %>% 
  dplyr::filter(category %in% final) %>% 
  dplyr::filter(treatment %in% live)

### Creating an object to store abundances only (remove first 6 columns with dplyr)
abundances_NMDS <- u_f %>% 
  dplyr::select(-c(1:13, 1688:1689))

#Get MDS stats using matrix with only numericals (without sites' column) 
set.seed(2)
NMDS.u_f <- metaMDS(log(abundances_NMDS+1), distance = "bray", k=2)  
NMDS.u_f

stressplot(NMDS.u_f)
plot(NMDS.u_f)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.u_f$points[,1]
MDS2 = NMDS.u_f$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, origin_site = u_f$origin_site, destination = u_f$destination, sample_id = u_f$sample_id)
NMDS

#renaming columns
#(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS, old=c("MDS1", "MDS2", "origin_site", "destination"), new=c("NMDS1","NMDS2", "Original_site", "Destination_site"))
NMDS

# re-order the factor levels before the plot
NMDS$Original_site
NMDS$Original_site <- factor(NMDS$Original_site, levels=c("Pruth Bay", "Pruth Pocket", "Choked Wolf", "Choked Sand Spit"))

# plot NMDS - Figure 6A
p2 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, shape = Destination_site, color = Original_site)) +
  geom_point(size = 7) +
  theme.nmds() +
  scale_shape_manual(values = c(16, 79, 15, 17)) +
  scale_colour_manual(values=c("dodgerblue", "blue4", "tomato", "darkred"))

p2

########################################################
### STERILE ORIGIN AND DESTINATION DAY 5 - FIGURE 6B ###
########################################################

# final sterile only
final <- ("final_eelgrass_swab")
sterile <- ("Sterilized seagrass")
s_f <- metadata_NMDS %>% 
  dplyr::filter(category %in% final) %>% 
  dplyr::filter(treatment %in% sterile)

### Creating an object to store abundances only (remove first 6 columns with dplyr)
abundances_NMDS <- s_f %>% 
  dplyr::select(-c(1:13, 1688:1689))

#Get MDS stats using matrix with only numericals (without sites' column) 
set.seed(2)
NMDS.s_f <- metaMDS(log(abundances_NMDS+1), distance = "bray", k=2)  
NMDS.s_f

stressplot(NMDS.s_f)
plot(NMDS.s_f)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.s_f$points[,1]
MDS2 = NMDS.s_f$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, origin_site = s_f$origin_site, destination = s_f$destination, sample_id = s_f$sample_id)
NMDS


#renaming columns
#(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS, old=c("MDS1", "MDS2", "origin_site", "destination"), new=c("NMDS1","NMDS2", "Original_site", "Destination_site"))
NMDS

# re-order the factor levels before the plot
NMDS$Original_site
NMDS$Original_site <- factor(NMDS$Original_site, levels=c("Pruth Bay", "Pruth Pocket", "Choked Wolf", "Choked Sand Spit"))

# NMDS plot - FIGURE 6B
p2 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, shape = Destination_site, color = Original_site)) +
  geom_point(size = 7) +
  theme.nmds() +
  scale_shape_manual(values = c(16, 79, 15, 17)) +
  scale_colour_manual(values=c("dodgerblue", "blue4", "tomato", "darkred"))

p2

###########################################################
### ARTIFICIAL ORIGIN AND DESTINATION DAY 5 - FIGURE 6C ###
###########################################################

# final artificial only
final <- ("final_eelgrass_swab")
asu <- ("Artificial seagrass")
asu_f <- metadata_NMDS %>% 
  dplyr::filter(category %in% final) %>% 
  dplyr::filter(treatment %in% asu)

### Creating an object to store abundances only (remove first 6 columns with dplyr)
abundances_NMDS <- asu_f %>% 
  dplyr::select(-c(1:13, 1688:1689))

#Get MDS stats using matrix with only numericals (without sites' column) 
set.seed(2)
NMDS.asu_f <- metaMDS(log(abundances_NMDS+1), distance = "bray", k=2)  
NMDS.asu_f

stressplot(NMDS.asu_f)
plot(NMDS.asu_f)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.asu_f$points[,1]
MDS2 = NMDS.asu_f$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, origin_site = asu_f$origin_site, destination = asu_f$destination, sample_id = asu_f$sample_id)
NMDS

#renaming columns
#(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS, old=c("MDS1", "MDS2", "origin_site", "destination"), new=c("NMDS1","NMDS2", "Original_site", "Destination_site"))
NMDS

# re-order the factor levels before the plot
NMDS$Original_site <- factor(NMDS$Original_site, levels=c("Pruth Bay", "Pruth Pocket", "Choked Wolf", "Choked Sand Spit"))
NMDS$Original_site

# plot NMDS - FIGURE 6C
p2 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, shape = Destination_site, color = Original_site)) +
  geom_point(size = 7) +
  theme.nmds() +
  scale_shape_manual(values = c(16, 79, 15, 17)) +
  scale_colour_manual(values=c("dodgerblue", "blue4", "tomato", "darkred"))

p2


#######################################################################################
#######################################################################################
#######################################################################################


#######################################################################
### FINAL NATURAL, ARTIFICIAL, STERILE, AND WATER DAY 5 - FIGURE S2 ###
#######################################################################

###select final swab samples
final <- c("final_eelgrass_swab")
final_only <- metadata_NMDS %>% 
  dplyr::filter(category %in% final )

###select final control, sterile, asu, water
all_final <- c("Sterilized seagrass", "Artificial seagrass", "Live seagrass")
final_seagrasses <- final_only %>% 
  dplyr::filter(treatment %in% all_final)


# selecting final water only
final_water <- c("water_final")
final_water_only <- metadata_NMDS %>% 
  dplyr::filter(category %in% final_water)

# combining all 4 treatments
acs <- bind_rows(list(final_seagrasses, final_water_only))

### Creating an object to store abundances only (remove first 6 columns with dplyr)
abundances_NMDS <- acs %>% 
  dplyr::select(-c(1:13, 1688:1689))

#Get MDS stats using matrix with only numericals (without sites' column) 
set.seed(2)
NMDS.acs <- metaMDS(log(abundances_NMDS+1), distance = "bray", k=2)  
NMDS.acs

stressplot(NMDS.acs)
plot(NMDS.acs)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.acs$points[,1]
MDS2 = NMDS.acs$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Destination_site = acs$destination, Treatment = acs$treatment)
NMDS

#renaming columns
#setnames(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS, old=c("MDS1", "MDS2", "Destination_site", "Treatment"), new=c("NMDS1","NMDS2", "Destination_site", "Treatment"))
NMDS

# re-order the factor levels before the plot
NMDS$Destination_site
NMDS$Destination_site <- factor(NMDS$Destination_site, levels=c("Pruth Bay", "Pruth Pocket", "Choked Wolf", "Choked Sand Spit"))
NMDS$Treatment
NMDS$Treatment <- factor(NMDS$Treatment, levels=c("Live seagrass", "Sterilized seagrass", "Artificial seagrass", "Water"))


# NMDS plot - Figure S2
p2 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, shape = Treatment, colour=Destination_site)) +
  geom_point(size = 5, alpha = 0.8) +
  theme.nmds() +
  scale_shape_manual(values = c(16, 15, 17, 8)) +
  #labs(title= "Final Water, ASU, live, sterile") +
  scale_colour_manual(values=c("dodgerblue", "blue4", "tomato", "darkred")) #+
#stat_ellipse(linetype = 1)

p2

#######################################################################################
#######################################################################################
#######################################################################################


##########################################################################
### INITIAL NATURAL AND FINAL NATURAL (REPLANT + TRANSPLANT) - FIGURE S3 #
##########################################################################

# natural only
cont <- ("Live seagrass")
cont2 <- (metadata_NMDS) %>% 
  dplyr::filter(treatment %in% cont)

### Creating an object to store abundances only (remove first 6 columns with dplyr)
abundances_NMDS <- cont2 %>% 
  dplyr::select(-c(1:13, 1688:1689))

#Get MDS stats using matrix with only numericals (without sites' column) 
set.seed(2)
NMDS.cont2 <- metaMDS(log(abundances_NMDS+1), distance = "bray", k=2)  
NMDS.cont2

stressplot(NMDS.cont2)
plot(NMDS.cont2)

#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS.cont2$points[,1]
MDS2 = NMDS.cont2$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, origin_site = cont2$origin_site, transplant_type = cont2$transplant_type_2, sample_id = cont2$sample_id)
NMDS

#renaming columns
#(data, old=c("old_name","another_old_name"), new=c("new_name", "another_new_name"))
library(data.table)
setnames(NMDS, old=c("MDS1", "MDS2", "origin_site", "transplant_type"), new=c("NMDS1","NMDS2", "Original_site", "Transplant_type"))
NMDS

# re-order the factor levels before the plot
NMDS$Original_site
NMDS$Original_site <- factor(NMDS$Original_site, levels=c("Pruth Bay", "Pruth Pocket", "Choked Wolf", "Choked Sand Spit"))
NMDS$Transplant_type
NMDS$Transplant_type <- factor(NMDS$Transplant_type, levels=c("Initial (Day 1)", "Replant (Day 5)", "Transplant same region (Day 5)","Transplant different region (Day 5)"))


# plot NMDS - Figure S3
p2 <- ggplot(NMDS, aes(x=NMDS1, y=NMDS2, shape = Transplant_type, color = Original_site)) +
  geom_point(size = 6, alpha = 0.9) +
  theme.nmds() +
  scale_shape_manual(values = c(16, 79, 15, 17)) +
  scale_colour_manual(values=c("dodgerblue", "blue4", "tomato", "darkred"))

p2


