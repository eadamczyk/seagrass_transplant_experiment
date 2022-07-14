# Project: Seagrass (Zostera marina) transplant experiment reveals core microbiome and resistance to environmental change 
# Code by: Emily Adamczyk
# Purpose: PERMANOVAs using all samples


####Libraries####
library(vegan)
library(readr)
library(reshape2)
library(dplyr)
library(devtools)
#devtools::install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(remotes)


#set working directory
setwd("/Users/Emily/Documents/seagrass_transplant_experiment/data_analysis/")

###########################  PERMANOVAS USING BRAY CURTIS  ############################


###Read table with labels of site or host type *this table also contain abundances
meta_otu <- read.csv("metadata_PERMANOVA_all_samples_20220427.csv", header=T)

#View(meta_otu)
meta_otu <- meta_otu[,-c(1,2)] #removing 2 unnecessary columns
names(meta_otu)

###########################################################################################
###########################################################################################
### making replants and transplant same region rows for natural, sterile and artificial ###
###########################################################################################
###########################################################################################


# note: throughout script, "natural" eelgrass samples are synonymous with "live" and "control" samples
# note: throughout script, "ASUs" (artificial seagrass units) are synonymous with "artificial" samples



### making replants and transplant same region rows for control (natural), sterile and asu (aritifical)

meta_otu$transplant_type <- NA # new column

######################
## natural eelgrass ##
######################

# live_eelgrass replant
meta_otu$transplant_type[as.character(meta_otu$origin_site) == as.character(meta_otu$destination) & meta_otu$treatment == "live_eelgrass" & meta_otu$category == "final_eelgrass_swab"] <- "live_eelgrass_replant"

### live_eelgrass transplant same region ###
# transplant same region (pruth bay to pruth pocket)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "pruth_bay" & as.character(meta_otu$destination) == "pruth_pocket" & meta_otu$treatment == "live_eelgrass" & meta_otu$category == "final_eelgrass_swab"] <- "live_eelgrass_transplant_same_region"

# transplant same region (pruth pocket to pruth bay)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "pruth_pocket" & as.character(meta_otu$destination) == "pruth_bay" & meta_otu$treatment == "live_eelgrass" & meta_otu$category == "final_eelgrass_swab"] <- "live_eelgrass_transplant_same_region"

# transplant same region (choked choked_wolf to choked sand spit)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "choked_wolf" & as.character(meta_otu$destination) == "choked_sand_spit" & meta_otu$treatment == "live_eelgrass" & meta_otu$category == "final_eelgrass_swab"] <- "live_eelgrass_transplant_same_region"

# transplant same region (choked sand spit to choked choked_wolf)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "choked_sand_spit" & as.character(meta_otu$destination) == "choked_wolf" & meta_otu$treatment == "live_eelgrass" & meta_otu$category == "final_eelgrass_swab"] <- "live_eelgrass_transplant_same_region"

############
### ASUs ###
############

# ASU replant
meta_otu$transplant_type[as.character(meta_otu$origin_site) == as.character(meta_otu$destination) & meta_otu$treatment == "asu" & meta_otu$category == "final_eelgrass_swab"] <- "asu_replant"

### asu transplant same region ###
# transplant same region (pruth bay to pruth pocket)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "pruth_bay" & as.character(meta_otu$destination) == "pruth_pocket" & meta_otu$treatment == "asu" & meta_otu$category == "final_eelgrass_swab"] <- "asu_transplant_same_region"

# transplant same region (pruth pocket to pruth bay)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "pruth_pocket" & as.character(meta_otu$destination) == "pruth_bay" & meta_otu$treatment == "asu" & meta_otu$category == "final_eelgrass_swab"] <- "asu_transplant_same_region"

# transplant same region (choked choked_wolf to choked sand spit)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "choked_wolf" & as.character(meta_otu$destination) == "choked_sand_spit" & meta_otu$treatment == "asu" & meta_otu$category == "final_eelgrass_swab"] <- "asu_transplant_same_region"

# transplant same region (choked sand spit to choked choked_wolf)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "choked_sand_spit" & as.character(meta_otu$destination) == "choked_wolf" & meta_otu$treatment == "asu" & meta_otu$category == "final_eelgrass_swab"] <- "asu_transplant_same_region"

#############
## STERILE ##
#############

# sterile_eelgrass replant
meta_otu$transplant_type[as.character(meta_otu$origin_site) == as.character(meta_otu$destination) & meta_otu$treatment == "sterile_eelgrass" & meta_otu$category == "final_eelgrass_swab"] <- "sterile_eelgrass_replant"

### sterile_eelgrass transplant same region ###
# transplant same region (pruth bay to pruth pocket)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "pruth_bay" & as.character(meta_otu$destination) == "pruth_pocket" & meta_otu$treatment == "sterile_eelgrass" & meta_otu$category == "final_eelgrass_swab"] <- "sterile_eelgrass_transplant_same_region"

# transplant same region (pruth pocket to pruth bay)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "pruth_pocket" & as.character(meta_otu$destination) == "pruth_bay" & meta_otu$treatment == "sterile_eelgrass" & meta_otu$category == "final_eelgrass_swab"] <- "sterile_eelgrass_transplant_same_region"

# transplant same region (choked choked_wolf to choked sand spit)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "choked_wolf" & as.character(meta_otu$destination) == "choked_sand_spit" & meta_otu$treatment == "sterile_eelgrass" & meta_otu$category == "final_eelgrass_swab"] <- "sterile_eelgrass_transplant_same_region"

# transplant same region (choked sand spit to choked choked_wolf)
meta_otu$transplant_type[as.character(meta_otu$origin_site) == "choked_sand_spit" & as.character(meta_otu$destination) == "choked_wolf" & meta_otu$treatment == "sterile_eelgrass" & meta_otu$category == "final_eelgrass_swab"] <- "sterile_eelgrass_transplant_same_region"

###################################################################################
###################################################################################

### making column for regions to compare among regions
meta_otu$new_regions <- NA

# Pruth region only
meta_otu$new_regions[as.character(meta_otu$destination) == "pruth_bay"] <- "pruth_region_only"
meta_otu$new_regions[as.character(meta_otu$destination) == "pruth_pocket"] <- "pruth_region_only"

# Choked region only
meta_otu$new_regions[as.character(meta_otu$destination) == "choked_wolf"] <- "choked_region_only"
meta_otu$new_regions[as.character(meta_otu$destination) == "choked_sand_spit"] <- "choked_region_only"


###################################################################################
###################################################################################
###########################  PERMANOVAS USING BRAY CURTIS  ########################
###################################################################################
###################################################################################

###########
# TABLE 1 #
###########

#################
# NATURAL Day 1 #
#################

#Filtering for natural (live) initial only among sites
meta_otu_c <- meta_otu %>% 
  filter(treatment == "live_eelgrass")
meta_otu_ic <-  meta_otu_c %>% 
  filter(category == "initial_eelgrass_swab")

### Create table containing only ASVs
bact_abund_ic <- meta_otu_ic %>%  
  dplyr::select(-c(1:12, 2457:2460))
names(bact_abund_ic)

###LOG-transformation
bact_log_ic <- log1p(bact_abund_ic)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.ic<-vegdist(bact_log_ic,m="bray")
bact.dist.ic

###Test for homogeneity in dispersions
homogeneity <-betadisper(bact.dist.ic,meta_otu_ic$origin_site)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)


#########################################################################
# PERMANOVA for "initial microbiota on natural seagrass leaves (Day 1)" #
#########################################################################

# comparing microbial communities among all sites on Day 1 for natural eelgrass
adonis2(bact.dist.ic~origin_site, data=meta_otu_ic, permutations = 999)

#######################################
# Pairwise PERMANOVA for TABLE 1 CODE #
#######################################

# pairwise comparison of microbial communities among all sites on Day 1 for natural eelgrass
pairwise.adonis(x = bact.dist.ic,
                factors = meta_otu_ic$origin_site, 
                sim.function="vegdist",
                sim.method="bray", 
                p.adjust.m="hoch")


###################################################################################
###################################################################################
###################################################################################

############
# TABLE S6 #
############

################################################################################
# artificial (asu), final sterile, final control (natural), and water on Day 5 #
################################################################################

# filtering for final artificial, sterile, natural, and water
final_swab <- meta_otu[which(meta_otu$category == "final_eelgrass_swab"),]
asu <- final_swab[which(final_swab$treatment == "asu"),]
control <- final_swab[which(final_swab$treatment == "live_eelgrass"),]
sterile <- final_swab[which(final_swab$treatment == "sterile_eelgrass"),]
water <- meta_otu[which(meta_otu$category == "recovery"),]

# making new df
asu_sterile_control <- rbind(asu, control, sterile, water)

### Create table containing only ASVs
bact_abund_asu_sterile_control <- asu_sterile_control %>% select(-c(1:12, 2457:2460))

###LOG-transformation
bact_log_asc <- log1p(bact_abund_asu_sterile_control)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.asc <- vegdist(bact_log_asc, m="bray")

###Test for homogeneity in dispersions 
homogeneity <-betadisper(bact.dist.asc, asu_sterile_control$destination)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)

#################
# TABLE S6 CODE #
#################

# pairwise PERMANOVA comparing all 3 eelgrass treatments and water on Day 5
pairwise.adonis(x = bact.dist.asc,
                factors = asu_sterile_control$treatment, 
                sim.function="vegdist", 
                sim.method="bray", 
                p.adjust.m="hoch")


###################################################################################
###################################################################################
###################################################################################

############
# Table S7 #
############

######################################################################################
# artificial (asu), final sterile, and final control (natural) on Day 5 - all sites  #
######################################################################################

# filtering for final artificial, sterile, natural
final_swab <- meta_otu[which(meta_otu$category == "final_eelgrass_swab"),]
asu <- final_swab[which(final_swab$treatment == "asu"),]
control <- final_swab[which(final_swab$treatment == "live_eelgrass"),]
sterile <- final_swab[which(final_swab$treatment == "sterile_eelgrass"),]

# making new df
asu_sterile_control_all_sites <- rbind(asu, control, sterile)

### Create table containing only ASVs
bact_abund_asu_sterile_control_all_sites <- asu_sterile_control_all_sites %>% select(-c(1:12, 2457:2460))

###LOG-transformation
bact_log_asc_all_sites <- log1p(bact_abund_asu_sterile_control_all_sites)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.asc_all_sites <- vegdist(bact_log_asc_all_sites, m="bray")

###Test for homogeneity in dispersions 
homogeneity <-betadisper(bact.dist.asc_all_sites, asu_sterile_control_all_sites$destination)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)


# PERMANOVA comparing all 3 eelgrass treatments on Day 5
adonis2(bact.dist.asc_all_sites~destination, data=asu_sterile_control_all_sites, permutations = 999)


#############################
# TABLE S7 CODE - ALL SITES #
#############################

# pairwise PERMANOVA comparing all 3 eelgrass treatments on Day 5
pairwise.adonis(x = bact.dist.asc_all_sites,
                factors = asu_sterile_control_all_sites$treatment, 
                sim.function="vegdist", 
                sim.method="bray", 
                p.adjust.m="hoch")



#############################################################################################
# artificial (asu), final sterile, and final control (natural) on Day 5 - Pruth region only #
#############################################################################################

#Filtering for ASU, final sterile, final control *** IN PRUTH ONLY ***
pruth <- meta_otu[which(meta_otu$new_regions == c("pruth_region_only")),]
final_swab <- pruth[which(pruth$category == "final_eelgrass_swab"),]
asu <- final_swab[which(final_swab$treatment == "asu"),]
control <- final_swab[which(final_swab$treatment == "live_eelgrass"),]
sterile <- final_swab[which(final_swab$treatment == "sterile_eelgrass"),]

# new df
asu_sterile_control_pruth <- rbind(asu, control, sterile)

### Create table containing only ASVs
bact_abund_asu_sterile_control_pruth <- asu_sterile_control_pruth %>% select(-c(1:12, 2457:2460))

###LOG-transformation
bact_log_asc_pruth <- log1p(bact_abund_asu_sterile_control_pruth)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.asc_pruth <- vegdist(bact_log_asc_pruth, m="bray")

###Test for homogeneity in dispersions = it has to be homogeneous (non-significant), it is an assumption of PERMANOVA
homogeneity <-betadisper(bact.dist.asc_pruth, asu_sterile_control_pruth$destination)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)

#############################
# TABLE S7 CODE - ALL SITES #
#############################

# pairwise PERMANOVA comparing all 3 eelgrass treatments on Day 5
pairwise.adonis(x = bact.dist.asc_pruth,
                factors = asu_sterile_control_pruth$treatment, 
                sim.function="vegdist", 
                sim.method="bray", 
                p.adjust.m="hoch")


##############################################################################################
# artificial (asu), final sterile, and final control (natural) on Day 5 - Choked region only #
##############################################################################################

#Filtering for ASU, final sterile, final control *** IN CHOKED ONLY ***
choked <- meta_otu[which(meta_otu$new_regions == c("choked_region_only")),]
final_swab <- choked[which(choked$category == "final_eelgrass_swab"),]
asu <- final_swab[which(final_swab$treatment == "asu"),]
control <- final_swab[which(final_swab$treatment == "live_eelgrass"),]
sterile <- final_swab[which(final_swab$treatment == "sterile_eelgrass"),]

# new df
asu_sterile_control_choked <- rbind(asu, control, sterile)

### Create table containing only ASVs
bact_abund_asu_sterile_control_choked <- asu_sterile_control_choked %>% select(-c(1:12, 2457:2460))

###LOG-transformation
bact_log_asc_choked <- log1p(bact_abund_asu_sterile_control_choked)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.asc_choked <- vegdist(bact_log_asc_choked, m="bray")

###Test for homogeneity in dispersions = it has to be homogeneous (non-significant), it is an assumption of PERMANOVA
homogeneity <-betadisper(bact.dist.asc_choked, asu_sterile_control_choked$destination)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)

#############################
# TABLE S7 CODE - ALL SITES #
#############################

# pairwise PERMANOVA comparing all 3 eelgrass treatments on Day 5
pairwise.adonis(x = bact.dist.asc_choked,
                factors = asu_sterile_control_choked$treatment, 
                sim.function="vegdist", 
                sim.method="bray", 
                p.adjust.m="hoch")

###################################################################################
###################################################################################
###################################################################################

###########
# TABLE 2 #
###########

########################################################
# replant and transplant to same region only for Day 5 #
#               ALL REGIONS                            #
########################################################

#Filtering for ASU (artificial), final sterile, final control (natural) 
final_swab <- meta_otu[which(meta_otu$category == "final_eelgrass_swab"),]
asu <- final_swab[which(final_swab$treatment == "asu"),]
control<- final_swab[which(final_swab$treatment == "live_eelgrass"),]
sterile <- final_swab[which(final_swab$treatment == "sterile_eelgrass"),]
replant_c <- control[which(control$transplant_type == "live_eelgrass_replant"),]
transplant_same_c <- control[which(control$transplant_type == "live_eelgrass_transplant_same_region"),]
replant_asu <- asu[which(asu$transplant_type == "asu_replant"),]
transplant_same_asu <- asu[which(asu$transplant_type == "asu_transplant_same_region"),]
replant_sterile <- sterile[which(sterile$transplant_type == "sterile_eelgrass_replant"),]
transplant_same_sterile <- sterile[which(sterile$transplant_type == "sterile_eelgrass_transplant_same_region"),]

# new df
asu_sterile_control_x <- rbind(replant_asu, transplant_same_asu, replant_sterile, transplant_same_sterile, replant_c, transplant_same_c)

### Create table containing only ASVs
bact_abund_asu_sterile_control <- asu_sterile_control_x %>% select(-c(1:12, 2457:2460))

###LOG-transformation
bact_log_asc <- log1p(bact_abund_asu_sterile_control)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.asc <- vegdist(bact_log_asc, m="bray")

###Test for homogeneity in dispersions
homogeneity <-betadisper(bact.dist.asc, asu_sterile_control_x$destination)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)

# two-factor permanova to see if treatment or destination site has more variance
adonis2(bact.dist.asc ~ treatment * destination, data=asu_sterile_control_x, permutations = 999, by_margin = T) # use by=margin otherwise it will understand as sequential and will try to explain destination first, then whatever is left from destination(residual) it would be used to explain treatment


# Pairwise PERMANOVA for comparison between treatments
pairwise.adonis(x = bact.dist.asc,
                factors = asu_sterile_control_x$treatment, 
                sim.function="vegdist", 
                sim.method="bray", 
                p.adjust.m="hoch")


###########
# TABLE 2 #
###########

########################################################
# replant and transplant to same region only for Day 5 #
#               PRUTH ONLY                             #
########################################################

#Filtering for ASU (artificial), final sterile, final control (natural) 
final_swab <- meta_otu[which(meta_otu$category == "final_eelgrass_swab"),]
asu <- final_swab[which(final_swab$treatment == "asu"),]
control<- final_swab[which(final_swab$treatment == "live_eelgrass"),]
sterile <- final_swab[which(final_swab$treatment == "sterile_eelgrass"),]
replant_c <- control[which(control$transplant_type == "live_eelgrass_replant"),]
transplant_same_c <- control[which(control$transplant_type == "live_eelgrass_transplant_same_region"),]
replant_asu <- asu[which(asu$transplant_type == "asu_replant"),]
transplant_same_asu <- asu[which(asu$transplant_type == "asu_transplant_same_region"),]
replant_sterile <- sterile[which(sterile$transplant_type == "sterile_eelgrass_replant"),]
transplant_same_sterile <- sterile[which(sterile$transplant_type == "sterile_eelgrass_transplant_same_region"),]

# new df
asu_sterile_control_x <- rbind(replant_asu, transplant_same_asu, replant_sterile, transplant_same_sterile, replant_c, transplant_same_c)

# filtering for pruth only
asu_sterile_control_pruth_pocket <- asu_sterile_control_x[which(asu_sterile_control_x$destination == "pruth_pocket"),]
asu_sterile_control_pruth_bay <- asu_sterile_control_x[which(asu_sterile_control_x$destination == "pruth_bay"),]

# pruth only df
pruth_region <- rbind(asu_sterile_control_pruth_pocket, asu_sterile_control_pruth_bay)

# filtering ASVs only
bact_abund_pruth_region <- pruth_region %>% select(-c(1:12, 2457:2460))

###LOG-transformation
bact_log_pr <- log1p(bact_abund_pruth_region)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.pr <- vegdist(bact_log_pr, m="bray")

###Test for homogeneity in dispersions
homogeneity <-betadisper(bact.dist.pr, pruth_region$destination)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)

# Pairwise PERMANOVA comparison between treatments for Pruth region only
pairwise.adonis(x = bact.dist.pr,
                factors = pruth_region$treatment, 
                sim.function="vegdist", 
                sim.method="bray", 
                p.adjust.m="hoch")



###########
# TABLE 2 #
###########

########################################################
# replant and transplant to same region only for Day 5 #
#               CHOKED ONLY                            #
########################################################

#Filtering for ASU (artificial), final sterile, final control (natural) 
final_swab <- meta_otu[which(meta_otu$category == "final_eelgrass_swab"),]
asu <- final_swab[which(final_swab$treatment == "asu"),]
control<- final_swab[which(final_swab$treatment == "live_eelgrass"),]
sterile <- final_swab[which(final_swab$treatment == "sterile_eelgrass"),]
replant_c <- control[which(control$transplant_type == "live_eelgrass_replant"),]
transplant_same_c <- control[which(control$transplant_type == "live_eelgrass_transplant_same_region"),]
replant_asu <- asu[which(asu$transplant_type == "asu_replant"),]
transplant_same_asu <- asu[which(asu$transplant_type == "asu_transplant_same_region"),]
replant_sterile <- sterile[which(sterile$transplant_type == "sterile_eelgrass_replant"),]
transplant_same_sterile <- sterile[which(sterile$transplant_type == "sterile_eelgrass_transplant_same_region"),]

# combining df
asu_sterile_control_x <- rbind(replant_asu, transplant_same_asu, replant_sterile, transplant_same_sterile, replant_c, transplant_same_c)

# filtering for choked only
asu_sterile_control_choked_wolf <- asu_sterile_control_x[which(asu_sterile_control_x$destination == "choked_wolf"),]
asu_sterile_control_choked_sand_spit <- asu_sterile_control_x[which(asu_sterile_control_x$destination == "choked_sand_spit"),]

# choked only df
choked_region <- rbind(asu_sterile_control_choked_wolf, asu_sterile_control_choked_sand_spit)

# filtering for ASVs only
bact_abund_choked_region <- choked_region %>% select(-c(1:12, 2457:2460))

###LOG-transformation
bact_log_cr <- log1p(bact_abund_choked_region)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.cr <- vegdist(bact_log_cr, m="bray")

###Test for homogeneity in dispersions
homogeneity <-betadisper(bact.dist.cr, choked_region$destination)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)

# Pairwise PERMANOVA comparison between treatments for Choked region only
pairwise.adonis(x = bact.dist.cr,
                factors = choked_region$treatment, 
                sim.function="vegdist", 
                sim.method="bray", 
                p.adjust.m="hoch")





###################################################################################
###################################################################################
###################################################################################

#############################################################################
# two-way PERMANOVAs comparing the effect of original and destination sites #
#    separate PERMANOVA for each eelgrass treatment using Day 5 samples     #
#############################################################################


###########
# TABLE 4 #
###########

#################
# NATURAL Day 5 #
#################

#Filtering for natural (live) final only among sites
meta_otu_f <- meta_otu %>% 
  filter(treatment == "live_eelgrass")
meta_otu_fc <-  meta_otu_f %>% 
  filter(category == "final_eelgrass_swab")

### Create table containing only ASVs
bact_abund_fc <- meta_otu_fc %>%  
  dplyr::select(-c(1:12, 2457:2460))
names(bact_abund_fc)

###LOG-transformation
bact_log_fc <- log1p(bact_abund_fc)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.fc<-vegdist(bact_log_fc,m="bray")
bact.dist.fc

###Test for homogeneity in dispersions
homogeneity <-betadisper(bact.dist.fc,meta_otu_fc$origin_site)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)


# two-factor permanova to see if origin site or destination site has more variance
adonis2(bact.dist.fc ~ origin_site * destination, data=meta_otu_fc, permutations = 999, by_margin = T) # use by=margin otherwise it will understand as sequential and will try to explain destination first, then whatever is left from destination(residual) it would be used to explain origin site



###########
# TABLE 4 #
###########

###########################
# STERILIZED/KILLED Day 5 #
###########################

#Filtering for sterilized/killed final only among sites
meta_otu_s <- meta_otu %>% 
  filter(treatment == "sterile_eelgrass")
meta_otu_fs <-  meta_otu_s %>% 
  filter(category == "final_eelgrass_swab")

### Create table containing only ASVs
bact_abund_fs <- meta_otu_fs %>%  
  dplyr::select(-c(1:12, 2457:2460))
names(bact_abund_fs)

###LOG-transformation
bact_log_fs <- log1p(bact_abund_fs)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.fs<-vegdist(bact_log_fs,m="bray")
bact.dist.fs

###Test for homogeneity in dispersions
homogeneity <-betadisper(bact.dist.fs,meta_otu_fs$origin_site)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)


# two-factor permanova to see if origin site or destination site has more variance
adonis2(bact.dist.fs ~ origin_site * destination, data=meta_otu_fs, permutations = 999, by_margin = T) # use by=margin otherwise it will understand as sequential and will try to explain destination first, then whatever is left from destination(residual) it would be used to explain origin site



###########
# TABLE 4 #
###########

####################
# ARTIFICIAL Day 5 #
####################

#Filtering for artificial eelgrass only among sites
meta_otu_a <- meta_otu %>% 
  filter(treatment == "asu")
meta_otu_fa <-  meta_otu_a %>% 
  filter(category == "final_eelgrass_swab")

### Create table containing only ASVs
bact_abund_fa <- meta_otu_fa %>%  
  dplyr::select(-c(1:12, 2457:2460))
names(bact_abund_fa)

###LOG-transformation
bact_log_fa <- log1p(bact_abund_fa)

###  Distance matrix * this is using Bray-Curtis  ###
bact.dist.fa<-vegdist(bact_log_fa,m="bray")
bact.dist.fa

###Test for homogeneity in dispersions
homogeneity <-betadisper(bact.dist.fa,meta_otu_fa$origin_site)
permutest(homogeneity, pairwise = TRUE)
plot(homogeneity)


# two-factor permanova to see if origin site or destination site has more variance
adonis2(bact.dist.fa ~ origin_site * destination, data=meta_otu_fa, permutations = 999, by_margin = T) # use by=margin otherwise it will understand as sequential and will try to explain destination first, then whatever is left from destination(residual) it would be used to explain origin site

