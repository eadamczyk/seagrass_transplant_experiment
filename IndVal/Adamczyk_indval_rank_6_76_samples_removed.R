# Project: Seagrass (Zostera marina) transplant experiment reveals core microbiome and resistance to environmental change 
# Code by: Emily Adamczyk
# Purpose: IndVal (indicator species analysis) for Rank 6 abundance for ALL SAMPLES (Table 3 - supplementary analyses document)


#### Libraries ####
library(tidyverse)
library(vegan)
#install.packages("labdsv")
library(labdsv)
#install.packages("indicspecies")
library(indicspecies)
library(data.table)

# PREPPING FOR RANK 6 #

#set working directory
setwd("/Users/Emily/Documents/seagrass_transplant_experiment/data_analysis/")

## INDVAL ##
# USING NON-RAREFIED RANK 6 DATA!!!

# reading in csv files
seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.otu <- read.csv("seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.otu.csv")
seqtab.nosingletons.nochim.1 <- fread("sequence_table.16s.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam <- read.csv("seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam.csv") 


# merging sam and otu tables together
names(seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam)
names(seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.otu)

# renaming sample_id
seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.otu <- seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.otu %>% 
  dplyr::rename(sample_id = X)

# removing unnecessary columns
seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam <- seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam[-c(1, 3:13, 15:16)]

# *** Note: Eelgrass "control" and "live" treatments are synonymous with "natural" (all refer to natural eelgrass)
# *** Note: Eelgrass "asu" (artificial seagrass unit) treatment is the same as artificial eelgrass

#######################################################################################
#######################################################################################
#######################################################################################

## creating transplant column ##
names(seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam)
seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$transplant_type = NA

# initial
seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$transplant_type[seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$treatment == "control" & seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$category == "initial_eelgrass_swab"] <- "initial"

# replant
seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$transplant_type[as.character(seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$origin_site) == as.character(seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$destination) & seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$treatment == "control" & seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$category == "final_eelgrass_swab"] <- "replant"

# transplant
seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$transplant_type[as.character(seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$origin_site) != as.character(seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$destination) & seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$treatment == "control" & seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam$category == "final_eelgrass_swab"] <- "transplant"


#######################################################################################
#######################################################################################
#######################################################################################

# joining otu and sam together
meta_otu <- left_join(seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.sam, seagrass_transplant_16s_NOT_RAREFIED_RANK_6_76removed.otu, by = "sample_id")

# setting as factor
meta_otu$transplant_type <- as.factor(meta_otu$transplant_type)

#######################################################################################
#######################################################################################
#######################################################################################


# RENAMING WATER VARIABLES #

# Pruth Bay
meta_otu$destination[as.character(meta_otu$origin_site) == "pruth_bay" & as.character(meta_otu$destination) == "na"] <- "pruth_bay"
# Pruth Pocket
meta_otu$destination[as.character(meta_otu$origin_site) == "pruth_pocket" & as.character(meta_otu$destination) == "na"] <- "pruth_pocket"
# Choked Sand Spit
meta_otu$destination[as.character(meta_otu$origin_site) == "sandspit" & as.character(meta_otu$destination) == "na"] <- "sandspit"
# Choked Wolf
meta_otu$destination[as.character(meta_otu$origin_site) == "wolf" & as.character(meta_otu$destination) == "na"] <- "wolf"
# water final
meta_otu$category[as.character(meta_otu$treatment) == "recovery" & as.character(meta_otu$category) == "water_filters_(0.22um)"] <- "water_final"
# water initial
meta_otu$category[as.character(meta_otu$treatment) == "deployment" & as.character(meta_otu$category) == "water_filters_(0.22um)"] <- "water_initial"
# water name
meta_otu <- meta_otu %>% 
  dplyr::mutate(treatment=recode(treatment, "recovery" = "water"),
                treatment=recode(treatment, "deployment" = "water"))

#######################################################################################
#######################################################################################
#######################################################################################

####################################
##          INDVAL                ##
# TABLE 3 (Supplementary analyses) #
####################################

##################################################################################################
### selecting for FINAL NATURAL AND FINAL WATER ONLY
##################################################################################################

natural_object <- c("control") # creating object to filter out more than 1 treatment (row) at a time

natural_only <- meta_otu  %>%    # filtering for natural only
  filter(treatment %in% natural_object) %>% 
  droplevels 

final_samples <- c("final_eelgrass_swab")

natural_final <- natural_only %>%    # filtering for final swabs only
  filter(category %in% final_samples) %>% 
  droplevels 

water <- "water_final"

water_final <- meta_otu  %>%    # filtering for final water only
  filter(category %in% water) %>% 
  droplevels 

# combining final natural and final water
nw <- bind_rows(list(natural_final, water_final))

### Creating an object to store abundances only (remove first 26 columns with dplyr)
names(nw)
nw_otu <- nw %>%
  dplyr::select(-(1:7)) #sometimes we need to use dplyr:: before functions because it crashes if used with other packages

### Factors
treatment_nw<-as.factor(nw$treatment)
class(treatment_nw)
levels(treatment_nw)


#### Multipatt analysis: final swabs for ASU and control treatments ####
multipatt.nw <- multipatt(nw_otu, nw$treatment, control = how(nperm=999))
summary(multipatt.nw)
out_z <- capture.output(summary(multipatt.nw, indvalcomp=TRUE))
out_z # A = specificity, B = fidelity



##################################################################################################
### selecting for FINAL STERILE AND FINAL WATER ONLY
##################################################################################################

sterilized_object <- c("sterile") # creating object to filter out more than 1 treatment (row) at a time

sterilized_only <- meta_otu  %>%    # filtering for sterilized only
  filter(treatment %in% sterilized_object) %>% 
  droplevels 

final_samples <- c("final_eelgrass_swab")

sterilized_final <- sterilized_only %>%    # filtering for final swabs only
  filter(category %in% final_samples) %>% 
  droplevels 


water <- "water_final"

water_final <- meta_otu  %>%    # filtering for final water only
  filter(category %in% water) %>% 
  droplevels 

# combining final sterilized and final water
sw <- bind_rows(list(sterilized_final, water_final))

### Creating an object to store abundances only (remove first 26 columns with dplyr)
names(sw)
sw_otu <- sw %>%
  dplyr::select(-(1:7)) #sometimes we need to use dplyr:: before functions because it crashes if used with other packages

### Factors
treatment_sw<-as.factor(sw$treatment)
class(treatment_sw)
levels(treatment_sw)


#### Multipatt analysis: final swabs for ASU and control treatments ####
multipatt.sw <- multipatt(sw_otu, sw$treatment, control = how(nperm=999))
summary(multipatt.sw)
out_z <- capture.output(summary(multipatt.sw, indvalcomp=TRUE))
out_z # A = specificity, B = fidelity


##################################################################################################
### selecting for FINAL ARTIFICIAL AND FINAL WATER ONLY
##################################################################################################

asu_object <- c("asu") # creating object to filter out more than 1 treatment (row) at a time

asu_only <- meta_otu  %>%    # filtering for asu only
  filter(treatment %in% asu_object) %>% 
  droplevels 

final_samples <- c("final_eelgrass_swab")

asu_final <- asu_only %>%    # filtering for final swabs only
  filter(category %in% final_samples) %>% 
  droplevels 

water <- "water_final"

water_final <- meta_otu  %>%    # filtering for final water only
  filter(category %in% water) %>% 
  droplevels 

# combining final asu and final water
asu_w <- bind_rows(list(asu_final, water_final))

### Creating an object to store abundances only (remove first 26 columns with dplyr)
names(asu_w)
asu_w_otu <- asu_w %>%
  dplyr::select(-(1:7)) #sometimes we need to use dplyr:: before functions because it crashes if used with other packages

### Factors
treatment_asu_w<-as.factor(asu_w$treatment)
class(treatment_asu_w)
levels(treatment_asu_w)


#### Multipatt analysis: final swabs for ASU and control treatments ####
multipatt.asu_w <- multipatt(asu_w_otu, asu_w$treatment, control = how(nperm=999))
summary(multipatt.asu_w)
out_z <- capture.output(summary(multipatt.asu_w, indvalcomp=TRUE))
out_z # A = specificity, B = fidelity


##################################################################################################
### selecting for FINAL ARTIFICIAL AND FINAL NATURAL ONLY
##################################################################################################

asu_live_object <- c("asu", "control") # creating object to filter out more than 1 treatment (row) at a time

asu_live <- meta_otu  %>%    # filtering for asu only
  filter(treatment %in% asu_live_object) %>% 
  droplevels 

final_samples <- c("final_eelgrass_swab")

asu_live_final <- asu_live %>%    # filtering for final swabs only
  filter(category %in% final_samples) %>% 
  droplevels 

### Creating an object to store abundances only (remove first 26 columns with dplyr)
names(asu_live_final)
asu_live_final_otu <- asu_live_final %>%
  dplyr::select(-(1:7)) #sometimes we need to use dplyr:: before functions because it crashes if used with other packages

### Factors
treatment_asu_live<-as.factor(asu_live_final$treatment)
class(treatment_asu_live)
levels(treatment_asu_live)


#### Multipatt analysis: final swabs for ASU and control treatments ####
multipatt.asu_live <- multipatt(asu_live_final_otu, asu_live_final$treatment, control = how(nperm=999))
summary(multipatt.asu_live)
out_z <- capture.output(summary(multipatt.asu_live, indvalcomp=TRUE))
out_z # A = specificity, B = fidelity


##################################################################################################
### selecting for FINAL ARTIFICIAL AND FINAL STERILE ONLY
##################################################################################################

asu_sterile_object <- c("asu", "sterile") # creating object to filter out more than 1 treatment (row) at a time

asu_sterile <- meta_otu  %>%    # filtering for asu only
  filter(treatment %in% asu_sterile_object) %>% 
  droplevels 

final_samples <- c("final_eelgrass_swab")

asu_sterile_final <- asu_sterile %>%    # filtering for final swabs only
  filter(category %in% final_samples) %>% 
  droplevels 

### Creating an object to store abundances only (remove first 26 columns with dplyr)
names(asu_sterile_final)
asu_sterile_final_otu <- asu_sterile_final %>%
  dplyr::select(-(1:7)) #sometimes we need to use dplyr:: before functions because it crashes if used with other packages

### Factors
treatment_asu_sterile<-as.factor(asu_sterile_final$treatment)
class(treatment_asu_sterile)
levels(treatment_asu_sterile)


#### Multipatt analysis: final swabs for ASU and control treatments ####
multipatt.asu_sterile <- multipatt(asu_sterile_final_otu, asu_sterile_final$treatment, control = how(nperm=999))
summary(multipatt.asu_sterile)
out_z <- capture.output(summary(multipatt.asu_sterile, indvalcomp=TRUE))
out_z # A = specificity, B = fidelity

