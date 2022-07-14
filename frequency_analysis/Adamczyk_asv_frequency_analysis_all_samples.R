# Project: Seagrass (Zostera marina) transplant experiment reveals core microbiome and resistance to environmental change 
# Code by: Emily Adamczyk
# Purpose: ASV frequency analysis - using RANK 6 non-rarefied ASVs and ALL SAMPLES. The read threshold is >= 10


####Libraries####
library(tidyverse)

#set working directory
setwd("/Users/Emily/Documents/seagrass_transplant_experiment/data_analysis")

### reading in metadata + asv table
# NOTE: USING RANK 6 NON-RAREFIED DATA
meta_otu <- read.csv("seagrass_transplant_16s_meta_otu_NOT_RAREFIED_RANK_6_all_samples.csv", header = T)

# reading in taxa table
taxa <- read.csv("seagrass_transplant_16s_NOT_RAREFIED_RANK_6_all_samples.tax.csv", header=T)
names(taxa)[1] <- "ASV" # renaming first column to ASV

#View(meta_otu)
meta_otu <- meta_otu[,-c(1)] #removing first column
names(meta_otu)


#####################
# CODE FOR TABLE S8 #
#####################


#################################################
## NATURAL EELGRASS - ALL SITES AND TIMEPOINTS ##
#################################################

#Filtering for natural - all sites
natural <- meta_otu %>% 
  filter(treatment == "control")

# making copy of natural
natural_1 <- natural

# for loop to make reads >= 10 as 1 and <= 9 as 0
# Martin's nested loop
for(r in 1:94) {
  for(c in 13:540) {
    if(natural_1[r,c] >= 10) {
      natural_1[r,c] = 1
    }
    else{
      natural_1[r,c] = 0
    }
  }
}

# making copy of natural_1 (natural_2)
natural_2 <- natural_1

# adding new row to natural_initial_2
natural_2[95,] <- 0

# running for loop to calculate frequency
# Thank you Martin!!!
for(c in 13:540) {
  for(r in 1:94) {
    natural_2[95,c] <-  natural_2[95,c] + natural_2[r,c]
  }
  natural_2[95,c] <-  natural_2[95,c]/94
}


# saving row 95 as its own df
natural_results <- natural_2[95, 13:540]

# transposing natural_results row so that it is a single column
natural_results2 <- as.data.frame(t(natural_results))
names(natural_results2)[1] <- "natural_frequency" # renaming 95 to natural frequency

# making row names a column
natural_results3 <- natural_results2                  # making copy of df 
natural_results3$ASV <- row.names(natural_results3)   # Apply row.names function

# merging dfs
natural_eelgrass <- merge(natural_results3, taxa, by = "ASV")

# making new column for order; family; genus
natural_eelgrass$order_family_genus <- paste(natural_eelgrass$Rank4, natural_eelgrass$Rank5, natural_eelgrass$Rank6, sep = "; ")

# removing unnecessary columns
natural_eelgrass.2 <- natural_eelgrass[-c(3:10)]

# reordering column so that frequencies are descending
natural_eelgrass.3 <- natural_eelgrass.2 %>% 
  arrange(desc(natural_frequency))

natural_eelgrass.3


#########################################################
## ARTIFICIAL EELGRASS - ALL SITES AND FINAL TIMEPOINT ##
#########################################################

#Filtering for artificial and sterile - all sites
final <- meta_otu %>% 
  filter(category == "final_eelgrass_swab")
asu <- final %>% 
  filter(treatment == "asu")


# making copy of asu
asu_1 <- asu

# for loop to make reads >= 2 asu 1 and <= 2 asu 0
# Martin's nested loop
for(r in 1:46) {
  for(c in 13:540) {
    if(asu_1[r,c] >= 10) {
      asu_1[r,c] = 1
    }
    else{
      asu_1[r,c] = 0
    }
  }
}

# making copy of asu_1 (asu_2)
asu_2 <- asu_1

# adding new row to asu_initial_2
asu_2[47,] <- 0

# running for loop to calculate frequency
# Thank you Martin!!!
for(c in 13:540) {
  for(r in 1:46) {
    asu_2[47,c] <-  asu_2[47,c] + asu_2[r,c]
  }
  asu_2[47,c] <-  asu_2[47,c]/46
}


# saving row 47 asu its own df
asu_results <- asu_2[47, 13:540]

# transposing asu_results row so that it is a single column
asu_results2 <- as.data.frame(t(asu_results))
names(asu_results2)[1] <- "artificial_frequency" # renaming 47 to frequency

# making row names a column
asu_results3 <- asu_results2                  # making copy of df 
asu_results3$ASV <- row.names(asu_results3)   # Apply row.names function

# merging dfs
artificial <- merge(asu_results3, taxa, by = "ASV")

# making new column for order; family; genus
artificial$order_family_genus <- paste(artificial$Rank4, artificial$Rank5, artificial$Rank6, sep = "; ")

# removing unnecessary columns
artificial.2 <- artificial[-c(3:10)]

# reordering column so that frequencies are descending
artificial.3 <- artificial.2 %>% 
  arrange(desc(artificial_frequency))

artificial.3

######################################################
## STERILE EELGRASS - FINAL TIMEPOINT AND ALL SITES ##
######################################################

#Filtering for sterile_final - all sites
sterile_final <- meta_otu %>% 
  filter(treatment == "sterile") %>% 
  filter(category == "final_eelgrass_swab")


# making copy of sterile_final
sterile_final_1 <- sterile_final

# for loop to make reads >= 2 sterile_final 1 and <= 2 sterile_final 0
# Martin's nested loop
for(r in 1:47) {
  for(c in 13:540) {
    if(sterile_final_1[r,c] >= 10) {
      sterile_final_1[r,c] = 1
    }
    else{
      sterile_final_1[r,c] = 0
    }
  }
}

# making copy of sterile_final_1 (sterile_final_2)
sterile_final_2 <- sterile_final_1

# adding new row to sterile_final_initial_2
sterile_final_2[48,] <- 0

# running for loop to calculate frequency
# Thank you Martin!!!
for(c in 13:540) {
  for(r in 1:47) {
    sterile_final_2[48,c] <-  sterile_final_2[48,c] + sterile_final_2[r,c]
  }
  sterile_final_2[48,c] <-  sterile_final_2[48,c]/47
}


# saving row 48 sterile_final its own df
sterile_final_results <- sterile_final_2[48, 13:540]

# transposing sterile_final_results row so that it is a single column
sterile_final_results2 <- as.data.frame(t(sterile_final_results))
names(sterile_final_results2)[1] <- "sterilized_frequency" # renaming 48 to frequency

# making row names a column
sterile_final_results3 <- sterile_final_results2                  # making copy of df 
sterile_final_results3$ASV <- row.names(sterile_final_results3)   # Apply row.names function

# merging dfs
sterile_final <- merge(sterile_final_results3, taxa, by = "ASV")

# making new column for order; family; genus
sterile_final$order_family_genus <- paste(sterile_final$Rank4, sterile_final$Rank5, sterile_final$Rank6, sep = "; ")

# removing unnecessary columns
sterile_final.2 <- sterile_final[-c(3:10)]

# reordering column so that frequencies are descending
sterile_final.3 <- sterile_final.2 %>% 
  arrange(desc(sterilized_frequency))

sterile_final.3

######################################################
## WATER - ALL TIMEPOINTS AND ALL SITES ##
######################################################

#Filtering for water - all sites and timepoints
water <- meta_otu %>% 
  filter(category == "water_filters_(0.22um)")


# making copy of water
water_1 <- water

# for loop to make reads >= 2 water 1 and <= 2 water 0
# Martin's nested loop
for(r in 1:24) {
  for(c in 13:540) {
    if(water_1[r,c] >= 10) {
      water_1[r,c] = 1
    }
    else{
      water_1[r,c] = 0
    }
  }
}

# making copy of water_1 (water_2)
water_2 <- water_1

# adding new row to water_initial_2
water_2[25,] <- 0

# running for loop to calculate frequency
# Thank you Martin!!!
for(c in 13:540) {
  for(r in 1:24) {
    water_2[25,c] <-  water_2[25,c] + water_2[r,c]
  }
  water_2[25,c] <-  water_2[25,c]/24
}


# saving row 25 water its own df
water_results <- water_2[25, 13:540]

# transposing water_results row so that it is a single column
water_results2 <- as.data.frame(t(water_results))
names(water_results2)[1] <- "water_frequency" # renaming 25 to frequency

# making row names a column
water_results3 <- water_results2                  # making copy of df 
water_results3$ASV <- row.names(water_results3)   # Apply row.names function

# merging dfs
water <- merge(water_results3, taxa, by = "ASV")

# making new column for order; family; genus
water$order_family_genus <- paste(water$Rank4, water$Rank5, water$Rank6, sep = "; ")

# removing unnecessary columns
water.2 <- water[-c(3:10)]

# reordering column so that frequencies are descending
water.3 <- water.2 %>% 
  arrange(desc(water_frequency))

water.3

################################################################################
### MAKING FREQUENCY TABLE FOR ALL NATURAL EELGRASS, STERILE, ASU, AND WATER ###
#                            TABLE S8                                          #
################################################################################

# joining together natural and sterilized
natural_sterilized <- merge(natural_eelgrass.3, sterile_final.3, by = c("ASV", "order_family_genus"))

# adding artificial
natural_sterilized_artificial <- merge(natural_sterilized, artificial.3, by = c("ASV", "order_family_genus"))

# adding water
all_treatments <- merge(natural_sterilized_artificial, water.3, by = c("ASV", "order_family_genus"))

# reordering column so that frequencies are descending
all_treatments.2 <- all_treatments[order(all_treatments$natural_frequency, decreasing = TRUE),]  

# saving df as CSV
write.csv(all_treatments.2, "frequency_table_all_samples.csv") # This is Table S8
