# Project: Seagrass (Zostera marina) transplant experiment reveals core microbiome and resistance to environmental change 
# Code by: Emily Adamczyk
# Purpose: Seagrass growth plot (FIGURE S4) and stats

# libraries
library(tidyverse)
library(vegan)

# setting working directory
setwd("/Users/Emily/Documents/seagrass_transplant_experiment/data_analysis")

# reading in data
data <- read.csv("pin_prick_data_second_blade_handpicked_20200601.csv")

# renaming origin_site
data1 <- data %>% 
  dplyr::rename(origin_site = ï..origin_site)

# converting factors to characters
data1$origin_site <- as.character(data1$origin_site)
data1$destination <- as.character(data1$destination)

# *** Note: Eelgrass "control" and "live" treatments are synonymous with "natural" (all refer to natural eelgrass)

##############################################################################
##############################################################################
##############################################################################

# making a column for replant and transplant
data1$transplant_type <- NA

# replant
data1$transplant_type[as.character(data1$origin_site) == as.character(data1$destination)] <- "replant"

# transplant
data1$transplant_type[as.character(data1$origin_site) != as.character(data1$destination)] <- "transplant"


##############################################################################
##############################################################################
##############################################################################

# making a ggplot theme
theme.ggplot <- function(){
  theme_bw()+
    theme (axis.title.x = element_text(size=20, margin = margin(t = 10, r = 0, b = 0, l = 0)), #font size of x title
           axis.title.y = element_text(size=20, margin = margin(t = 0, r = 10, b = 0, l = 0)), #font size of y title
           axis.text.x = element_text(angle = 45, hjust = 1),
           axis.text = element_text(size = 16), #font size of numbers in axis
           panel.grid.major = element_blank(), #remove major grid
           panel.grid.minor = element_blank(), #remove minor grid
           axis.line = element_line(colour = "black"), #draw line in the axis
           panel.border = element_blank(), #remove lines outside the graph
           legend.title=element_text(size = 16), #remove legend title
           #legend.direction = "vertical", #direction
           #legend.justification = c(1, 1), 
           legend.position = "right", #legend is top right
           legend.key.size = unit(2.0, 'lines'), #spacing between legends
           legend.text = element_text(size = 18), #font size of legend
           plot.title = element_text(hjust = 0.5, size = 20),  #center plot title and set font size
           strip.text.x = element_text(size = 18),
           strip.background = element_rect(fill = "white")) # panel titles
}


##############################################################################
##############################################################################
##############################################################################


# creating pruth and choked regions
data1$Origin_region <- NA

# pruth origin_site
data1$Origin_region[as.character(data1$origin_site) == "pruth_bay"] <- "Pruth region"
data1$Origin_region[as.character(data1$origin_site) == "pruth_pocket"] <- "Pruth region"

# choked origin_site
data1$Origin_region[as.character(data1$origin_site) == "wolf"] <- "Choked Passage region"
data1$Origin_region[as.character(data1$origin_site) == "sand_spit"] <- "Choked Passage region"

data1$Transplanted_region <- NA

# pruth destination
data1$Transplanted_region[as.character(data1$destination) == "pruth_bay"] <- "Pruth region"
data1$Transplanted_region[as.character(data1$destination) == "pruth_pocket"] <- "Pruth region"

# choked destination
data1$Transplanted_region[as.character(data1$destination) == "wolf"] <- "Choked Passage region"
data1$Transplanted_region[as.character(data1$destination) == "sand_spit"] <- "Choked Passage region"


##############################################################################
##############################################################################
##############################################################################

###########################################
# GROWTH PLOT BASED ON REGION - FIGURE S4 #
###########################################

# making treatments panels #
data1$treatment[as.character(data1$treatment) == "control"] <- "Natural seagrass"
data1$treatment[as.character(data1$treatment) == "sterile"] <- "Sterilized seagrass"

# factors
data1$treatment = factor(data1$treatment, levels=c("Natural seagrass", "Sterilized seagrass"))

ggplot(data1, aes(x = Transplanted_region, y = second_youngest_blade_growth, fill = Origin_region)) +
  facet_wrap(. ~ treatment, scale="free") +   
  stat_boxplot(geom ='errorbar') + 
  geom_boxplot() +
  scale_y_continuous(limits = c(1, 30)) +
  scale_fill_manual(values=(c("indianred2", "dodgerblue")), 
                    labels = c(expression("Choked Passage","Pruth"))) +
  labs(y = "Second youngest leaf growth (cm)") +
  theme.ggplot()

##############################################################################
##############################################################################
##############################################################################

##############################
# SEAGRASS GROWTH STATISTICS #
##############################

#########
# ANOVA #
#########

# comparing growth of natural vs sterile shoot, regardless of region
# using one=way anova b/c i have univariate data (growth rates are only in 1 column)
shoots <- lm(second_youngest_blade_growth~treatment, data = data1)
anova(shoots)

#########################
# MEANS, SD, AND T-TEST #
#########################

# NATURAL SHOOTS

# looking at growth for natural shoots from Choked region only
data1 %>% 
  dplyr::filter(Origin_region == "Choked Passage region") %>% 
  filter(treatment == "Natural seagrass") %>% 
  summarize(mean_choked_growth = mean(second_youngest_blade_growth), choked_sd = sd(second_youngest_blade_growth))

# looking at growth for natural shoots from Pruth only
data1 %>% 
  dplyr::filter(Origin_region == "Pruth region") %>% 
  filter(treatment == "Natural seagrass") %>% 
  summarize(mean_pruth_growth = mean(second_youngest_blade_growth), pruth_sd = sd(second_youngest_blade_growth))


natural <- data1 %>% 
  filter(treatment == "Natural seagrass")
t.test(second_youngest_blade_growth~Origin_region, data = natural)


# STERILIZED/KILLED SHOOTS

# looking at growth for sterilized/killed shoots from Choked region only
data1 %>% 
  dplyr::filter(Origin_region == "Choked Passage region") %>% 
  filter(treatment == "Sterilized seagrass") %>% 
  summarize(mean_choked_growth = mean(second_youngest_blade_growth), choked_sd = sd(second_youngest_blade_growth))

data1 %>% 
  dplyr::filter(Origin_region == "Pruth region") %>% 
  filter(treatment == "Sterilized seagrass") %>% 
  summarize(mean_pruth_growth = mean(second_youngest_blade_growth), pruth_sd = sd(second_youngest_blade_growth))

sterilized <- data1 %>% 
  filter(treatment == "Sterilized seagrass")
t.test(second_youngest_blade_growth~Origin_region, data = sterilized)

