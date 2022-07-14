# Project: Seagrass (Zostera marina) transplant experiment reveals core microbiome and resistance to environmental change 
# Code by: Emily Adamczyk
# Purpose: taxa summary plots using all samples (Fig. 5 and Fig. S1)


####Libraries####
library("dada2")
library("phyloseq")
library("tidyverse")
library("reshape2")
library("stringr")
library("data.table")
library("ape")
library("viridis")
library("ShortRead")
library("Biostrings")
library("seqinr")
library("data.table")
library("readr")
library("qualpalr")
library("svglite")
library("devtools")
library("microbiome")

# *** Note: Eelgrass "control" and "live" treatments are synonymous with "natural" (all refer to natural eelgrass)
# *** Note: Eelgrass "asu" (artificial seagrass unit) treatment is the same as artificial eelgrass

# setting working directory
setwd("C:/Users/Emily/Documents/seagrass_transplant_experiment/data_analysis")

# reading in phyloseq object
taxonomy_plot_obj <- readRDS("Adamczyk_taxa_summary_plots_all_samples.RDS")

# selecting rank 6 only (genus level)
taxonomy_plot_obj <- taxonomy_plot_obj  %>%
  tax_glom(taxrank = "Rank6")


####################################################################
####################################################################
####################################################################

# reading in colour table
colour_table <- read.csv("taxa_summary_color_chart.csv")
colour_table <- colour_table %>% rename("ï..Taxa" = "Taxa")

####################################################################
####################################################################
####################################################################

#########################
##### NATURAL DAY 1 #####
####  RANK 6 LEVEL  #####
####    FIGURE 5     ####
#########################

# filtering for natural day 1 eelgrass
control_only <- subset_samples(taxonomy_plot_obj, treatment == "control")
initial_control_only <- subset_samples(control_only, category == "initial_eelgrass_swab")

#creating otu table from phyloseq object
abund_table_i_c <- otu_table(initial_control_only)

# getting metadata from phyolseq object
meta_table_i_c <- sample_data(initial_control_only)

# apply proportional normalisation (relative abundance)
i_c <- abund_table_i_c/rowSums(abund_table_i_c)

# sort from most to least abundant
i_c <- i_c[,order(colSums(i_c), decreasing = TRUE)]
rowSums(i_c)

#Extract list of top N Taxa 
N<-15# keep this in mind if i need to add more taxa later depending on functional redundancy/functional traits
taxa_list<-colnames(i_c)[1:N]

#remove "__Unknown__" and add it to others
taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
N<-length(taxa_list)

# taxa_list are 15 most abundant ASVs

#Generate a new table with everything added to Others
i_c <- t(i_c)
colSums(i_c)
new_i_c<-data.frame(i_c[1:N,]) #get most abundant taxa
new_i_c <- as.data.frame(t(new_i_c)) #transpose and cast as data frame
new_i_c$Others <- colSums(i_c[(N+1):length(row.names(i_c)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset

#fix column labels
#pull taxa table subset by 21 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
alltaxa_ic <- as.data.frame(unclass(tax_table(initial_control_only)))
taxa_N <- alltaxa_ic[which(row.names(alltaxa_ic) %in% taxa_list),]

colnames(new_i_c) #taxa in order of abundance

taxa_N # taxa_N are the names for top 15 abundant ASVs

# this is for RANK 6 (without ASV numbers)
colnames(new_i_c) <- c("Saprospiraceae; uncultured", #10
                       "Methylophilaceae; Methylotenera", #4
                       "Rhodobacteraceae; Sulfitobacter", #8
                       "Thiohalorhabdaceae; Granulosicoccus", #26 
                       "Methylophilaceae; Methylophilus", #7
                       "Flavobacteriaceae; Winogradskyella", #58  
                       "Rhodobacteraceae; Pseudooctadecabacter", #43
                       "Rhodobacteraceae; Loktanella",  #14
                       "Flavobacteriaceae; Polaribacter",#9
                       "Pirellulaceae; Blastopirellula",   #135
                       "Flavobacteriaceae; Maribacter",  #35
                       "Thiotrichaceae; Cocleimonas", #36
                       "Microtrichaceae; uncultured", #90
                       "Rickettsiaceae; Rickettsia", #152
                       "Saprospiraceae; Lewinella", #102  <-  THIS IS NEW??>
                       "Others")


## save taxa names as df
taxa_names_ic = as.data.frame(colnames(new_i_c))
dput(as.character(taxa_names_ic$Taxa)) 
colnames(taxa_names_ic)[1] <- "Taxa"

#Add metadata for plotting to new_x
# THIS IS FOR ORIGIN SITE ONLY
new_i_c$origin_site <- sample_data(initial_control_only)$origin_site[which(row.names(sample_data(initial_control_only)) %in% row.names(new_i_c))]
new_i_c$sample_id <- row.names(new_i_c)
new_ic_melt <- reshape2::melt(new_i_c, 
                              id.vars = c("sample_id", "origin_site"), 
                              variable = "Taxa", na.rm=TRUE)


#joining colour table with initial control table
taxaxf = as.data.frame(taxa_names_ic$Taxa)
colnames(taxaxf)[1] <- "Taxa"
colxf2018 <- inner_join(taxaxf, colour_table)
colxf2018 <- unique(colxf2018)


# Assign taxa names to colors
x.colours <- colxf2018$colour
names(x.colours) <- colxf2018$Taxa
head(x.colours)
scales::show_col(x.colours)


  ### Taxa plot (Natural Day 1)
ggplot(new_ic_melt, aes(sample_id, value, fill=Taxa)) +
    geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid") + #height of the column equal the value, stacked
    facet_grid(. ~ origin_site, drop=TRUE,scale="free",space="free_x") +
    scale_fill_manual(values=x.colours) +
    theme_bw() + 
    ylab("Relative abundance") +
    scale_y_continuous(expand = c(0,0)) +
    theme(strip.background = element_rect(fill="white"),
          strip.text.x = element_text(size = 16))+
    theme(panel.spacing = unit(0.1, "lines")) +
    theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.y=element_text(size = 14), #change font size of numbers
                  axis.title.y=element_text(size = 18), #change font size of y title
                  panel.grid.major.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  legend.text=element_text(size=13),
                  legend.title = element_text(size = 15)) +
              guides(fill = guide_legend(ncol = 1))


####################################################################
####################################################################
####################################################################
  

#########################
##### NATURAL DAY 5 #####
####  RANK 6 LEVEL  #####
####    FIGURE 5     ####
#########################

# selecting Day 5 natural seagrass
control_only <- subset_samples(taxonomy_plot_obj, treatment == "control")
final_control_only <- subset_samples(control_only, category == "final_eelgrass_swab")

  #creating otu table from phyloseq object
  abund_table_f_c <- otu_table(final_control_only)
  
  # getting metadata from phyolseq object
  meta_table_f_c <- sample_data(final_control_only)
  
  # apply proportional normalisation (relative abundance)
  f_c <- abund_table_f_c/rowSums(abund_table_f_c)
  
  # sort from most to least abundant
  f_c <- f_c[,order(colSums(f_c), decreasing = TRUE)]
  rowSums(f_c)
  
  #Extract list of top N Taxa 
  N<-15 # keep this in mind if i need to add more taxa later depending on functional redundancy/functional traits
  taxa_list<-colnames(f_c)[1:N]
  #remove "__Unknown__" and add it to others
  taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
  N<-length(taxa_list)
  # taxa_list are 20 most abundant ASVs
  
  #Generate a new table with everything added to Others
  f_c <- t(f_c)
  colSums(f_c)
  new_f_c<-data.frame(f_c[1:N,]) #get most abundant taxa
  new_f_c <- as.data.frame(t(new_f_c)) #transpose and cast as data frame
  new_f_c$Others <- colSums(f_c[(N+1):length(row.names(f_c)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset
  
  
  #fix column labels
  #pull taxa table subset by 21 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
  alltaxa_fc <- as.data.frame(unclass(tax_table(final_control_only)))
  taxa_N <- alltaxa_fc[which(row.names(alltaxa_fc) %in% taxa_list),]
  colnames(new_f_c)
  taxa_N # taxa_N are the names for top 15 abundant ASVs
  # FINAL Natural
  
  colnames(new_f_c) <- c("Saprospiraceae; uncultured",	#10				
                         "Rhodobacteraceae; Sulfitobacter", #8					
                         "Methylophilaceae; Methylotenera", #4					
                         "Methylophilaceae; Methylophilus",	#7				
                         "Thiohalorhabdaceae; Granulosicoccus", #26 					
                         "Rhodobacteraceae; Pacificibacter", #15					
                         "Rhodobacteraceae; Loktanella",	 #14				
                         "Thiotrichaceae; Cocleimonas",	#36	
                         "Flavobacteriaceae; Maribacter", #35
                         "Thiotrichaceae; Leucothrix",	#28
                         "Verrucomicrobiales; DEV007; uncultured", #118	
                         "Flavobacteriaceae; Winogradskyella",  #58	
                         "Pirellulaceae; Blastopirellula",	#135
                         "Sphingobacteriales; KD3-93; uncultured", #70
                         "Methylophagaceae; uncultured", #123
                         "Others")	
  
  ## save taxa names as df
  taxa_names_fc = as.data.frame(colnames(new_f_c))
  colnames(taxa_names_fc)[1] <- "Taxa"
  
  #Add metadata for plotting to new_x
  # THIS IS FOR DESTINATION SITE
  new_f_c$destination <- sample_data(final_control_only)$destination[which(row.names(sample_data(final_control_only)) %in% row.names(new_f_c))]
  new_f_c$sample_id <- row.names(new_f_c)
  new_fc_melt <- reshape2::melt(new_f_c, id.vars = c("sample_id", "destination"), variable = "Taxa", na.rm=TRUE)
  
    #joining colour table with initial control table
  taxaxf = as.data.frame(taxa_names_fc$Taxa)
  colnames(taxaxf)[1] <- "Taxa"
  colxf2018 <- inner_join(taxaxf, colour_table)
  colxf2018 <- unique(colxf2018)
  
  # Assign taxa names to colors
  x.colours <- colxf2018$colour
  names(x.colours) <- colxf2018$Taxa
  head(x.colours)
  scales::show_col(x.colours)
  
  
### Taxa plot
ggplot(new_fc_melt, aes(sample_id, value, fill=Taxa)) +
    geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid") + #height of the column equal the value, stacked
    facet_grid(. ~ destination, drop=TRUE,scale="free",space="free_x") + 
    scale_fill_manual(values=x.colours) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() + 
    ylab("Relative abundance") +
    xlab("Transplant Location") +
    theme(strip.background = element_rect(fill="white"),
          strip.text.x = element_text(size = 16))+
    theme(panel.spacing = unit(0.1, "lines")) +
    theme (axis.ticks.x=element_blank(),
           axis.text.x = element_blank(),
           axis.text.y=element_text(size = 14), #change font size of numbers
           axis.title.y=element_text(size = 18), #change font size of y title
           axis.title.x=element_blank(),
           panel.grid.major.x=element_blank(),
           panel.grid.minor.x=element_blank(),
           legend.text=element_text(size=13),
           legend.title = element_text(size = 15)) +
    guides(fill = guide_legend(ncol = 1))


####################################################################
####################################################################
####################################################################


#########################
### STERILIZED DAY 5 ####
####  RANK 6 LEVEL  #####
####    FIGURE 5     ####
#########################

# filtering for sterilized eelgrass only
  sterile_only <- subset_samples(taxonomy_plot_obj, treatment == "sterile")
  final_sterile_only <- subset_samples(sterile_only, category == "final_eelgrass_swab")
  
  #creating otu table from phyloseq object
  abund_table_f_s <- otu_table(final_sterile_only)
  
  # getting metadata from phyolseq object
  meta_table_f_s <- sample_data(final_sterile_only)
  
  # apply proportional normalisation (relative abundance)
  f_s <- abund_table_f_s/rowSums(abund_table_f_s)
  
  # sort from most to least abundant
  f_s <- f_s[,order(colSums(f_s), decreasing = TRUE)]
  rowSums(f_s)
  
  #Extract list of top N Taxa 
  N<-15 # keep this in mind if i need to add more taxa later depending on functional redundancy/functional traits
  taxa_list<-colnames(f_s)[1:N]
  #remove "__Unknown__" and add it to others
  taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
  N<-length(taxa_list)
  # taxa_list are 15 most abundant ASVs
  
  #Generate a new table with everything added to Others
  f_s <- t(f_s)
  colSums(f_s)
  new_f_s<-data.frame(f_s[1:N,]) #get most abundant taxa
  new_f_s <- as.data.frame(t(new_f_s)) #transpose and cast as data frame
  new_f_s$Others <- colSums(f_s[(N+1):length(row.names(f_s)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset
  
  
  #fix column labels
  #pull taxa table subset by 21 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
  alltaxa_fs <- as.data.frame(unclass(tax_table(final_sterile_only)))
  taxa_N <- alltaxa_fs[which(row.names(alltaxa_fs) %in% taxa_list),]
  colnames(new_f_s)
  # taxa_N are the names for top 15 abundant ASVs
  
  taxa_N ###run taxa_N to see which taxa are the ones described by the number codes in colnames(new_x)
  
  # RANK 6 - FINAL STERILE
  colnames(new_f_s) <- c("Rhodobacteraceae; Sulfitobacter", #8
                         "Rhodobacteraceae; Pacificibacter", #15
                         "Saprospiraceae; uncultured", #10
                         "Rhodobacteraceae; Loktanella", #14
                         "Alteromonadaceae; Alteromonas",#16
                         "Methylophilaceae; Methylophilus", #7
                         "Thiotrichaceae; Leucothrix", #28
                         "Colwelliaceae; Colwellia", #50
                         "Methylophilaceae; Methylotenera", #4
                         "Thiotrichaceae; Cocleimonas", #36
                         "Flavobacteriaceae; Polaribacter", #9
                         "Methylophagaceae; uncultured", #123
                         "Flavobacteriaceae; Maribacter", #35
                         "Sphingobacteriales; KD3-93; uncultured", #70
                         "Ectothiorhodospiraceae; uncultured", #3
                         "Others")
  
  ## save taxa names as df
  taxa_names_fs = as.data.frame(colnames(new_f_s))
  colnames(taxa_names_fs)[1] <- "Taxa"
  
  #Add metadata for plotting to new_x
  # THIS IS FOR DESTINATION ONLY
  new_f_s$destination <- sample_data(final_sterile_only)$destination[which(row.names(sample_data(final_sterile_only)) %in% row.names(new_f_s))]
  new_f_s$sample_id<- row.names(new_f_s)
  new_fs_melt <- reshape2::melt(new_f_s, id.vars = c("sample_id", "destination"), variable = "Taxa", na.rm=TRUE)
  
  
  #joining colour table with initial control table
  taxaxf = as.data.frame(taxa_names_fs$Taxa)
  colnames(taxaxf)[1] <- "Taxa"
  colxf2018 <- inner_join(taxaxf, colour_table)
  colxf2018 <- unique(colxf2018)
  
  # Assign taxa names to colors
  x.colours <- colxf2018$colour
  names(x.colours) <- colxf2018$Taxa
  head(x.colours)
  scales::show_col(x.colours)
  
  ### Taxa plot
ggplot(new_fs_melt, aes(sample_id, value, fill=Taxa)) +
    geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid") + #height of the column equal the value, stacked
    facet_grid(. ~ destination, drop=TRUE,scale="free",space="free_x") + 
    scale_fill_manual(values=x.colours) +
    theme_bw() + 
    ylab("Relative abundance") +
    xlab("Transplant Location") +
    scale_y_continuous(expand = c(0,0))+
    theme(strip.background = element_rect(fill="white"),
          strip.text.x = element_text(size = 14))+
    theme(panel.spacing = unit(0.1, "lines")) +
    theme (axis.ticks.x=element_blank(),
           axis.text.x = element_blank(),
           axis.text.y=element_text(size = 14), #change font size of numbers
           axis.title.y=element_text(size = 18), #change font size of y title
           panel.grid.major.x=element_blank(),
           panel.grid.minor.x=element_blank(),
           legend.text=element_text(size=13),
           legend.title = element_text(size = 15),
           plot.title = element_text(hjust = 0.5, size = 20, face="bold")) +#center plot title and set font size
    guides(fill = guide_legend(ncol = 1))


####################################################################
####################################################################
####################################################################


#########################
### ARTIFICIAL DAY 5 ####
####  RANK 6 LEVEL  #####
####    FIGURE 5     ####
#########################

# selecting artificial eelgrass only
asu_only <- subset_samples(taxonomy_plot_obj, treatment == "asu")
final_asu_only <- subset_samples(asu_only, category == "final_eelgrass_swab")

#creating otu table from phyloseq object
abund_table_f_asu <- otu_table(final_asu_only)

# getting metadata from phyolseq object
meta_table_f_asu <- sample_data(final_asu_only)

# apply proportional normalisation (relative abundance)
f_asu <- abund_table_f_asu/rowSums(abund_table_f_asu)

# sort from most to least abundant
f_asu <- f_asu[,order(colSums(f_asu), decreasing = TRUE)]
rowSums(f_asu)

#Extract list of top N Taxa 
N<-15 # keep this in mind if i need to add more taxa later depending on functional redundancy/functional traits
taxa_list<-colnames(f_asu)[1:N]
#remove "__Unknown__" and add it to others
taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
N<-length(taxa_list)
# taxa_list are 15 most abundant ASVs

#Generate a new table with everything added to Others
f_asu <- t(f_asu)
colSums(f_asu)
new_f_asu<-data.frame(f_asu[1:N,]) #get most abundant taxa
new_f_asu <- as.data.frame(t(new_f_asu)) #transpose and cast as data frame
new_f_asu$Others <- colSums(f_asu[(N+1):length(row.names(f_asu)),]) #making final column which is a sum of relative abundances for the rest of the taxa in the dataset

#fix column labels
#pull taxa table subset by 21 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
alltaxa_asu <- as.data.frame(unclass(tax_table(final_asu_only)))
taxa_N <- alltaxa_asu[which(row.names(alltaxa_asu) %in% taxa_list),]
colnames(new_f_asu)

taxa_N # taxa_N are the names for top 15 abundant ASVs

# RANK 6 - ASU final
colnames(new_f_asu) <- c("Saprospiraceae; uncultured", #10
                         "Rhodobacteraceae; Sulfitobacter",  #8
                         "Ectothiorhodospiraceae; uncultured", #3
                         "Rhodobacteraceae; Loktanella", #14
                         "Thiotrichaceae; Leucothrix", #28
                         "Colwelliaceae; Colwellia",  #50
                         "Alteromonadaceae; Alteromonas", #16
                         "Rhodobacteraceae; Pacificibacter", #15
                         "Pirellulaceae; Blastopirellula", #135
                         "Flavobacteriaceae; Maribacter", #35
                         "Rhodobacteraceae; Amylibacter",#83
                         "Flavobacteriaceae; Winogradskyella", #58
                         "Oleiphilaceae; Oleiphilus", #79
                         "Thiohalorhabdaceae; Granulosicoccus", #26
                         "Thiotrichaceae; Cocleimonas", #36
                         "Others")


## save taxa names as df
taxa_names_f_asu = as.data.frame(colnames(new_f_asu))
colnames(taxa_names_f_asu)[1] <- "Taxa"

#Add metadata for plotting to new_x
# THIS IS FOR DESTINATION ONLY
new_f_asu$destination <- sample_data(final_asu_only)$destination[which(row.names(sample_data(final_asu_only)) %in% row.names(new_f_asu))]
new_f_asu$sample_id <- row.names(new_f_asu)
new_asu_melt <- reshape2::melt(new_f_asu, id.vars = c("sample_id", "destination"), variable = "Taxa", na.rm=TRUE)

#joining colour table with initial control table
taxaxf = as.data.frame(taxa_names_f_asu$Taxa)
colnames(taxaxf)[1] <- "Taxa"
colxf2018 <- inner_join(taxaxf, colour_table)
colxf2018 <- unique(colxf2018)

# Assign taxa names to colors
x.colours <- colxf2018$colour
names(x.colours) <- colxf2018$Taxa
head(x.colours)
scales::show_col(x.colours)


### Taxa plot
ggplot(new_asu_melt, aes(sample_id, value, fill=Taxa)) +
  geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid") + #height of the column equal the value, stacked
  facet_grid(. ~ destination, drop=TRUE,scale="free",space="free_x") + 
  scale_fill_manual(values=x.colours) +
  scale_y_continuous(expand = c(0,0))+
  theme_bw() + 
  ylab("Relative abundance") +
  theme(strip.background = element_rect(fill="white"),
        strip.text.x = element_text(size = 12))+
  theme(panel.spacing = unit(0.1, "lines")) +
  theme (axis.ticks.x=element_blank(),
         #axis.text.x = element_text(size = 12, angle = 90),
         axis.text.x = element_blank(),
         axis.text.y=element_text(size = 14), #change font size of numbers
         axis.title.y=element_text(size = 18), #change font size of y title
         axis.title.x=element_blank(), #change font size of x title
         panel.grid.major.x=element_blank(),
         panel.grid.minor.x=element_blank(),
         legend.text=element_text(size=13),
         legend.title = element_text(size = 15),
         plot.title = element_text(hjust = 0.5, size = 20, face="bold")) +#center plot title and set font size
  guides(fill = guide_legend(ncol = 1))


####################################################################
####################################################################
####################################################################


#########################
### STERILIZED DAY 1 ####
####  RANK 6 LEVEL  #####
####    FIGURE S1    ####
#########################

# sterilized Day 1 eelgrass only
  sterile_only <- subset_samples(taxonomy_plot_obj, treatment == "sterile")
  initial_sterile_only <- subset_samples(sterile_only, category == "initial_eelgrass_swab")

    #creating otu table from phyloseq object
  abund_table_i_s <- otu_table(initial_sterile_only)
  
  # getting metadata from phyolseq object
  meta_table_i_s <- sample_data(initial_sterile_only)
  
  # apply proportional normalisation (relative abundance)
  i_s <- abund_table_i_s/rowSums(abund_table_i_s)
  
  # sort from most to least abundant
  i_s <- i_s[,order(colSums(i_s), decreasing = TRUE)]
  rowSums(i_s)
  
  #Extract list of top N Taxa 
  N<-15 # keep this in mind if i need to add more taxa later depending on functional redundancy/functional traits
  taxa_list<-colnames(i_s)[1:N]
  #remove "__Unknown__" and add it to others
  taxa_list<-taxa_list[!grepl("Unassigned",taxa_list)]
  N<-length(taxa_list)
  # taxa_list are 15 most abundant ASVs
  
  #Generate a new table with everything added to Others
  i_s <- t(i_s)
  colSums(i_s)
  new_i_s<-data.frame(i_s[1:N,]) #get most abundant taxa
  new_i_s <- as.data.frame(t(new_i_s)) #transpose and cast as data frame
  new_i_s$Others <- colSums(i_s[(N+1):length(row.names(i_s)),]) #making initial column which is a sum of relative abundances for the rest of the taxa in the dataset
  
  #fix column labels
  #pull taxa table subset by 21 most prevalent taxa #visually inspect, then make a vector of N=22 for the new more informative column labels
  alltaxa_fs <- as.data.frame(unclass(tax_table(initial_sterile_only)))
  taxa_N <- alltaxa_fs[which(row.names(alltaxa_fs) %in% taxa_list),]
  colnames(new_i_s)
  # taxa_N are the names for top 15 abundant ASVs
  
  taxa_N ###run taxa_N to see which taxa are the ones described by the number codes in colnames(new_x)
  
  ## MAKE SURE THIS MATCHES UP WITH TAXA N
  colnames(new_i_s) <- c("Alteromonadaceae; Alteromonas", #16
                         "Pseudomonadaceae; Pseudomonas",  #52
                         "Marinobacteraceae; Marinobacter", #47
                         "Verrucomicrobiales; DEV007; uncultured", #118
                         "Kordiimonadaceae; Kordiimonas", #110
                         "Halomonadaceae; Halomonas", #111
                         "Marinomonadaceae; Marinomonas", #67
                         "Pseudoalteromonadaceae; Pseudoalteromonas", #39
                         "Microtrichaceae; uncultured", #90
                         "Rhodobacteraceae; Sulfitobacter", #8
                         "Methylophagaceae; uncultured", #123
                         "Ectothiorhodospiraceae; uncultured", #3
                         "Saprospiraceae; uncultured",#10
                         "Thalassospiraceae; Thalassospira", #280
                         "Rhodobacteraceae; Loktanella", #14
                         "Others")
  

  ## save taxa names as df
  taxa_names_is = as.data.frame(colnames(new_i_s))
  colnames(taxa_names_is)[1] <- "Taxa"
  
#Add metadata for plotting to new_x
  # THIS IS FOR ORIGIN ONLY
  new_i_s$origin_site <- sample_data(initial_sterile_only)$origin_site[which(row.names(sample_data(initial_sterile_only)) %in% row.names(new_i_s))]
  new_i_s$sample_id<- row.names(new_i_s)
  new_is_melt <- melt(new_i_s, id.vars = c("sample_id", "origin_site"), variable = "Taxa", na.rm=TRUE)
  
  #joining colour table with initial control table
  taxaxf = as.data.frame(taxa_names_is$Taxa)
  colnames(taxaxf)[1] <- "Taxa"
  colxf2018 <- inner_join(taxaxf, colour_table)
  colxf2018 <- unique(colxf2018)
  
  # Assign taxa names to colors
  x.colours <- colxf2018$colour
  names(x.colours) <- colxf2018$Taxa
  head(x.colours)
  scales::show_col(x.colours)

  
  
  ### Taxa plot
  ggplot(new_is_melt, aes(sample_id, value, fill=Taxa)) +
    geom_bar(stat="identity", width=1, color="black", position = "stack", linetype = "solid") + #height of the column equal the value, stacked
    facet_grid(. ~ origin_site, drop=TRUE,scale="free",space="free_x") + 
    #ggtitle("Sterile/killed eelgrass (initial swab)") + 
    scale_fill_manual(values=x.colours) +
    scale_y_continuous(expand = c(0,0))+
    theme_bw() + 
    ylab("Relative abundance") +
    xlab("Transplant Location") +
    theme(strip.background = element_rect(fill="white"),
          strip.text.x = element_text(size = 16))+
    theme(panel.spacing = unit(0.1, "lines")) +
    theme (axis.ticks.x=element_blank(),
           #axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
           axis.text.y=element_text(size = 14), #change font size of numbers
           axis.title.y=element_text(size = 18), #change font size of y title
           #axis.title.x=element_text(size = 18), #change font size of x title
           axis.text.x = element_blank(),
           axis.title.x = element_blank(),
           panel.grid.major.x=element_blank(),
           panel.grid.minor.x=element_blank(),
           legend.text=element_text(size=13),
           legend.title = element_text(size = 15),
           plot.title = element_text(hjust = 0.5, size = 20, face="bold")) +#center plot title and set font size
    guides(fill = guide_legend(ncol = 1))

  

  