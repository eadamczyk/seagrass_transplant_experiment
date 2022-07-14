# Project: Seagrass (Zostera marina) transplant experiment reveals core microbiome and resistance to environmental change 
# Code by: Emily Adamczyk
# Purpose: Processing seagrass transplant experiment sequences and rarefying 


####Libraries####
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2")
library(tidyverse)
library(dada2)
library(phyloseq)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(ape)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)


################################
#       DADA2 PROCESSING       #
# ORIGINAL CODE BY EVAN MORIEN #
################################


####Environment Setup####
theme_set(theme_bw())
setwd("C:/Users/Emily/Documents/seagrass_transplant_experiment/data_analysis/")

####File Path Setup####
#this is so dada2 can quickly iterate through all the R1 and R2 files in your read set
path <- "C:/Users/Emily/Documents/seagrass_transplant_experiment/sequence_data/fastq_files" # directory containing the fastq files
list.files(path)
fnFs <- sort(list.files(path, pattern="r1.fq.gz", full.names = TRUE)) #change the pattern to match all your R1 files
fnRs <- sort(list.files(path, pattern="r2.fq.gz", full.names = TRUE)) #same for this pattern, for R2 files
sample.names <- sapply(strsplit(basename(fnFs), "_r1.fq.gz"), `[`, 1) #change the delimiter in quotes and the number at the end of this command to decide how to split up the file name, and which element to extract for a unique sample name
sample.names <- gsub("^_", "", sample.names) #remove leading underscore #aesthetic change only

####fastq Quality Plots####
plotQualityProfile(fnFs[100:110]) #this plots the quality profiles for each sample, if you have a lot of samples, it's best to look at just a few of them, the plots take a minute or two to generate even only showing 10-20 samples.
plotQualityProfile(fnRs[287:293])

####Primer Removal####
####identify primers####
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME to your reverse primer sequence
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[112]]), #add the index of the sample you'd like to use for this test (your first sample may be a blank/control and not have many sequences in it, be mindful of this)
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[112]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[112]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[112]]))
which(sample.names == "unknown")

####OPTIONAL!!!!####
#REV <- REV.orients[["RevComp"]] #IMPORTANT!!! change orientation ONLY IF NECESSARY. see the dada2 ITS_workflow guide section "Identify Primers" for details. it is linked at the top of this guide.

#### primer removal ####
cutadapt <- "/anaconda2/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

#Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#sanity check, should report zero for all orientations and read sets
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[112]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[112]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[112]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[112]]))

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "r1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "r2", full.names = TRUE))

####filter and trim reads####
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

####trim & filter####
#filter and trim command. dada2 can canonically handle lots of errors, I am typically permissive in the maxEE parameter set here, in order to retain the maximum number of reads possible. error correction steps built into the dada2 pipeline have no trouble handling data with this many expected errors.
#it is best, after primer removal, to not truncate with 18s data, or with data from any region in which the length is broadly variable. you may exclude organisms that have a shorter insert than the truncation length (definitely possible, good example is giardia). defining a minimum sequence length is best.
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, truncLen=c(0,0), minLen = c(150,120),
                     maxN=c(0,0), maxEE=c(8,10), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE,
                     compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
View(retained)

####learn error rates####
#the next three sections (learn error rates, dereplication, sample inference) are the core of dada2's sequence processing pipeline. read the dada2 paper and their online documentation (linked at top of this guide) for more information on how these steps work
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE) #assess this graph. it shows the error rates observed in your dataset. strange or unexpected shapes in the plot should be considered before moving on.

####dereplication####
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names #this is just to ensure that all your R objects have the same sample names in them
names(derepFs) <- sample.names
names(derepRs) <- sample.names

####sample inference####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]
dadaRs[[1]]

####OPTIONAL: remove low-sequence samples before merging####
#a "subscript out of bounds" error at the next step (merging) may indicate that you aren't merging any reads in one or more samples.
#NB, NOT getting this error doesn't necessarily mean that all of your samples end up with more than 0 merged reads, as i found out while processing a large 18s dataset. your guess is as good as mine as to why this error does or does not appear, but filtering out the samples that cause it is necessary for completion of the pipeline.
#samples_to_keep <- as.numeric(out[,"reads.out"]) > 100 #example of simple method used above after the filter and trim step. if you already did this but still got an error when merging, try the steps below
#getN <- function(x) sum(getUniques(x)) #keeping track of read retention, number of unique sequences after ASV inference
#track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
#samples_to_keep <- track[,4] > 50 #your threshold. try different ones to get the lowest one that will work. #this method accounts for dereplication/ASVs left after inference
#samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)] #record names of samples you have the option of removing

####merge paired reads####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

####construct sequence table####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#[1]  293 7515

####View Sequence Length Distribution Post-Merging####
#most useful with merged data. this plot will not show you much for forward reads only, which should have a uniform length distribution.
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab)))) #tabulate sequence length distribution
plot(x=length.histogram[,1], y=length.histogram[,2]) #view length distribution plot

####remove low-count singleton ASVs####
#create phyloseq otu_table
otus <- otu_table(t(seqtab), taxa_are_rows = TRUE)

#some metrics from the sequence table
otu_pres_abs <- otus
otu_pres_abs[otu_pres_abs >= 1] <- 1 #creating a presence/absence table
otu_pres_abs_rowsums <- rowSums(otu_pres_abs) #counts of sample per ASV
length(otu_pres_abs_rowsums) #how many ASVs
length(which(otu_pres_abs_rowsums == 1)) #how many ASVs only present in one sample
#[1] 3942

#what are the counts of each ASV
otu_rowsums <- rowSums(otus) #raw counts per ASV
otu_singleton_rowsums <- as.data.frame(otu_rowsums[which(otu_pres_abs_rowsums == 1)]) #raw read counts in ASVs only presesnt in one sample
hist(otu_singleton_rowsums[,1], breaks=500, xlim = c(0,200), xlab="# Reads in ASV") #histogram plot of above
length(which(otu_singleton_rowsums <= 50)) #how many ASVs are there with N reads or fewer? (N=50 in example)
#[1] 3224

#IF you want to filter out rare variants (low-read-count singleton ASVs) you can use phyloseq's "transform_sample_counts" to create a relative abundance table, and then filter your ASVs by choosing a threshold of relative abundance: otus_rel_ab = transform_sample_counts(otus, function(x) x/sum(x))
dim(seqtab) # sanity check
#[1]  293 7515
dim(otus) # (this should be the same as last command, but the dimensions reversed)
otus_rel_ab <- transform_sample_counts(otus, function(x) x/sum(x)) #create relative abundance table
df <- as.data.frame(unclass(otus_rel_ab)) #convert to plain data frame
df[is.na(df)] <- 0 #if there are samples with no merged reads in them, and they passed the merge step (a possiblity, converting to a relative abundance table produes all NaNs for that sample. these need to be set to zero so we can do the calculations in the next steps.)
otus_rel_ab.rowsums <- rowSums(df) #compute row sums (sum of relative abundances per ASV. for those only present in one sample, this is a value we can use to filter them for relative abundance on a per-sample basis)
a <- which(as.data.frame(otu_pres_abs_rowsums) == 1) #which ASVs are only present in one sample
b <- which(otus_rel_ab.rowsums <= 0.001) #here is where you set your relative abundance threshold #which ASVs pass our filter for relative abundance
length(intersect(a,b)) #how many of our singleton ASVs fail on this filter
#[1] 3513
rows_to_remove <- intersect(a,b) #A also in B (we remove singleton ASVs that have a lower relative abundance value than our threshold)
otus_filt <- otus[-rows_to_remove,] #filter OTU table we created earlier
dim(otus_filt) #how many ASVs did you retain?
#[1] 4002  293
seqtab.nosingletons <- t(as.matrix(unclass(otus_filt))) #convert filtered OTU table back to a sequence table matrix to continue with dada2 pipeline

####remove chimeras####
#here we remove "bimeras" or chimeras with two sources. look at "method" to decide which type of pooling you'd like to use when judging each sequence as chimeric or non-chimeric
seqtab.nosingletons.nochim <- removeBimeraDenovo(seqtab.nosingletons, method="pooled", multithread=TRUE, verbose=TRUE) #this step can take a few minutes to a few hours, depending on the size of your dataset
dim(seqtab.nosingletons.nochim)
#293 3715
sum(seqtab.nosingletons.nochim)/sum(seqtab.nosingletons) #proportion of nonchimeras #it should be relatively high after filtering out your singletons/low-count ASVs, even if you lose a lot of ASVs, the number of reads lost should be quite low
#[1] 0.9819295

####track read retention through steps####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nosingletons), rowSums(seqtab.nosingletons.nochim))
# If processing only a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nosingletons", "nochimeras", "percent_singletons", "percent_chimeras")
rownames(track) <- sample.names

####save output from sequnce table construction steps####
write.table(data.frame("row_names"=rownames(track),track),"read_retention.16s.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nosingletons.nochim),seqtab.nosingletons.nochim),"sequence_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")

#if you must save your sequence table and load it back in before doing taxonomy assignments, here is how to reformat the object so that dada2 will accept it again
seqtab.nosingletons.nochim <- fread("sequence_table.16s.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with the row names in it
seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix

####assign taxonomy####
#note, this takes ages if you have a large dataset. saving the sequences as a fasta file (with writeFasta) and using QIIME's taxonomy assignment command will save you time, and is only slightly less accurate than the dada2 package's taxonomy assignment function.
taxa <- assignTaxonomy(seqtab.nosingletons.nochim, "~/Desktop/lab_member_files/taxonomy_databases/silva_for_dada2/v132_for_parfreylab/16s/silva_132.16s.99_rep_set.dada2.fa.gz", multithread=TRUE)

#NA taxa are hard to separate later if they have no label. apply "Unassigned" label here now.
unique(taxa[,1]) #possible labels here: eukaryotic, archaeal, bacterial, and "NA" taxa. 
NAs <- is.na(taxa[,1]) #test for NA
NAs <- which(NAs == TRUE) #get indices of NA values
taxa[NAs,1] <- "Unassigned" #apply new label to identified indices
colnames(taxa) <- c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6", "Rank7", "Accession")

#identifying eukaryotic sequences using the dada2 dev's version of SILVA 132 (faster taxonomy assignment, remove those ASVs that have a phylumn assignment in the eukaryotic domain)
# taxa_dada2 <- assignTaxonomy(seqtab.nosingletons.nochim, "~/Desktop/lab_member_files/taxonomy_databases/silva_for_dada2/v132_by_dada2_devs/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# 
# #removing high-confidence eukaryotic taxa (i won't remove all sequences that only get an assignment at rank1 for 18s, some of them are real, i know from using blast. a better database would label them.)
# a <- which(taxa_dada2[,1] == "Eukaryota") #which have a eukaryotic assignment
# b <- which((is.na(taxa_dada2[,2]) == F)) #which taxa have an assignment at phylum level (we can be more confident that these are really in the kingdoms they say they are in)
# c <- intersect(a,b)
# 
# d <- which((is.na(taxa[,2]) == T)) #which taxa do not have a phylum assignment from silva 16s
# eukaryotic_taxa <- Reduce(intersect, list(c,d)) #which taxa share these two conditions, i think we can safely remove them
# 
# taxa <- taxa[-eukaryotic_taxa,] #we didn't actually end up removing any taxa with this step

####saving taxonomy data####
write.table(data.frame("row_names"=rownames(taxa),taxa),"taxonomy_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")



#reading in taxonomy table
taxa <- fread("taxonomy_table.16s.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(taxa) <- taxa[,1] #set row names
taxa <- taxa[,-1] #remove column with the row names in it
taxa <- as.matrix(taxa) #cast the object as a matrix


#if you must save your sequence table and load it back in before doing taxonomy assignments, here is how to reformat the object so that dada2 will accept it again
seqtab.nosingletons.nochim <- fread("sequence_table.16s.merged.txt", sep="\t", header=T, colClasses = c("row_names"="character"), data.table=FALSE)
row.names(seqtab.nosingletons.nochim) <- seqtab.nosingletons.nochim[,1] #set row names
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,-1] #remove column with the row names in it
seqtab.nosingletons.nochim <- as.matrix(seqtab.nosingletons.nochim) #cast the object as a matrix



## hand off to PhyloSeq ##

#load sample data
rawmetadata <- read_delim(file = file.path("C:/Users/Emily/Documents/Hakai Eelgrass Microbiome Transplant 2017/sequence_data/metadata_seagrass_emily_20190425.with_barcodes.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not "escape" quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries

# joining "origin_site" and "destination" into one new column called "origin_and_destination"
rawmetadata <- rawmetadata %>%  
  dplyr::mutate("origin_and_destination" = paste(origin_site, destination, sep = "_")) #paste is really only good for character values

# renaming #SampleID to sample_id
names(rawmetadata)[names(rawmetadata) == "#SampleID"] <- "sample_id"

#set row names as the sample IDs for the metadata
row.names(rawmetadata) <- rawmetadata$sample_id

# setting rawmetadata as data frame
rawmetadata <- as.data.frame(rawmetadata)

#record samples absent in either metadata or OTU table
notinmeta <- setdiff(row.names(seqtab.nosingletons.nochim), row.names(rawmetadata)) #only the "unknown" sample (reads not assigned to a barcode by deML - we could have removed these before running dada2, but it doesn"t matter much) are not in the metadata. this is fine.
notinraw <- setdiff(row.names(rawmetadata), row.names(seqtab.nosingletons.nochim))

#create phyloseq object from "seqtab.nosingletons.nochim", "rawmetadata", and "taxa"
ps.dada2_join <- phyloseq(otu_table(seqtab.nosingletons.nochim, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(taxa))

# now replace the long ASV names (the actual sequences) with human-readable names, and save the new names and sequences as a .fasta file in your project working directory
my_otu_table <- t(as.data.frame(unclass(otu_table(ps.dada2_join)))) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ntaxa(ps.dada2_join)), sep="") #create new names
#write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "seagrass_transplant_expt.16s_ASV_sequences.fasta") #save sequences with new names in fasta format
taxa_names(ps.dada2_join) <- ASV.num #rename your sequences in the phyloseq object

# removing Rhea's samples (12 of them)
remove <- c("r300_rhea", "r302_rhea", "r303_rhea", "r304_rhea", "r305_rhea", "r307_rhea", "r308_rhea", "r312_rhea", "r313_rhea", "r314_rhea", "r315_rhea", "r319_rhea")
keep <- sample_names(ps.dada2_join)[-which(sample_names(ps.dada2_join) %in% remove)]
ps.dada2_join <- prune_samples(keep, ps.dada2_join)

#saveRDS(ps.dada2_join, "seagrass_transplant_expt.full_dataset.phyloseq_format.RDS") # saving full dataset

# removing samples that were run at 32 PCR cycles (76 samples)

##### note: skip this step if retaining all samples
remove2 <- c("pb10_pruth_bay_sterile_i_p1_h8", 
             "pb11_wolf_sterile_i_p2_h1", 
             "pb13_sandspit_asu_f_p1_g4", 
             "pb13_sandspit_control_f_p1_a7", 
             "pb13_sandspit_sterile_i_p2_f1", 
             "pb14_pruth_bay_asu_f_p1_d1", 
             "pb14_pruth_bay_control_f_p1_d3", 
             "pb14_pruth_bay_sterile_i_p2_a6", 
             "pb15_wolf_sterile_i_p1_a6", 
             "pb1_pruth_pocket_control_f_p1_g3", 
             "pb1_pruth_pocket_sterile_i_p2_b2", 
             "pb2_sandspit_control_i_p1_b7", 
             "pb4_wolf_control_i_p1_g6", 
             "pb5_pruth_pocket_sterile_f_p2_g6", 
             "pb5_pruth_pocket_sterile_i_p1_f8", 
             "pb6_pruth_pocket_asu_f_p1_b3", 
             "pb9_sandspit_asu_f_p2_e8", 
             "pb9_sandspit_sterile_f_p2_d8", 
             "pb9_sandspit_sterile_i_p1_a5", 
             "pp10_wolf_sterile_i_p2_b10", 
             "pp11_pruth_bay_sterile_f_p1_h3", 
             "pp11_pruth_bay_sterile_i_p2_a1", 
             "pp14_sandspit_control_i_p1_h1", 
             "pp14_sandspit_sterile_i_p1_g1", 
             "pp15_wolf_sterile_i_p1_c1", 
             "pp1_wolf_control_i_p1_a8", 
             "pp2_pruth_pocket_asu_f_p2_h8", 
             "pp3_pruth_pocket_sterile_i_p1_d2", 
             "pp4_pruth_pocket_asu_f_p1_h11", 
             "pp4_pruth_pocket_control_i_p1_c2", 
             "pp5_sandspit_control_f_p1_b8", 
             "pp6_sandspit_control_f_p2_c9", 
             "pp6_sandspit_sterile_i_p2_f3", 
             "pp8_pruth_bay_asu_f_p1_a3", 
             "pp8_pruth_bay_control_i_p2_f10", 
             "pp8_pruth_bay_sterile_i_p1_e1", 
             "pruth_bay_water_deployment_p1_b2", 
             "pruth_bay_water_recovery_p2_e4", 
             "pruth_pocket_water_deployment_p1_f1", 
             "pruth_pocket_water_deployment_p2_d3", 
             "pruth_pocket_water_recovery_p1_c11", 
             "pruth_pocket_water_recovery_p2_d4", 
             "sandspit_water_deployment_p2_a3", 
             "ss10_pruth_pocket_sterile_i_p1_a4", 
             "ss11_pruth_bay_asu_f_p1_g2", 
             "ss11_pruth_bay_control_i_p1_a10", 
             "ss11_pruth_bay_sterile_f_p1_f3", 
             "ss11_pruth_bay_sterile_i_p1_d4", 
             "ss12_wolf_control_f_p1_h5", 
             "ss12_wolf_control_i_p1_b11", 
             "ss13_wolf_sterile_i_p1_g9", 
             "ss14_sandspit_asu_f_p1_c5", 
             "ss14_sandspit_control_i_p2_e5", 
             "ss14_sandspit_sterile_i_p1_h9", 
             "ss15_sandspit_sterile_i_p2_g3", 
             "ss2_pruth_pocket_asu_f_p2_g11", 
             "ss3_wolf_sterile_i_p2_e3", 
             "ss7_pruth_bay_sterile_i_p2_h3", 
             "ss8_sandspit_control_i_p1_d11", 
             "ss9_pruth_pocket_control_i_p1_c4", 
             "ss9_pruth_pocket_sterile_i_p2_d1", 
             "wo10_sandspit_asu_f_p2_h9", 
             "wo10_sandspit_control_i_p1_b1", 
             "wo10_sandspit_sterile_i_p2_f2", 
             "wo11_pruth_bay_asu_f_p1_c7", 
             "wo11_pruth_bay_control_i_p1_a11", 
             "wo11_pruth_bay_sterile_i_p1_d8", 
             "wo13_pruth_pocket_sterile_i_p2_g4", 
             "wo15_pruth_pocket_control_i_p1_h4", 
             "wo15_pruth_pocket_sterile_i_p1_e10", 
             "wo1_wolf_sterile_i_p1_f2", 
             "wo4_sandspit_sterile_i_p1_a1", 
             "wo5_wolf_control_i_p1_h2", 
             "wo8_wolf_control_i_p2_g2", 
             "wo8_wolf_sterile_i_p2_h7", 
             "wolf_water_deployment_p1_f9",
             "plate_1.2_blank_1",
             "plate_1.2_blank_2",
             "plate_2.2_blank_1",
             "plate_2.2_blank_2"
)


keep2 <- sample_names(ps.dada2_join)[-which(sample_names(ps.dada2_join) %in% remove2)]
ps.dada2_join <- prune_samples(keep2, ps.dada2_join)


# looking at total number of ASVs
dim(tax_table(ps.dada2_join)) # there are 3715

# viewing table
#View(as.data.frame(tax_table(ps.dada2_join)))

# Set aside unassigned taxa #with dada2, there may not be any unassigned taxa as dada2"s RDP classifier usually ends up classifying everything. You may want to adjust this command to set aside ASVs that have no assignment beyond Rank1 instead of no assignment at Rank1.
ps.dada2_join.unassigned <- ps.dada2_join %>%
  subset_taxa(Rank1 == "Unassigned") # there are no Rank1 unassigned 
# Remove unassigned taxa - not necessary because I do not have any Rank1 unassigned taxa
# ps.dada2_join <- ps.dada2_join %>%
# subset_taxa(Rank1 != "Unassigned")

#View(as.data.frame(sort_samples(ps.dada2_join)))

# Remove ASV counts of 2 or less from OTU table #this is a de-noising method.
otu <- as.data.frame(otu_table(ps.dada2_join))
otu_table(ps.dada2_join)[otu <= 2] <- 0

# viewing table
#View(as.data.frame(tax_table(ps.dada2_join)))

#after denoising you can also remove ASVs that are in groups you wish to exclude (i.e. mammalia, embryophyta, etc.)
#to do this, just determine which rank the clade is captured by, and filter like so:
# looking at number of chloroplast ASVss
ps.dada2_chloroplast <- ps.dada2_join %>%
  subset_taxa(Rank4 == "Chloroplast"| is.na(Rank4))
#View(as.data.frame(tax_table(ps.dada2_chloroplast)))

# looking at number of mitochondrial ASVs
ps.dada2_mito <- ps.dada2_join %>%
  subset_taxa(Rank5 == "Mitochondria" | is.na(Rank5))
#View(as.data.frame(tax_table(ps.dada2_mito)))

# merging chloroplast and mitochondria dataframes together
ps.dada2_chloro_mito <- merge_phyloseq(ps.dada2_chloroplast, ps.dada2_mito)
#View(as.data.frame(tax_table(ps.dada2_chloro_mito)))

dim(tax_table(ps.dada2_chloro_mito)) #356 if not using `is.na`, 931 if using `is.na` for chloroplast or mito

chl_mito_otu <- as.data.frame(unclass(otu_table(ps.dada2_mito)))


# Remove mitochondrial and chloroplast ASVs
ps.dada2_join <- ps.dada2_join %>%
  subset_taxa(Rank4 != "Chloroplast" | is.na(Rank4)) %>% #you can chain as many of these subset_taxa calls as you like into this single command using the R pipe (%>%)
  subset_taxa(Rank5 != "Mitochondria" | is.na(Rank5))

#View(as.data.frame(tax_table(ps.dada2_join)))

sort(sample_sums(ps.dada2_join))
summary(sort(sample_sums(ps.dada2_join)))
dim(tax_table(ps.dada2_join)) #931 chloroplast or mito

#View(as.data.frame(sample_sums(ps.dada2_join)))

# Remove samples with less than N reads (N=500 in example. adjust per experiment.)
plot(sort(sample_sums(ps.dada2_join)), ylim=c(0, 2000)) 
# looking at sample_sums to see what the samples are
sort(sample_sums(ps.dada2_join))
#looks like there is an easy choice here at around 900 reads, (is the group of samples with ~1000-2000 reads a specific set of samples? do they share attributes?) 
# decided to remove samples with less than 900 reads
ps.dada2_join <- prune_samples(sample_sums(ps.dada2_join) >= 900, ps.dada2_join)

# removing that outlier sandspit sample
#View(as.data.frame(sample_data(ps.dada2_join)))

# removing the following 2 outliers for when 76 samples were filtered
ps.dada2_join <- subset_samples(ps.dada2_join, sample_id != "wo8_wolf_control_f_p1_c8")
ps.dada2_join <- subset_samples(ps.dada2_join, sample_id != "wo9_pruth_pocket_control_f_p1_a9")


#ps.dada2_join <- subset_samples(ps.dada2_join, sample_id != "pb9_sandspit_asu_f_p2_e8") # only removing this outlier when including all samples in the dataset
max(sample_sums(ps.dada2_join))

summary(sample_sums(ps.dada2_join))

summarize_phyloseq(ps.dada2_join)

dim(sample_sums(ps.dada2_join))



## All ASVs are assigned at Rank 1
tax.table <- as.data.frame(unclass(tax_table(ps.dada2_join)))
length(which(tax.table[,1] == "Bacteria"))
length(which(tax.table[,1] == "Archaea"))
length(which(is.na(tax.table[,2]) == TRUE))
dim(tax.table)

# Rarefying filtered data (data with large outlier removed)

#found rarefaction curve function here: https://rdrr.io/github/gauravsk/ranacapa/src/R/ggrare.R #download this script
#to use this you have to load the ggrare() function in from ggrare.R
source("~/Hakai Eelgrass Microbiome Transplant 2017/data_analysis/ggrare_function.R") # load in ggrare function
#here are the usage instructions: https://rdrr.io/github/gauravsk/ranacapa/man/ggrare.html
#and here is an example of usage in this workflow

# plotting the limits of the number of reads per sample
plot(sort(sample_sums(ps.dada2_join)), ylim = c(0,2000), xlim = c(0,120))
abline(h = 900)

#you can use the "which" command like the example below to see which samples you"ll lose for a given cutoff value
sort(sample_sums(ps.dada2_join))
which(sample_sums(ps.dada2_join) < 900) #rarefying to 900

# making a plot using ggrare to take a look at rarefied data
p <- ggrare(ps.dada2_join, step = 900, color = "treatment", se = FALSE)


#### change metadata to factors where necessary####
#changing destination and origin site to factors
sample_data(ps.dada2_join)$destination <- factor(sample_data(ps.dada2_join)$destination, levels = c("pruth_bay", "pruth_pocket", "sandspit", "wolf", "na"))
sample_data(ps.dada2_join)$origin_site <- factor(sample_data(ps.dada2_join)$origin_site)


sort(sample_sums(ps.dada2_join))

#####################
#### rarefy data ####
#####################

#note: skip this step if not rarefying data

set.seed(24) #you must set a numerical seed like this for reproducibility
ps.dada2_join.rarefied <- rarefy_even_depth(ps.dada2_join, sample.size = 900) #rarefying samples so that if each sample had 900 reads (based on sample with lowest read, which is 900, then we could see how the abundance of different ASVs varies among samples in a standardized way


saveRDS(ps.dada2_join.rarefied, "ps.dada2_join.rarefied_900_20220108.RDS")

sort(sample_sums(ps.dada2_join.rarefied))


#############################################################
###########################################################
## saving the 900 rarefied data WITH 76 SAMPLES REMOVED as csv documents ##
hakai_2017_16s_900_20220108.otu <- as.data.frame(otu_table(ps.dada2_join.rarefied))
hakai_2017_16s_900_20220108.tax <- as.data.frame(tax_table(ps.dada2_join.rarefied))
hakai_2017_16s_900_20220108.sam <- as.data.frame(sample_data(ps.dada2_join.rarefied))


write.csv(hakai_2017_16s_900_20220108.otu, "hakai_2017_16s_900_RAREFIED_20220108.otu.csv")
write.csv(hakai_2017_16s_900_20220108.tax, "hakai_2017_16s_900_RAREFIED_20220108.tax.csv")
write.csv(hakai_2017_16s_900_20220108.sam, "hakai_2017_16s_900_RAREFIED_20220108.sam.csv")


## creating metadata file to use for Indval, nMDS, PERMANOVA, essentially anything except for the taxonomy plot because those plots should not use rarefied data ##

#############################################################################
# CREATING METADATA FILE FOR RAREFIED DATA 

# reading in sample_data file
metadata <- read.csv("hakai_2017_16s_900_RAREFIED_20220108.sam.csv")

View(metadata)
# selecting relevant columns only
metadata1 <- metadata[c(1,3, 5:9, 14, 17:20)]

# reading in otu table
otu <- read.csv("hakai_2017_16s_900_RAREFIED_20220108.otu.csv")

# joining metadata1 to otu via left join
meta_otu <- left_join(metadata1, otu, by = "X")

# remaning X as sample_id
meta_otu <- meta_otu %>%
  dplyr::rename(sample_id = X)

View(meta_otu)

# saving meta_otu as csv

write.csv(meta_otu, "hakai_2017_16s_meta_otu_900_RAREFIED_76removed_20220108.csv")

#############################################################
###########################################################
## saving the 900 rarefied data WITH ALL SAMPLES as csv documents ##
hakai_2017_16s_900_all_samples_20220108.otu <- as.data.frame(otu_table(ps.dada2_join.rarefied))
hakai_2017_16s_900_all_samples_20220108.tax <- as.data.frame(tax_table(ps.dada2_join.rarefied))
hakai_2017_16s_900_all_samples_20220108.sam <- as.data.frame(sample_data(ps.dada2_join.rarefied))


write.csv(hakai_2017_16s_900_all_samples_20220108.otu, "hakai_2017_16s_900_RAREFIED_all_samples_20220108.otu.csv")
write.csv(hakai_2017_16s_900_all_samples_20220108.tax, "hakai_2017_16s_900_RAREFIED_all_samples_20220108.tax.csv")
write.csv(hakai_2017_16s_900_all_samples_20220108.sam, "hakai_2017_16s_900_RAREFIED_all_samples_20220108.sam.csv")


## creating metadata file to use for Indval, nMDS, PERMANOVA - everything except for the taxonomy plot because those plots should not use rarefied data ##

#############################################################################
# CREATING METADATA FILE FOR RAREFIED DATA 

# reading in sample_data file
metadata <- read.csv("hakai_2017_16s_900_RAREFIED_all_samples_20220108.sam.csv")

View(metadata)
# selecting relevant columns only
metadata1 <- metadata[c(1,3, 5:9, 14, 17:20)]

# reading in otu table
otu <- read.csv("hakai_2017_16s_900_RAREFIED_all_samples_20211124.otu.csv")

# joining metadata1 to otu via left join
meta_otu <- left_join(metadata1, otu, by = "X")

# remaning X as sample_id
meta_otu <- meta_otu %>%
  dplyr::rename(sample_id = X)

View(meta_otu)

# saving meta_otu as csv

write.csv(meta_otu, "hakai_2017_16s_meta_otu_900_RAREFIED_all_samples_20220108.csv")



#############################################################
#### 76 SAMPLES REMOVED ####################################
###########################################################
## saving data WITHOUT rarefaction and AT RANK 6 as csv documents to use for frequency analysis  ##

# first selecting rank 6 only (genus level)
rank6 <- ps.dada2_join %>%
  tax_glom(taxrank = "Rank6")

hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.otu <- as.data.frame(otu_table(rank6))
hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.tax <- as.data.frame(tax_table(rank6))
hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.sam <- as.data.frame(sample_data(rank6))

write.csv(hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.otu, "hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.otu.csv")
write.csv(hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.tax, "hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.tax.csv")
write.csv(hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.sam, "hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.sam.csv")


# CREATING METADATA FILE 

# reading in sample_data file
metadata <- read.csv("hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.sam.csv")

#View(metadata)
# selecting relevant columns only
metadata1 <- metadata[c(1,3, 5:9, 14, 17:20)]

# reading in otu table
otu <- read.csv("hakai_2017_16s_NOT_RAREFIED_RANK_6_76removed_20220110.otu.csv")

# joining metadata1 to otu via left join
meta_otu <- left_join(metadata1, otu, by = "X")

# remaning X as sample_id
meta_otu <- meta_otu %>%
  dplyr::rename(sample_id = X)

#View(meta_otu)

# saving meta_otu as csv

write.csv(meta_otu, "hakai_2017_16s_meta_otu_NOT_RAREFIED_RANK_6_76removed_20220110.csv")


############################################################
### ALL SAMPLES INCLUDED ###################################
###########################################################
## saving data WITHOUT rarefaction and AT RANK 6 as csv documents to use for frequency analysis  ##

# first selecting rank 6 only (genus level)
rank6 <- ps.dada2_join %>%
  tax_glom(taxrank = "Rank6")

hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.otu <- as.data.frame(otu_table(rank6))
hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.tax <- as.data.frame(tax_table(rank6))
hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.sam <- as.data.frame(sample_data(rank6))

write.csv(hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.otu, "hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.otu.csv")
write.csv(hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.tax, "hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.tax.csv")
write.csv(hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.sam, "hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.sam.csv")


# CREATING METADATA FILE 

# reading in sample_data file
metadata <- read.csv("hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.sam.csv")

#View(metadata)
# selecting relevant columns only
metadata1 <- metadata[c(1,3, 5:9, 14, 17:20)]

# reading in otu table
otu <- read.csv("hakai_2017_16s_NOT_RAREFIED_RANK_6_all_samples_20220110.otu.csv")

# joining metadata1 to otu via left join
meta_otu <- left_join(metadata1, otu, by = "X")

# remaning X as sample_id
meta_otu <- meta_otu %>%
  dplyr::rename(sample_id = X)

#View(meta_otu)

# saving meta_otu as csv

write.csv(meta_otu, "hakai_2017_16s_meta_otu_NOT_RAREFIED_RANK_6_all_samples_20220110.csv")
