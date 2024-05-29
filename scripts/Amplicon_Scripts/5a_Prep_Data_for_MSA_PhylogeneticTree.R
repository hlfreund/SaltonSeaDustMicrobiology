#### Load Packages ####

suppressPackageStartupMessages({ # load packages quietly
  library(devtools)
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  #library(scales)
  library(grid)
  library(ape)
  library(plyr)
  library(dplyr)
  library(viridis)
  #library(ggbiplot)
  library(readxl)
  #library(metagenomeSeq)
  library(DESeq2)
  library(dplyr)
  library(magrittr)
  library(MASS)
  library(dendextend)
  library(tidyr)
  library(reshape)
  library(reshape2)
  #library(wesanderson)
  #library(nationalparkcolors)
  library(shades)
  #library(ALDEx2)
  library(rstatix)
  #library(decontam)
  #library(ggvegan)
  #library(microbiome)
  #library(pairwiseAdonis)
  library(corrplot)
  #library(fst)
  library(plotly)
  library(htmlwidgets)
  library(DECIPHER)
  library(phangorn)
})

# NOTE: Run this script before generating Multiple Sequence Alignment (MSA) and Phylogenetic Tree!
## point of this script is to generate a new FASTA file that contains only the ASVs you are studying (aka no contaminant or eukaryotic ASVs included)

#### Load Global Env to Import Count/ASV Tables ####

# load your ready-to-go Rdata file that you output from 1a_Prep_Amplicon_Data_All.R
## aka the Rdata file you load at the start of every analysis script you run - should include clean, decontaminated ASV table!
load("SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file

# load your DADA2 output for that specific run - must have finished running DADA2 pipeline to use this!
load("mydada_16S.V3V4.Rdata")

#### Pull out ASV Sequences of Interest for MSA & Phylogenetic Tree ####
## pulling sequences from DADA2 results only for ASVs that we have in our cleaned up bac.ASV_table
## use this new set of ASVs of interest to construct phylogenetic tree, import the tree with ape package, and calculate Unifrac distance with rbiom
# giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab.nochim) # pull out ASV sequences
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character") # pull out ASV headers

# add ASV in front of each # of each sequence
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste("ASV", i, sep="_")
}

# create table of ASVs and their respective sequences so that we can filter out ASVs from samples we aren't looking at
asv_all <- as.data.frame(cbind(asv_headers, asv_seqs))
length(unique(asv_all$asv_headers))

# only pull out ASVs and their sequences for ASVs in our taxa table for only samples of interest
new.asv_all<-asv_all[asv_all$asv_headers %in% colnames(bac.ASV_table[,-1]),]

# insert ">" before ASV IDs before we make the new FASTA file
new.asv_all$asv_headers<-gsub("ASV_",">ASV_",new.asv_all$asv_headers)
new.asv_all$asv_headers[1:5]

# create and save the new fasta file - will use this to generate Multiple Sequence Alignment & build phylogenetic tree
new.asv_fasta<-c(rbind(new.asv_all$asv_headers,new.asv_all$asv_seqs))
write(new.asv_fasta, "SSD_16SV3V4_ASVs_Updated.fa")

# NOTE: Now you can use SSD_16SV3V4_ASVs_Updated.fa to generate MSA and build phylogenetic tree in script 5b_MultipleSequenceAlignment_and_PhylogeneticTree.sh
