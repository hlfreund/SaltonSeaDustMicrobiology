#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=10G  # max you can get in batch nodes
#SBATCH --mem=400G # < if you are threading, ask for whole memory and not memory per CPU
#SBATCH --time=3-00:00:00     # 3 days
#SBATCH --output=16S.V3V4_SSD_PhyloTrees.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="16S.V3V4_SSD_PhyloTree"
#SBATCH -p aronsonlab
# you can use any of the following: intel, batch, highmem, gpu

# You must have generated a clean, updated FASTA file only containing the ASVs you are studying, aka no contaminant or eukaryotic ASVs
## if you do not have this yet, go back and run 5a_Prep_Data_for_MSA_PhylogeneticTree.R

# First load all modules you will need to generate multiple sequence alignment (MSA) and build phylogenetic tree!
module load ssu-align hmmer fasttree iqtree

# Step1: run the multiple sequence alignment of 16S V3V4 rRNA with ssu-align
ssu-align 16S.V3V4_ASVs_dada2_clean.fa SSD_16S.V3V4_MSA_ssualign

# Step 2: convert Stockholm aligned file (.stk) into a multifasta using esl-reformat which is part of the easel toolkit within hmmer
esl-reformat afa ./SSD_16S.V3V4_MSA_ssualign/SSD_16SV3V4_ssualign.bacteria.stk > ./SSD_16S.V3V4_MSA_ssualign/SSD_16SV3V4_ssualign.bacteria.fasaln

# Step 3: run FastTree to build phylogenetic tree from 16S V3V4 MSA
FastTreeMP -nt -gtr -gamma < ./SSD_16S.V3V4_MSA_ssualign/SSD_16SV3V4_ssualign.bacteria.fasaln > ./SSD_16S.V3V4_MSA_ssualign/SSD_16SV3V4_ssualign.bacteria.FastTree.tre

# Step 4 (Optional) : run iqtree to build phylogenetic tree
iqtree2 -s ./SSD_16S.V3V4_MSA_ssualign/SSD_16SV3V4_ssualign.bacteria.fasaln -nt AUTO -B 1000 --alrt 1000

# After you have finished running this, you can now generate Unifrac distance!
