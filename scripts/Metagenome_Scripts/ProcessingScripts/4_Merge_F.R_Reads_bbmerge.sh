#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --mem-per-cpu=500G
#SBATCH --mem=400G
#SBATCH --time=2-00:00:00     # 2 days
#SBATCH --output=Merge_NormReads_bbnorm_5.2.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Merging Normalized Reads 5/2/23"
#SBATCH -p aronsonlab
# you can use any of the following: intel, batch, highmem, gpu

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format
path=$("/rhome/hfreu002/SaltonSea/Metagenomes/SeqCenter_3.30.2023")

module load BBMap

if [[ ! -d ./Merged_Normalized_Seqs ]]; then
    mkdir ./Merged_Normalized_Seqs
fi

for i in ./Normalized_Seqs/*_R1_norm.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_R*}
    
    if [[ ! -f ./Merged_Normalized_Seqs/${SAMPLE}_norm_merged.fastq ]]; then # if the file doesn't exist in this directory, then proceed w/...
        
        bbmerge-auto.sh -Xmx400g in1=./Normalized_Seqs/${SAMPLE}_R1_norm.fastq in2=./Normalized_Seqs/${SAMPLE}_R2_norm.fastq out=./Merged_Normalized_Seqs/${SAMPLE}_norm_merged.fastq outu=./Merged_Normalized_Seqs/${SAMPLE}_norm_unmerged.fastq ihist=./Merged_Normalized_Seqs/${SAMPLE}_ihist_${today}.txt ecct extend2=20
    
        #cp ${SAMPLE}_norm_merged.fastq ./Merged_Normalized_Seqs
    
    fi

done

### *** confirm what the input files need to be and adjust file names for each in= ^^^

## Notes
# Basic merging (link = https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
#bbmerge.sh in=reads.fq out=merged.fq outu=unmerged.fq ihist=ihist.txt
#This will merge the reads by overlap. If no best overlap is found, the pair will go to outu; otherwise, the reads will be merged and sent to out. After finishing, an insert size histogram will be written to ihist.txt. This can be produced even if “out” or “outu” are not specified.

# Merging of nonoverlapping reads using kmers:
# bbmerge-auto.sh in=reads.fq out=merged.fq outu=unmerged.fq ihist=ihist.txt ecct extend2=20 iterations=5

# This will attempt to merge each pair by overlap. If unsuccessful, both reads will be error-corrected using Tadpole, and then merging will be tried again. If still unsuccessful, both reads will be extended by 20bp, then merging will be attempted again. This will repeat up to 5 times, or until neither of the reads can be extended any more due to a branch or dead-end in the kmer graph. If the reads are not merged, all of the changes are undone and the original pair will be sent to outu. “extend2=20 iterations=5” will extend each read by up to 100bp, which increases the maximum insert sizes that can be merged by 200bp. So, for example, a 2x150bp library can normally only merge inserts up to around 290bp; this would extend that capability to 490bp, and the middle would be filled in with assembled bases.

#     /Volumes/HLF_SSD/Aronson_Lab_Data/bbmap/bbmerge-auto.sh -Xmx500g in1=${SAMPLE}_L001_R1_norm.fastq in2=${SAMPLE}_L001_R2_norm.fastq out=${SAMPLE}_norm_merged.fastq outu=${SAMPLE}_norm_unmerged.fastq ihist=${SAMPLE}_ihist_${today}.txt ecct extend2=20
