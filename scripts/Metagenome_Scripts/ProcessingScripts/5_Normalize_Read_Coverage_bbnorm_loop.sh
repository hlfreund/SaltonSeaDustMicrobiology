#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=12 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=500G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=400G
#SBATCH --time=7-00:00:00     # 7 days, 0 hrs
#SBATCH --output=Normalize_Read_Coverage_4.25.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Normalize Read Coverage"
#SBATCH -p highmem

# you can use any of the following: intel, batch, highmem, gpu

Path="/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023"

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

module load BBMap/38.95

if [[ ! -d ./Normalized_Seqs ]]; then
    mkdir ./Normalized_Seqs
fi

for i in ./Trimmed_Seqs/*_R1_clean.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_R*}
    echo "Normalizing trimmed reads from" ${SAMPLE}_R1_clean.fastq ${SAMPLE}_R2_clean.fastq
    
    if [[ ! -f ./Normalized_Seqs/${SAMPLE}_R1_norm.fastq ]] && [[ ! -f ./Normalized_Seqs/${SAMPLE}_R2_norm.fastq ]]; then
        bbnorm.sh in1=${Path}/Trimmed_Seqs/${SAMPLE}_R1_clean.fastq in2=${Path}/Trimmed_Seqs/${SAMPLE}_R2_clean.fastq out1=${Path}/Normalized_Seqs/${SAMPLE}_R1_norm.fastq out2=${Path}/Normalized_Seqs/${SAMPLE}_R2_norm.fastq target=100 min=5
        # This will run 2-pass normalization to produce an output file of reads with an average depth of 100x. Reads with an apparent depth of under 5x will be presumed to be errors and discarded.
    
    
        #cp ${SAMPLE}_R1_norm.fastq ${SAMPLE}_R2_norm.fastq ${Path}/Normalized_Seqs/

    fi
        
done

# ** Normalization happens AFTER trimming reads
# Link = https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbnorm-guide/

#bbnorm.sh in=reads.fq out=normalized.fq target=100 min=5

## Notes
# BBNorm is designed to normalize coverage by down-sampling reads over high-depth areas of a genome, to result in a flat coverage distribution. This process can dramatically accelerate assembly and render intractable datasets tractable, and often improve assembly quality. It can also do depth-binning, kmer frequency histogram generation, error-correction, error-marking, and genome-size estimation. BBNorm has 4 particularly notable features:
# 1) It stores kmers in a probabilistic data structure called a count-min sketch. This means it will never run out of memory, or swap to disk, on any dataset. Rather, as the number of unique kmers increases, accuracy gradually declines.
# 2) It has numerous features such as multipass normalization, which reduce the average error rate in the normalized output; whereas standard normalization enriches for reads containing errors.
# 3) It is extremely fast and easy-to-use compared to other normalization programs.
# 4) It supports unlimited kmer lengths.

## When Not To Use BBNorm:
# For normalization, BBNorm is mainly intended for use in assembly, and with short reads. Normalization is often useful if you have too much data (for example, 600x average coverage when you only want 100x) or uneven coverage (amplified single-cell, RNA-seq, viruses, metagenomes, etc). It is not useful if you have smooth coverage and approximately the right amount of data, or too little data. BBNorm cannot inflate low coverage (bring 15x coverage up to 100x), only reduce it. Never normalize read data prior to a quantitative analysis (like ChIP-seq, RNA-seq for expression profiling, etc); if you assemble normalized data, and want to use mapping to determine coverage, map the non-normalized reads. Also, do not normalize data prior to mapping for variant discovery; it will cause bias. If you need to reduce data volume in any of these scenarios, use subsampling rather than normalization. Do not attempt to normalize high-error-rate data from platforms such as PacBio or Nanopore; it is designed for relatively-low-error-rate, short, fixed-length reads such as Illumina and Ion Torrent.
# Also, error-correction is not advisable when you are looking for rare variants. It should generally be fine with relatively high-depth coverage of heterozygous mutations in a diploid (where you expect a 50/50 allele split), but with low-depth coverage (like 5x), or very lopsided distributions (like a 1/100 allele split), it may correct the minority allele into the majority allele, so should be used with caution.





