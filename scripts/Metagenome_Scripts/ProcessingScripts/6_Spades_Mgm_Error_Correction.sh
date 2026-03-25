#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=200G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=400G
#SBATCH --time=14-00:00:00     # 14 days, 0 hrs
#SBATCH --output=read_error_correction_spades_5.29.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Error correcting clean reads with Spades 5/29/2023"
#SBATCH -p aronsonlab

# you can use any of the following: intel, batch, highmem, gpu

Path="/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023"

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

module load spades

# First error correct normalized sequences
#for i in ./Normalized_Seqs/*_R1_norm.fastq;
#do
#    f=$(basename $i)
#    SAMPLE=${f%_R*}
#    
#    if [[ ! -d ./${SAMPLE}_norm_error_corrected ]]; then
#        spades.py -1 ./Normalized_Seqs/${SAMPLE}_R1_norm.fastq -2 ./Normalized_Seqs/${SAMPLE}_R2_norm.fastq -o ${SAMPLE}_norm_error_corrected -t 4 --meta --only-error-correction
#        
#    fi
    
#    #for j in ${Path}/${SAMPLE}_norm_error_corrected/corrected/*.fastq.gz;
#    #    rename norm_EC_EC.00.0_0.cor norm_EC $j
#    #    #cp *.fastq.gz ${Path}/
#    #done
#    
#    # to include merged reads in addition to F & R reads, include --merged ${SAMPLE}_norm_merged.fastq
#    #echo -e "Finished \aerror correction with ${SAMPLE}_R1 ${SAMPLE}_R2"
#
#done

# Then error correct trimmed, non-normalized sequences
for i in ./Trimmed_Seqs/*_R1_clean.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_R*}
    
    if [[ ! -d ./${SAMPLE}_error_corrected ]]; then
        spades.py -1 ./Trimmed_Seqs/${SAMPLE}_R1_clean.fastq -2 ./Trimmed_Seqs/${SAMPLE}_R2_clean.fastq -o ${SAMPLE}_error_corrected -t 4 --meta --only-error-correction
        
    fi
    
    #for j in ${Path}/${SAMPLE}_error_corrected/corrected/*.fastq.gz;
    #    rename clean_EC_EC.00.0_0.cor clean_EC $j
    #    #cp *.fastq.gz ${Path}/
    #done
    
    # to include merged reads in addition to F & R reads, include --merged ${SAMPLE}_norm_merged.fastq
    #echo -e "Finished \aerror correction with ${SAMPLE}_R1 ${SAMPLE}_R2"

done


# Copy files to home directory -- change path in the cp command
#for i in *_error_corrected;
#do
#    f=$(basename $i)
#    SAMPLE=${f%_error_corrected*}
#
    #if [[ ! -d ./${SAMPLE}_error_corrected ]]; then
    #    spades.py -1 ${SAMPLE}_R1_norm.fastq -2 ${SAMPLE}_R2_norm.fastq -o ${SAMPLE}_error_corrected -t 4 --meta --only-error-correction
        
   # fi

#   cp ${SAMPLE}_error_corrected/corrected/*.fastq.gz /bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_7.12.2022/Sequences_Analysis/
    
    
    # to include merged reads in addition to F & R reads, include --merged ${SAMPLE}_norm_merged.fastq
    #echo -e "Finished \aerror correction with ${SAMPLE}_R1 ${SAMPLE}_R2"

#done

#rename _unpaired.00.0_0.cor _unpaired_EC *_unpaired.00.0_0.cor.fastq.gz

## Spades assembly - https://github.com/ablab/spades#meta
## Notes for Spades...
## spades.py -1 read1.fq -2 read2.fq --merged merged.fq -o spades_test <<< -s is the flag to indicate merged reads as a "library"
## Error corrected sequences wil be found in output_directory/corrected/ – files from read error correction

## Spades Output
#The full list of <output_dir> content is presented below:

#scaffolds.fasta – resulting scaffolds (recommended for use as resulting sequences) *****
#contigs.fasta – resulting contigs
#assembly_graph.fastg – assembly graph
#contigs.paths – contigs paths in the assembly graph
#scaffolds.paths – scaffolds paths in the assembly graph
#before_rr.fasta – contigs before repeat resolution
#corrected/ – files from read error correction
#   configs/ – configuration files for read error correction
#   corrected.yaml – internal configuration file
#   Output files with corrected reads
#params.txt – information about SPAdes parameters in this run
#spades.log – SPAdes log
#dataset.info – internal configuration file
#input_dataset.yaml – internal YAML data set file
#K<##>/– directory containing intermediate files from the run with K=<##>. These files should not be used as assembly results; use resulting contigs/scaffolds in files mentioned above.

## HPCC (Cluster systems) note
##SBATCH --mem-per-cpu=500G # --mem=900gb --> how to request for total allocation of mem rather than mem per cpu; include in job submission command


