#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=200G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=400G
##SBATCH --time=4-00:00:00     # 4 days, 0 hrs
#SBATCH --output=read_error_correction_spades_6.6.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Error correcting clean reads with Spades 6/6/2023"
#SBATCH -p aronsonlab

# you can use any of the following: intel, batch, highmem, gpu

Path="/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023"

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

module load spades

# First error correct trimmed, non-normalized sequences
for i in ./Normalized_Seqs/*_R1_norm.fastq;
do
    f=$(basename $i)
    SAMPLE=${f%_R*}
    
    if [[ ! -d ./${SAMPLE}_error_corrected ]]; then
        spades.py -1 ./Normalized_Seqs/${SAMPLE}_R1_norm.fastq -2 ./Normalized_Seqs/${SAMPLE}_R2_norm.fastq -o ${SAMPLE}_norm_error_corrected -t 4 --meta --only-error-correction
        
    fi
    
    # to include merged reads in addition to F & R reads, include --merged ${SAMPLE}_norm_merged.fastq
    #echo -e "Finished \aerror correction with ${SAMPLE}_R1 ${SAMPLE}_R2"

done

# Then rename all files; first normalized error corrected
for i in *_norm_error_corrected/corrected/*.fastq.gz;
do
    #echo $(basename $i)
        
    rename norm.00.0_0.cor norm_EC $i
    rename unpaired.00.0_0.cor unpaired_norm_EC $i
        #cp *.fastq.gz ${Path}/
done

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

for i in *_error_corrected/corrected/*.fastq.gz;
do
    #echo $(basename $j)
    gunzip $i
    cp $i /rhome/hfreu002/SaltonSea/Metagenomes/SeqCenter_3.30.2023/ContigAssembly/
        

done
