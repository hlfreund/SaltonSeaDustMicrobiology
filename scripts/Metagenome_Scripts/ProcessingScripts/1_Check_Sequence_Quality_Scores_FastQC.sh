#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100G  # max you can get in batch nodes
#SBATCH --time=6:00:00     # 6 hours
#SBATCH --output=Check_Seq_Quality_Scores_4.2.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Check Paired Sequence Quality Scores 4/2/2023"
#SBATCH -p aronsonlab # This is the default partition, you can use any of the following: intel, batch, highmem, gpu

module load fastqc

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

# Use FastQC: command fastqc, then fastq or fasta files, then -o /ouputdirectory/
if [[ ! -d ./FastQC_Results_PreTrim ]]; then
    mkdir FastQC_Results_PreTrim
fi

for i in *_R1.fastq.gz;
do
    file=$(basename $i)
    SAMPLE=${file%_R*}
    echo ${SAMPLE} "Checking quality of sequences with FastQC"
    
    fastqc ${SAMPLE}_R1.fastq.gz -o FastQC_Results_PreTrim
    fastqc ${SAMPLE}_R2.fastq.gz -o FastQC_Results_PreTrim
    
done
