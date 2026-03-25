#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=500G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=700G # run this line if you are threading; request bulk memory
#SBATCH --time=10-00:00:00     # 7 days, 0 hrs
#SBATCH --output=SSD_ContigAssembly_megahit.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Assembling MGMs with merged + unmerged reads"
#SBATCH -p highmem
# you can use any of the following: intel, batch, highmem, gpu

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

module load megahit

# move the files you want to use for assembly to directory where you'll run megahit
# this way we can use the list variable without having to worry about the directory
F_list=$(ls *_R1_clean_EC.fastq | sed -e '$ ! s/$/,/g' | tr -d '\n') # list of F read files, separated by commas except no comma on the end of the list
R_list=$(ls *_R2_clean_EC.fastq | sed -e '$ ! s/$/,/g' | tr -d '\n') # list of R read files, separated by commas except no comma on the end of the list
M_list=$(ls *_norm_merged.fastq | sed -e '$ ! s/$/,/g' | tr -d '\n') # list of merged read files, separated by commas except no comma on the end of the list

if [[ ! -d ./Contig_Assembly_megahit ]]; then
    megahit -1 ${F_list} -2 ${R_list} -r ${M_list} -t 12 -o Contig_Assembly_megahit
fi

    
    # to include merged reads in addition to F & R reads, include --merged ${SAMPLE}_norm_merged.fastq
    #echo -e "Finished \aerror correction with ${SAMPLE}_R1 ${SAMPLE}_R2"

# sed tip at link below
## https://www.quora.com/What-does-sed-w*-d-mean-in-Unix#:~:text='d'%20in%20sed%20means%20delete,for%20regex%20in%20pattern%20space.

### MEGAHIT assembly tips
# for ultra complex metagenomics data such as soil, a larger kmin, say 27, is recommended to reduce the complexity of the de Bruijn graph.
### Quality trimming is also recommended
# --min-count 2 ----> (kmin+1)-mer with multiplicity lower than d (default 2, specified by --min-count option) will be discarded;
### recommend using the default value 2 for metagenomics assembly

## MORE MEGAHIT TIPS: https://sites.google.com/site/wiki4metagenomics/tools/assembly/megahit
