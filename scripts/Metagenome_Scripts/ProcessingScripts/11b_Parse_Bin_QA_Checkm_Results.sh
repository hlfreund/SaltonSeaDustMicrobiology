#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=400G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
##SBATCH --time=01:00:00     # 1 day
#SBATCH --output=Parse_Bin_QA_Checkm_Results_6.29.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Parse_Bin_QA_Checkm_Results_6.29.23"
#SBATCH -p aronsonlab
# you can use any of the following: intel, batch, highmem, gpu

# Need to look through CheckM results to examine completeness and contamination of (assembled) genome bins
# Separate out bins based on completeness and contamination

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

# double check the dates given on the names of the files before you start!
for i in *_checkm_qa.tsv; # change to your output file names from when you ran CheckM
do
    file=$(basename $i)
    SAMPLE=${file%_checkm*} # double

    echo "${SAMPLE}: separating bins based on quality -- CheckM results"
    
    #sed -i 's/_contigs.fasta.metabat-bins-_-t_4_-m_1500./_bin_/g' ${i}
    if [[ ! -f ${SAMPLE}_good_bins.tsv ]]; then
        awk -F '\t' 'BEGIN {OFS="\t"} { if ( NR == 1 ) print "Sample_ID", "Bin_Num" ,"Taxa_Level", "Marker_Lineage", "Lineage_ID", "Genome_Num", "Completeness", "Contamination", "Strain_Heterogeneity", "GC_Content"; else if ($6 >=80 && $7 < 5) print $1,$2,$3,$6,$7,$8,$19}' ${i} | sed 's/ /\t/g' | sed 's/_bin/\tBin/g' | sed 's/__/\t/g' > ${SAMPLE}_good_bins.tsv
    fi
    
    ## Bins over >80% complete, <5% contamination for annotation
    
    if [[ ! -f ${SAMPLE}_med_bins.tsv ]]; then
        awk -F '\t' 'BEGIN {OFS="\t"} { if ( NR == 1 ) print "Sample_ID", "Bin_Num" ,"Taxa_Level", "Marker_Lineage", "Lineage_ID", "Genome_Num", "Completeness", "Contamination", "Strain_Heterogeneity", "GC_Content"; else if ($6 >=50 && $7 < 5) print $1,$2,$3,$6,$7,$8,$19}' ${i} | sed 's/ /\t/g' | sed 's/_bin/\tBin/g' | sed 's/__/\t/g' > ${SAMPLE}_med_bins.tsv
    fi
    ## Bins over >50% complete, <5% contamination for annotation

    
done

for SAMPDIR in $(ls -d *_bins);
do
    SAMPLE=$(basename $SAMPDIR _bins)
    
    if [[ ! -d ${SAMPLE}_good_bins ]] && [[ ! -d ./${SAMPLE}_med_bins ]]; then
        mkdir ${SAMPLE}_good_bins
        mkdir ${SAMPLE}_med_bins
    fi
    
    #echo $ffa
    echo "Copying quality bins (>80% completeness, <5% contamination) into ${SAMPLE}_good_bins"
    
    for f in $(awk 'FNR == 1 {next} NR > 1{ print $1 }' ${SAMPLE}_good_bins.tsv);
    do
         ffa=$(echo $f.fa)
         #mkdir ${SAMPLE}_quality_bins
         cp "${SAMPDIR}/${ffa}" "./${SAMPLE}_good_bins"
         
         # if you want to keep track of the date in which you ran this step, just add ${today} to your ${SAMPLE}_good_bins directory name

         #echo $ffa
    done
    
    echo "Copying medium quality bins (>50% completeness, <5% contamination) into ${SAMPLE}_med_bins"
    
    for f in $(awk 'FNR == 1 {next} NR > 1{ print $1 }' ${SAMPLE}_med_bins.tsv);
    do
         ffa=$(echo $f.fa)
         #mkdir ${SAMPLE}_quality_bins
         cp "${SAMPDIR}/${ffa}" "./${SAMPLE}_med_bins"
         
         # if you want to keep track of the date in which you ran this step, just add ${today} to your ${SAMPLE}_good_bins directory name

         #echo $ffa
    done
    
done
   
   # Rename actual bins so you know which sample they came from
for dir in *_good_bins;
do
    file=$(basename $dir)
    SAMPLE=${file%_good*}

    for bin in ${dir}/*;
    do

        bin_name=$(basename $bin)
        #SAMPLE=${file%_good*}
        rename bin ${SAMPLE}_bin ./${dir}/${bin_name}
    done

done

for dir in *_med_bins;
do
    file=$(basename $dir)
    SAMPLE=${file%_med*}

    for bin in ${dir}/*;
    do

        bin_name=$(basename $bin)
        #SAMPLE=${file%_med*}
        rename bin ${SAMPLE}_bin ./${dir}/${bin_name}
    done

done

## Thresholds taken from the following papers:
# Chen, L. X., Anantharaman, K., Shaiber, A., Murat Eren, A., & Banfield, J. F. (2020). Accurate and complete genomes from metagenomes. Genome Research. Cold Spring Harbor Laboratory Press. https://doi.org/10.1101/gr.258640.119
# Bowers, R., Kyrpides, N., Stepanauskas, R. et al. Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea. Nat Biotechnol 35, 725–731 (2017). https://doi.org/10.1038/nbt.3893

# AWK notes
# -F ' ' == field separator of input doc is a space
# NR == 1 == if first line is 1, do _____
