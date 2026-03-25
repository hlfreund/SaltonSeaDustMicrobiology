#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=1 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=400G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=700G # < if you are threading, ask for whole memory and not memory per CPU
#SBATCH --time=5-00:00:00     # 5 days
#SBATCH --output=GTDBtk_Taxonomic_Annotation_8.17.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="GTDBtk_taxa_annotation_8.17.23"
#SBATCH -p highmem

# you can use any of the following: intel, batch, highmem, gpu

today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

module load gtdbtk/2.1.0 # load the module on your computing cluster system (ignore this if running locally)

# Make directory to store results
if [[ ! -d ./Taxo_Annotation ]]; then
    mkdir ./Taxo_Annotation
    mkdir ./Taxo_Annotation/GoodBins_Results
    #mkdir ./Taxo_Annotation/MedBins_Results
fi

for dir in *_good_bins;
do
    file=$(basename $dir)
    SAMPLE=${file%*}
    echo $SAMPLE

    for bin in ${dir}/*;
    do

        bin_fasta=$(basename $bin)
        bin_num=${bin_fasta%*.fa}
        if [[ ! -d ./gtdbtk_${bin_num} ]]; then
            gtdbtk classify_wf --genome_dir $dir -x fa --out_dir taxo.annot_${SAMPLE} --prefix taxo_annot
            mv taxo.annot_${SAMPLE} ./Taxo_Annotation/GoodBins_Results
        fi
        #SAMPLE=${file%_good*}
        #echo $bin_fasta $bin_num
    done

done

#for SAMPDIR in $(ls -d *_good_bins)
#do
#    SAMPLE=$(basename $SAMPDIR _good_bins)
#    if [[ ! -d ./gtdbtk_${SAMPLE} ]]; then
#        gtdbtk classify_wf --genome_dir $SAMPDIR -x fa --out_dir gtdbtk_${SAMPLE} --prefix taxo_annot
#       mv gtdbtk_${SAMPLE} ./Taxo_Annotation/GoodBins_Results
#    fi
 #echo $SAMPLE
#done
#for SAMPDIR in $(ls -d *_med_bins)
#do
#    SAMPLE=$(basename $SAMPDIR _med_bins)
#    if [[ ! -d gtdbtk_med.bin_${SAMPLE} ]]; then
#        gtdbtk classify_wf --genome_dir $SAMPDIR -x fa --out_dir gtdbtk_med.bin_${SAMPLE} --prefix taxo_annot
#        mv gtdbtk_med.bin_${SAMPLE} ./Taxo_Annotation/MedBins_Results
#    fi
# #echo $SAMPLE
#done

#find ./*_bins_metabat | while read i;
#do
#x=$(echo $i)
#y=${x%_Bin*}
#SAMPLE=$(echo ${y##*/})
##echo $SAMPLE
#gen_dir=$(echo $(dirname ${i}))
##echo ${gen_dir}
#gtdbtk classify_wf --genome_dir ${gen_dir} -x fa --out_dir gtdbtk_${SAMPLE}_${today}

#done
#gtdbtk classify_wf --genome_dir EA_Pool-SeaWater_1a3_S27_bins_metabat -x fa --out_dir gtdbtk_SeaWater_1a3_S27_${today}


for SAMPDIR in $(ls -d *_good_bins)
do
    SAMPLE=$(basename $SAMPDIR _good_bins)
    echo $SAMPLE
done
