#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=900G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=800G
#SBATCH --time=10-00:00:00     # 10 days, 0 hrs
#SBATCH --output=Check_QA_assembled_contigs_6.28.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Check_QA_assembled_contigs_6.28.23"
#SBATCH -p highmem

# you can use any of the following: intel, batch, highmem, gpu

#cd /bigdata/aronsonlab/shared/HannahFreund/HannahFTemp/ # change to your desired directory


today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

#module load checkm/1.1.3 # load the module on your computing cluster system (ignore this if running locally)
#

conda activate MGM_Quality

if  [[ ! -d ./Checkm_Results ]]; then
    mkdir Checkm_Results
fi

if  [[ ! -f ./Bacteria_marker_file ]] && [[ ! -f ./Archaea_marker_file ]]; then
    checkm taxon_set domain Bacteria Bacteria_marker_file
    checkm taxon_set domain Archaea Archaea_marker_file
fi


for SAMPDIR in $(ls -d *_bins);
do
    SAMPLE=$(basename $SAMPDIR _bins)
    
    if [[ ! -f ./${SAMPLE}_checkm_lineage/${SAMPLE}_checkm_lineage_results.txt ]]; then
        checkm lineage_wf -t 8 -x fa ./${SAMPLE}_bins ./${SAMPLE}_checkm_lineage --file ${SAMPLE}_checkm_lineage_results.txt --tab_table
        ## checkm lineage_wf -t {#threads} -x {file format, aka fa} {BIN_DIRECTORY} {OUTPUT_DIRECTORY} -f {output file} --tab_tablels -
        mv ${SAMPLE}_checkm_lineage_results.txt ./Checkm_Results
    fi
    
    if [[ ! -f ./${SAMPLE}_checkm_qa.tsv ]]; then
        checkm qa -o 2 -f ./${SAMPLE}_checkm_qa.tsv --tab_table ./${SAMPLE}_checkm_lineage/lineage.ms ./${SAMPLE}_checkm_lineage/
        # checkm qa --output_format 2 -f {output_file} --table_table {/path/to/marker/file(either lineage.ms if used checkm lineage_wf, or a taxon specific ms used for analyze command} {/path/to/checkm_lineage (or taxonomy) output}
        
        cp ${SAMPLE}_checkm_qa.tsv ./Checkm_Results
    fi
    
    if [[ ! -f ./${SAMPLE}_checkm_coverage.tsv ]]; then
        checkm coverage ./${SAMPLE}_bins ${SAMPLE}_checkm_coverage.tsv ${SAMPLE}_contigs_sorted.bam
        ## NOTE: Need to have original BAM file, sorted BAM file, and indexed bam.bai file all in same directory when running coverage command (https://github.com/Ecogenomics/CheckM/issues/71)
        mv ${SAMPLE}_checkm_coverage.tsv ./Checkm_Results
    fi
    
    if [[ ! -d ./${SAMPLE}_bacteria_markers ]] && [[ ! -d ./${SAMPLE}_archaea_markers ]]; then
        checkm analyze Bacteria_marker_file ./${SAMPLE}_bins -x fa -t 4 ./${SAMPLE}_bacteria_markers
        checkm qa Bacteria_marker_file ./${SAMPLE}_bacteria_markers -o 2 -f ${SAMPLE}_bacteria_marker_results.txt --tab_table
        

        checkm analyze Archaea_marker_file ./${SAMPLE}_bins -x fa -t 4 ./${SAMPLE}_archaea_markers
        checkm qa Archaea_marker_file ./${SAMPLE}_archaea_markers -o 2 -f ${SAMPLE}_archaea_marker_results.txt --tab_table
        
        mv ${SAMPLE}_bacteria_marker_results.txt ${SAMPLE}_archaea_marker_results.txt ${SAMPLE}_archaea_markers ${SAMPLE}_bacteria_markers ./Checkm_Results
    fi
done

conda deactivate


## Info on CheckM
## https://github.com/Ecogenomics/CheckM/wiki/Workflows#lineage-specific-workflow
## https://github.com/Ecogenomics/CheckM/wiki/
## https://www.hadriengourle.com/tutorials/meta_assembly/ -- using it
## https://bitbucket.org/berkeleylab/metabat/wiki/Best%20Binning%20Practices -- using it after MetaBAT

