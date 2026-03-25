#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=400G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=800G # < if you are threading, ask for whole memory and not memory per CPU
#SBATCH --time=1-12:00:00  # 6 days, 0 hrs
#SBATCH --output=KOFamScan_Function_Annotation_Contigs_8.2.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="KOFamScan_Function_Annotation_Contigs_8.2.23"
#SBATCH -p highmem
# you can use any of the following: intel, batch, highmem, gpu

conda activate FunctionalAnnotation
#module load kofam_scan/1.3.0
# more /opt/linux/centos/7.x/x86_64/pkgs/kofam_scan/1.3.0/bin/config.yml -- where config.yml is stored
# For KOFam_Scan to work, you need the following installed: hmmer; parallel; ruby

# Create variable for date + time
day=$(date '+%m-%d-%y')
now=$(date '+%H:%M:%S')
 # %F is formal date (year-month-day)
 # %m is month; %d is day, %y is year, %H is hour, %M is minute, %S is second
 # %m is month, but %M is minute!
 
# Download KOfam database of HMM profiles
#if [[ ! -d /bigdata/aronsonlab/shared/RefDBs/kofam ]]; then
#    wget ftp://ftp.genome.jp/pub/db/kofam/*
#
#    gunzip /bigdata/aronsonlab/shared/RefDBs/kofam/ko_list.gunzip
#    tar xf /bigdata/aronsonlab/shared/RefDBs/kofam/profiles.tar.gz
#
#fi

# Make directory to store results
if [[ ! -d ./KOFamScan_Contig_Results ]]; then
    mkdir ./KOFamScan_Contig_Results
fi

# input into KOFam_Scan will be .faa files from Prokka or Prodigal (FASTA file of translated CDSs found in the contigs into Protein sequences)

#exec_annotation -o SSD_CoAssembly_functions_kfs.txt --profile=/bigdata/operations/pkgadmin/srv/projects/db/KEGG/97.0/profiles --ko_list=/bigdata/operations/pkgadmin/srv/projects/db/KEGG/97.0/ko_list --tmp-dir=SSD_CoAssembly_Contigs_KOFamScan -E 1.0e-5 -f detail-tsv --create-alignment SSD_Contigs_proteins.faa

# Contigs First
#for FILE in ./Contig_Protein_Seqs/*_proteins.faa;
#do
#    f=$(basename $FILE)
#    SAMPLE=${f%_proteins*}
#
#    if [[ ! -f ${SAMPLE}_functions_kfs.txt ]]; then
#        exec_annotation -o ${SAMPLE}_functions_kfs.txt --profile=/bigdata/operations/pkgadmin/srv/projects/db/KEGG/97.0/profiles --ko_list=/bigdata/operations/pkgadmin/srv/projects/db/KEGG/97.0/ko_list --tmp-dir=${SAMPLE}_KOFamScan -E 1.0e-5 -f detail-tsv --create-alignment ${FILE}
#    fi
#
#    if [[ ! -f ${SAMPLE}_functions_kfs_mapper.txt ]]; then
#        exec_annotation -o ${SAMPLE}_functions_kfs_mapper.txt --profile=/bigdata/operations/pkgadmin/srv/projects/db/KEGG/97.0/profiles --ko_list=/bigdata/operations/pkgadmin/srv/projects/db/KEGG/97.0/ko_list --tmp-dir=${SAMPLE}_KOFamScan.map -E 1.0e-5 -f mapper --create-alignment ${FILE}
#    fi
#
#    cp -R ${SAMPLE}_KOFamScan ${SAMPLE}_KOFamScan.map ${SAMPLE}_functions_kfs_mapper.txt ${SAMPLE}_functions_kfs.txt ./KOFamScan_Contig_Results/
#done

# Bins second

# Make directory to store results
if [[ ! -d ./KOFamScan_Bin_Results ]]; then
    mkdir ./KOFamScan_Bin_Results
fi

for FILE in ./Bin_Protein_Seqs/*_proteins.faa;
do
    f=$(basename $FILE)
    BIN=${f%_proteins*}
    
    if [[ ! -f ${BIN}_functions_kfs.txt ]]; then
        exec_annotation -o ${BIN}_functions_kfs.txt --profile=/bigdata/operations/pkgadmin/srv/projects/db/KEGG/97.0/profiles --ko_list=/bigdata/operations/pkgadmin/srv/projects/db/KEGG/97.0/ko_list --tmp-dir=${BIN}_KOFamScan -E 1.0e-5 -f detail-tsv --create-alignment ${FILE}
    fi
    
    if [[ ! -f ${BIN}_functions_kfs_mapper.txt ]]; then
        exec_annotation -o ${BIN}_functions_kfs_mapper.txt --profile=/bigdata/operations/pkgadmin/srv/projects/db/KEGG/97.0/profiles --ko_list=/bigdata/operations/pkgadmin/srv/projects/db/KEGG/97.0/ko_list --tmp-dir=${BIN}_KOFamScan.map -E 1.0e-5 -f mapper --create-alignment ${FILE}
    fi
    
    cp -R ${SAMPLE}_KOFamScan ${SAMPLE}_KOFamScan.map ${SAMPLE}_functions_kfs_mapper.txt ${SAMPLE}_functions_kfs.txt ./KOFamScan_Bin_Results/
done
conda deactivate

# Info on arguments: https://github.com/takaram/kofam_scan

    #exec_annotation -o POW_1_1a_S28_Bin_12.txt -E 1.0e-5 -f detail-tsv POW_1_1a_S28_Bin_12.faa

    #exec_annotation -o SeaWater_1a3_S27_kfs_functions.txt --tmp-dir=SeaWater_1a3_S27_tmp -E 1.0e-5 -f detail-tsv SeaWater_1a3_S27_proteins.faa

# Notes
# KOfam - HMM profiles for KEGG/KO with predefined score thresholds
# KOfamScan - Software to search KOfam
