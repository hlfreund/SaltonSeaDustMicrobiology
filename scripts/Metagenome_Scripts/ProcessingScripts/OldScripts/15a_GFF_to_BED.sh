#!/bin/bash -l

#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=500G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=10G
#SBATCH --time=1-00:00:00     # 1 days, 0 hrs
#SBATCH --output=GFF_to_BED_8.21.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="GFF_to_BED_file_8.21.23"
#SBATCH -p batch

# Convert GFF files to BED files
# Run this step before featureCounts; need BED file for featureCounts

for i in *_genes.gff; # change to your output file names from when you ran Prodigal
do
    file=$(basename $i)
    SAMPLE=${file%_genes*}
    
    
    if [[ ! -f ${SAMPLE}_genes2.bed ]]; then
        awk -F '\t' '/#/ {next} {OFS = FS} {print $1=$1 "_" (++count[$1]),$4,$5,$1}' ${i} > ${SAMPLE}_genes2.bed
        mv *_genes2.bed ./GeneCoverage_Results/BED_Reformatted
    fi
    
    # $1=$1 "_" (++count[$1]) -> count every unique ID and add counter after "_" and original ID name in column 1 ($1)
done

# Format of BED file (https://genome.ucsc.edu/FAQ/FAQformat.html):
## has at least 3-4 columns, up to 12
## col1: name of contig; col2: start position of unique feature in contig; col3: end position of unique feature in contig; col4: name or ID of unique feature
# Format of GFF file (https://genome.ucsc.edu/FAQ/FAQformat.html#format3):
## 9 columns
## col1: seq or contig name; col2: program that generated feature; col3: name of the feature; col4: feature start; col5: feature end; col6: feature score; col7: strand origin; col8: reading frame; col9: group (?)

