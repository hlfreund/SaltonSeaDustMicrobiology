#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=400G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=700G # < if you are threading, ask for whole memory and not memory per CPU
#SBATCH --time=7-00:00:00  # 7 days, 0 hrs
#SBATCH --output=eggNOG_mapper_Function_Annotation_6.30.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="eggNOG_mapper_Function_Annotation_6.30.23"
#SBATCH -p highmem
# you can use any of the following: intel, batch, highmem, gpu

module load eggnog-mapper
module load hmmer

# Make directory to store results
if [[ ! -d ./eggNOG_FxnAnnotation_Results ]]; then
    mkdir eggNOG_FxnAnnotation_Results
fi

# input into eggNOG-mapper will be .faa files from Prokka (Protein FASTA file of the genome bins)
emapper.py -m mmseqs --itype metagenome -i SSD_Contigs_proteins.faa -o SSD_CoAssembly --cpu 1 --genepred prodigal --temp_dir ./eggNOG_FxnAnnotation_Results/SSD_CoAssembly_EN_tmp --output_dir ./eggNOG_FxnAnnotation_Results

for FILE in ./Bin_Protein_Seqs/*_proteins.faa;
do
    f=$(basename $FILE)
    BIN=${f%_proteins*}
    
    if [[ ! -d ./eggNOG_FxnAnnotation_Results/${BIN}_EN_tmp ]]; then
        mkdir eggNOG_FxnAnnotation_Results/${BIN}_EN_tmp
    fi
    emapper.py -m mmseqs --itype metagenome -i ${FILE} -o ${BIN}_EN --cpu 1 --genepred prodigal --temp_dir ./eggNOG_FxnAnnotation_Results/${BIN}_EN_tmp --output_dir ./eggNOG_FxnAnnotation_Results
done

# More info about program usage here: https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.2-to-v2.1.4#Basic_usage

