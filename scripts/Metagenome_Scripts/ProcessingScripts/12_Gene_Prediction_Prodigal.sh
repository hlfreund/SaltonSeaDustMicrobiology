#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=1 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=400G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
##SBATCH --mem=400G # < if you are threading, ask for whole memory and not memory per CPU
##SBATCH --time=20:00:00     # 20 hrs
#SBATCH --output=Prodigal_Functional_Annotation_6.29.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Prodigal_functional_annotation_6.29.23"
#SBATCH -p aronsonlab

# you can use any of the following: intel, batch, highmem, gpu

module load prodigal
#today=$(date "+%m.%d.%y") # date is the command to get today's date, and the "+%m_%d_%y" will print it in month_day_year format

if [[ ! -d ./Prodigal_Results_Contigs ]]; then
    mkdir ./Prodigal_Results_Contigs
fi

if [[ ! -d ./Contig_Protein_Seqs ]]; then
    mkdir ./Contig_Protein_Seqs
fi

# Functionally Annotate Contigs
#for i in *_contigs.fasta;
#do
#    file=$(basename $i)
#    SAMPLE=${file%_contigs*}
#
#    # Run Prodigal
#
#    if [[ ! -f ./${SAMPLE}_proteins.faa ]]; then
#
#        prodigal -i ${file} -o ${SAMPLE}_genes.gff -f gff -a ${SAMPLE}_proteins.faa -d ${SAMPLE}_genes.fasta -s ${SAMPLE}_pot.genes.txt -p meta
#
#        cp ${SAMPLE}_genes.gff ${SAMPLE}_proteins.faa ${SAMPLE}_genes.fasta ${SAMPLE}_pot.genes.txt ./Prodigal_Results_Contigs/
#
#        cp ${SAMPLE}_proteins.faa ./Contig_Protein_Seqs
#
#    fi
#
#done

# Functionally Annotate Contigs
if [[ ! -f ./SSD_Contigs_proteins.faa ]]; then

    prodigal -i SSD_Contigs_CleanName.fa -o SSD_Contigs_genes.gff -f gff -a SSD_Contigs_proteins.faa -d SSD_Contigs_genes.fasta -s SSD_Contigs_pot.genes.txt -p meta
    
    cp SSD_Contigs_genes.gff SSD_Contigs_proteins.faa SSD_Contigs_genes.fasta SSD_Contigs_pot.genes.txt ./Prodigal_Results_Contigs/
        
    cp SSD_Contigs_proteins.faa ./Contig_Protein_Seqs
    
fi

# Functionally Annotate High Quality MAG Bins
if [[ ! -d ./Prodigal_Results_Bins ]]; then
    mkdir ./Prodigal_Results_Bins
fi

if [[ ! -d ./Bin_Protein_Seqs ]]; then
    mkdir ./Bin_Protein_Seqs
fi


## code below will only annotate bins that were previously filed into directories with name "${SAMPLE}_good_bins"
for dir in *_good_bins;
do
    file=$(basename $dir)
    SAMPLE=${file%_good*}
    
    for bin in ${dir}/*.fa;
    do
        
        bin_name=$(basename $bin)
        name=${bin_name%.fa}
        
        if [[ ! -f ./${name}_proteins.faa ]]; then
            # Run Prodigal
        
            prodigal -i ${bin} -o ${name}_genes.gff -f gff -a ${name}_proteins.faa -d ${name}_genes.fasta -s ${name}_pot.genes.txt -p meta
    
            cp ${name}_genes.gff ${name}_proteins.faa ${name}_genes.fasta ${name}_pot.genes.txt ./Prodigal_Results_Bins/
            
            cp ${name}_proteins.faa ./Bin_Protein_Seqs
        fi
        
    done

done

##https://github.com/hyattpd/Prodigal/wiki/Gene-Prediction-Modes#anonymous-mode
# Note: If running on assembled contigs, use anonymous mode; if running on MAG bins individually, use normal mode
# More on the modes here https://github.com/hyattpd/Prodigal/wiki/gene-prediction-modes
#
# Understanding the output: https://github.com/hyattpd/Prodigal/wiki/Understanding-the-Prodigal-Output
# The GFF3 format requires an "id" in the first field. Prodigal pulls the first word from the FASTA header and uses that as its id. This id is not guaranteed to be unique (the first word of various headers in the file could be identical), so we recommend the user rely on the "ID" field in the semicolon-delimited string instead.

# The fields in the semicolon-delimited string are as follows:

# ID: A unique identifier for each gene, consisting of the ordinal ID of the sequence and an ordinal ID of that gene within the sequence (separated by an underscore). For example, "4_1023" indicates the 1023rd gene in the 4th sequence in the file.
# partial: An indicator of if a gene runs off the edge of a sequence or into a gap. A "0" indicates the gene has a true boundary (a start or a stop), whereas a "1" indicates the gene is "unfinished" at that edge (i.e. a partial gene). For example, "01" means a gene is partial at the right boundary, "11" indicates both edges are incomplete, and "00" indicates a complete gene with a start and stop codon.
# start_type: The sequence of the start codon (usually ATG, GTG, or TTG). If the gene has no start codon, this field will be labeled "Edge".
# stop_type: The sequence of the stop codon (usually TAA, TGA, or TAG). If the gene has no stop codon, this field will be labeled "Edge".
# rbs_motif: The RBS motif found by Prodigal (e.g. "AGGA" or "GGA", etc.)
# rbs_spacer: The number of bases between the start codon and the observed motif.
# gc_cont: The GC content of the gene sequence.
# gc_skew: The GC skew of the gene sequence.
# conf: A confidence score for this gene, representing the probability that this gene is real, i.e. 78.3% means Prodigal believes that gene is real 78.3% of the time and a false positive 21.7% of the time.
# score: The total score for this gene.
# cscore: The hexamer coding portion of the score, i.e. how much this gene looks like a true protein.
# sscore: A score for the translation initiation site for this gene; it is the sum of the following three fields.
# rscore: A score for the RBS motif of this gene.
# uscore: A score for the sequence surrounding the start codon.
# tscore: A score for the start codon type (ATG vs. GTG vs. TTG vs. Nonstandard).
# mscore: A score for the remaining signals (stop codon type and leading/lagging strand information).

