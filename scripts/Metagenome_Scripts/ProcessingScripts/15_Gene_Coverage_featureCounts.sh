#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 # must match the # of threads (-t #)
##SBATCH --mem-per-cpu=400G # * if threading, do not let this line run (use ##). Cannot ask for too much memory per cpu!
#SBATCH --mem=400G # < if you are threading, ask for whole memory and not memory per CPU
##SBATCH --time=20:00:00     # 20 hrs
#SBATCH --output=Gene_Coverage_featureCounts_8.21.23.stdout
#SBATCH --mail-user=hfreu002@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="Gene_Coverage_featureCounts_8.21.23"
#SBATCH -p aronsonlab

module load subread/2.0.3

if [[ ! -d ./Gene_Fxn_Files/Gene_Coverage/GeneCoverage_featureCounts ]]; then
    mkdir Gene_Fxn_Files/Gene_Coverage/GeneCoverage_featureCounts
fi
outdir="/bigdata/aronsonlab/shared/SaltonSea/Metagenomes/SeqCenter_3.30.2023/Gene_Fxn_Files/Gene_Coverage/GeneCoverage_featureCounts"
## Coverage by gene position in BED file

# featureCounts on genes found in bins
for i in ./*_genes.gff;
do
    file=${i##*/}
    SAMPLE=${file%_genes*}
    echo $file
    
    if [[ ${SAMPLE} == *"bin."* ]]; then
            
        OGSample=${SAMPLE%_bin*}
        # use original bam file for entire sample to determine coverage in individual bins
        featureCounts -p -O -t CDS -g ID -o ${outdir}/${SAMPLE}_fC_gene_coverage.tsv -a ${SAMPLE}_genes.gff -T 4 ${OGSample}_contigs_sorted.bam

    fi
    
done

# featureCounts on genes found in co-assembled contigs
for i in ./*_contigs_sorted.bam;
do
    file=${i##*/}
    SAMPLE=${file%_contigs*}
    echo $SAMPLE
    
    featureCounts -p -O -t CDS -g ID -o ${outdir}/${SAMPLE}_fC_gene_coverage.tsv -a SSD_Contigs_genes.gff -T 4 ${SAMPLE}_contigs_sorted.bam
    
    
done

# featureCounts on genes found in bins and contigs (single-assembly so contigs per sample, not co-assembly)
#for i in ./*_genes.gff;
#do
#    file=${i##*/}
#    SAMPLE=${file%_genes*}
#    echo $file
#
#    if [[ ${SAMPLE} == *"bin."* ]]; then
#
#        OGSample=${SAMPLE%_bin*}
#        # use original bam file for entire sample to determine coverage in individual bins
#        featureCounts -p -O -t CDS -g ID -o ${outdir}/${SAMPLE}_fC_gene_coverage.tsv -a ${SAMPLE}_genes.gff -T 4 ${OGSample}_contigs_sorted.bam
#
#    elif [[ ${SAMPLE} != *"bin."* ]]; then
#
#        featureCounts -p -O -t CDS -g ID -o ${outdir}/${SAMPLE}_fC_gene_coverage.tsv -a ${SAMPLE}_genes.gff -T 4 ${SAMPLE}_contigs_sorted.bam
#
#    fi
#
#done



# featureCounts website: https://subread.sourceforge.net/featureCounts.html#:~:text=featureCounts%20is%20a%20highly%20efficient,and%20genomic%20DNA%2Dseq%20reads.

# tutorial link describing featureCounts: https://staff.fnwi.uva.nl/a.u.s.heintzbuschart/08_calcFunc.html#running-featurecounts-for-all-genes
