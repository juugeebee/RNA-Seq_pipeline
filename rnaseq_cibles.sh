#!/usr/bin/env bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

echo ""
echo "rnaseq_cibles.sh start"
echo ""


genome_dir='/media/jbogoin/Data1/References/RNA-seq/hg19/STAR'
ref='/media/jbogoin/Data1/References/fa_hg19/rna-seq/GRCh37.primary_assembly.genome.fa'
gtf_file='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.annotation.gtf'
refflat='/media/jbogoin/Data1/References/RNA-seq/hg19/refFlat_hg19.txt'


gtf_gene='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.genes.gtf'
# Obtenu en utilisant le script collapse_annotation.py sur gtf_annotation


ng_target='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions_hg19.bed'
ng_target_il='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions_hg19.interval_list'


# TRIMMING DES ADATATEURS

# cd Fastq_cat
# echo "TRIMMER"
# echo ""
# bash ~/SCRIPTS/RNA-Seq/agent_trimmer.sh
# cd ..


# conda activate rnaseq


# ## ALIGNEMENT
# ## STAR (Spliced Transcripts Alignment to a Reference)


# echo "ALIGNEMENT"
# echo ""


# ## GENERATING GENOME INDEXES
# STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $genome_dir\
#  --genomeFastaFiles $ref --sjdbGTFfile $gtf_file --sjdbOverhang 72


# ## RUNNING MAPPING JOB
# cd Fastq_trimmed

for R1 in *_R1_001.fastq.gz; 
    do R2=${R1/_R1_/_R2_}; 
    SAMPLE=${R1%%_*}; 
    FLOWCELL="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 3)"; 
    DEVICE="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 1 | cut -d "@" -f 2)"; 
    BARCODE="$(zcat $R1 | head -1 | awk '{print $2}' | cut -d ":" -f 4)"; 
    STAR --runThreadN 16 --genomeDir $genome_dir \
        --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
        --readFilesCommand zcat --readFilesIn $R1 $R2 \
        --outSAMattrRGline ID:${DEVICE}.${FLOWCELL}.${SAMPLE} PL:ILLUMINA PU:${FLOWCELL}.${BARCODE} LB:SureSelect-XT-HS2mRNA-Library_${SAMPLE}_${BARCODE} SM:${SAMPLE} \
        --outFileNamePrefix ${SAMPLE} \
        --quantMode TranscriptomeSAM GeneCounts \
        --twopassMode Basic; 
done

conda deactivate


## CREATION DES INDEXS

for i in *Aligned.sortedByCoord.out.bam; do samtools index -@ 16 $i; done


## QC


conda activate gatk4


echo "CollectHsMetrics"
echo ""


# Couverture totale a chaque position du bed
for i in *Aligned.sortedByCoord.out.bam; 
    do sample=${i%Aligned.sortedByCoord.out.bam}; 
    gatk CollectHsMetrics \
    -I $i \
    -O ${sample}.hsMetrics.txt \
    -R $ref \
    --BAIT_INTERVALS $ng_target_il \
    --TARGET_INTERVALS $ng_target_il \
    --PER_TARGET_COVERAGE ${sample}.hsMetrics_pertargetcoverage.txt;
done


mkdir ../QC_hg19
mv *.hsMetrics.txt ../QC_hg19
mv *_pertargetcoverage.txt ../QC_hg19
cd ../QC_hg19


conda deactivate
# conda activate rnaseq

# echo "MULTIQC"

# multiqc -f .

# conda deactivate


mkdir -p hsMetrics
mv *.hsMetrics.txt hsMetrics/
mkdir -p pertargetcoverage
mv *_pertargetcoverage.txt pertargetcoverage/
# mkdir -p RNA-SeQC
# mv *RNA-SeQC RNA-SeQC/


mkdir -p ../BAM_hg19
mv `ls . | grep -v "\.gz$"` ../BAM_hg19


echo ""
echo "rnaseq_cilbes.sh job done!"
echo ""
