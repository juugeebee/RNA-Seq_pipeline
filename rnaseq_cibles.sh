#!/usr/bin/env bash

source ~/miniconda3/etc/profile.d/conda.sh


echo ""
echo "rnaseq_cibles.sh start"
echo ""


genome_dir='/media/jbogoin/Data1/References/RNA-seq/hg19/STAR'
ref='/media/jbogoin/Data1/References/fa_hg19/rna-seq/GRCh37.primary_assembly.genome.fa'
gtf_file='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.annotation.gtf'
refflat='/media/jbogoin/Data1/References/RNA-seq/hg19/refFlat_hg19.txt'


gtf_gene='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.genes.gtf'
# Obtenu en utilisant le script collapse_annotation.py sur gtf_annotation

# NG
# target='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions_hg19.bed'
# target_il='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions_hg19.interval_list'

# OA
# target='/media/jbogoin/Data1/References/cibles_panel_OA/CAPONCOV3_design1.bed'
# target_il='/media/jbogoin/Data1/References/cibles_panel_OA/CAPONCOV3_design1.interval_list'

#globines
target='/media/jbogoin/Data1/References/fa_hg19/rna-seq/globines_hg19.bed'
target_il='/media/jbogoin/Data1/References/fa_hg19/rna-seq/globines_hg19.interval_list'


#***********************************************************************#
# TRIMMING DES ADATATEURS

# echo "TRIMMER"
# echo ""

cd Fastq
# conda activate rnaseq
# bash ~/SCRIPTS/RNA-Seq/agent_trimmer.sh
# cd ..


#***********************************************************************
echo "FastQC"
echo ""

conda activate fastqc

mkdir -p ../QC/fastqc
for R1 in *_R1_001.fastq.gz; do R2=${R1/_R1/_R2}; fastqc -o ../QC/fastqc -f fastq $R1 $R2; done

conda deactivate


#***********************************************************************#
# ALIGNEMENT
# STAR (Spliced Transcripts Alignment to a Reference)


echo ""
echo "ALIGNEMENT"
echo ""


conda activate rnaseq


# ***********************************************************************#
# GENERATING GENOME INDEXES
# STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $genome_dir\
#   --genomeFastaFiles $ref --sjdbGTFfile $gtf_file --sjdbOverhang 72


#***********************************************************************#
# RUNNING MAPPING JOB
# cd Fastq_trimmed

for R1 in *_R1_001.fastq.gz; 
    do R2=${R1/_R1_/_R2_}; 
    SAMPLE=${R1%%_*}; 
    FLOWCELL="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 3)"; 
    DEVICE="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 1 | cut -d "@" -f 2)"; 
    BARCODE="$(zcat $R1 | head -1 | awk '{print $2}' | cut -d ":" -f 4)"; 
    STAR --runThreadN 12 --genomeDir $genome_dir \
        --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
        --readFilesCommand zcat --readFilesIn $R1 $R2 \
        --outSAMattrRGline ID:${DEVICE}.${FLOWCELL}.${SAMPLE} PL:ILLUMINA PU:${FLOWCELL}.${BARCODE} LB:SureSelect-XT-HS2mRNA-Library_${SAMPLE}_${BARCODE} SM:${SAMPLE} \
        --outFileNamePrefix ${SAMPLE} \
        --quantMode TranscriptomeSAM GeneCounts \
        --twopassMode Basic; 
done


## CREATION DES INDEXS
for i in *Aligned.sortedByCoord.out.bam; do samtools index -@ 16 $i; done


############# QC #############

echo "QC"
echo ""



#***********************************************************************#
echo "RNASEQC"
echo ""


for i in *Aligned.sortedByCoord.out.bam; 
   do sample=${i/Aligned.sortedByCoord.out.bam/}; 
   rnaseqc $gtf_gene $i ${sample}_RNA-SeQC --sample=${sample} --stranded='rf'; 
done

conda deactivate



#***********************************************************************#
echo "CollectHsMetrics"
echo ""

conda activate gatk4


#***********************************************************************#
# Couverture totale a chaque position du bed
echo "pertargetcoverage"
echo ""



for i in *Aligned.sortedByCoord.out.bam; 
    do sample=${i%Aligned.sortedByCoord.out.bam}; 
    gatk CollectHsMetrics \
    -I $i \
    -O ${sample}.hsMetrics.txt \
    -R $ref \
    --BAIT_INTERVALS $target_il \
    --TARGET_INTERVALS $target_il \
    --PER_TARGET_COVERAGE ${sample}.pertargetcoverage.txt;
done


#***********************************************************************#
echo "CollectRnaSeqMetrics"
echo ""

for i in *Aligned.sortedByCoord.out.bam; 
   do sample=${i/Aligned.sortedByCoord.out.bam/}; 
   gatk CollectRnaSeqMetrics -I $i -O ${sample}.RNAseqMetrics.txt \
   --REF_FLAT $refflat \
   -R $ref -STRAND SECOND_READ_TRANSCRIPTION_STRAND \
   --RIBOSOMAL_INTERVALS '/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.rRNA.transcripts.interval_list'; 
done


conda deactivate


#***********************************************************************#
echo "salmon"
echo ""

conda activate salmon

mkdir ../QC

# INDEX
#salmon index -t '/media/jbogoin/Data1/References/RNA-seq/hg19/gencode.v41lift37.transcripts.fa' \
#-i '/media/jbogoin/Data1/References/RNA-seq/hg19/gencode.v41lift37.transcripts-salmon.idx'

# COUNT
for R1 in *_R1_001.fastq.gz; 
   do R2=${R1/_R1/_R2};
   sample=${R1/_R1_001.fastq.gz/};
   salmon quant -i '/media/jbogoin/Data1/References/RNA-seq/hg19/gencode.v41lift37.transcripts-salmon.idx' \
   -l ISR \
   -1 $R1 -2 $R2 \
   --validateMappings \
   -p 24 \
   -o ../QC/salmon/$sample;
done

conda deactivate 


#***********************************************************************
### CLEANING


mv *_RNA-SeQC ../QC
mv *.RNAseqMetrics.txt ../QC
mv *.hsMetrics.txt ../QC
mv *_pertargetcoverage.txt ../QC

cd ../QC



#***********************************************************************#
echo "multiqc"
echo ""

conda activate rnaseq

multiqc -f .

conda deactivate


mkdir -p hsMetrics
mv *hsMetrics.txt hsMetrics/
mkdir -p pertargetcoverage
mv *pertargetcoverage.txt pertargetcoverage/
mkdir -p RNA-SeQC
mv *RNA-SeQC RNA-SeQC/
mkdir -p RnaSeqMetrics
mv *.RNAseqMetrics.txt RnaSeqMetrics


mkdir -p ../BAM

cd ../Fastq
# cd Fastq_trimmed
mv `ls . | grep -v "\.gz$"` ../BAM

cd ..


#***********************************************************************#
echo "fraser"
echo ""

bash ~/SCRIPTS/RNA-Seq/DROP/drop_hg19.sh


echo ""
echo "rnaseq_cibles.sh job done!"
echo ""