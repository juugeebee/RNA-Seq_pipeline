#!/usr/bin/env bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

echo ""
echo "rnaseq_with_umi.sh start"
echo ""


genome_dir='/media/jbogoin/Data1/References/RNA-seq/STAR'
ref='/media/jbogoin/Data1/References/fa_hg19/rna-seq/GRCh37.primary_assembly.genome.fa'
gtf_file='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.annotation.gtf'
refflat='/media/jbogoin/Data1/References/RNA-seq/refFlat_hg19.txt'

gtf_gene='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.genes.gtf'
# Obtenu en utilisant le script collapse_annotation.py sur gtf_annotation

ng_target='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions.bed'
ng_target_il='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions.interval_list'


# TRIMMING DES ADATATEURS

cd Fastq_cat
echo "TRIMMER"
echo ""
bash ~/SCRIPTS/RNA-Seq/agent_trimmer.sh
cd ..

conda deactivate


## ALIGNEMENT
## BWA-MEM

echo "1er Alignement avec BWA-mem"
echo ""

conda activate fastq_bam_env

cd Fastq_trimmed

gzip -d *.txt.gz


for R1 in *_R1.fastq.gz; 
   do R2=${R1/_R1/_R2}; \
   SAMPLE=${R1%%_*}; \
   FLOWCELL="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 3)"; \
   DEVICE="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 1 | cut -d "@" -f 2)"; \
   BARCODE="$(zcat $R1 | head -1 | awk '{print $2}' | cut -d ":" -f 4)"; \
   RG=$(echo "\"@RG\tID:${DEVICE}.${FLOWCELL}.${SAMPLE}\tPU:${FLOWCELL}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:${SAMPLE}-${BARCODE}\""); \
   MAPPING_CMD=$(echo "bwa mem -M -Y -t 24 -R $RG $ref $R1 $R2 -C | samtools sort -@ 6 -o ${SAMPLE}.bam -"); \
   eval $MAPPING_CMD;
done

conda deactivate

mkdir -p ../BAM_bwa-mem
mv *.bam ../BAM_bwa-mem
cd ../BAM_bwa-mem


## CREATION DES INDEXS                                                                                                                                
for i  in *.bam; do samtools index -@ 16 $i; done


## Dedup avec traitement des UMI

echo "CReaK"
echo ""

conda activate agent_env

for i in *.bam;
    do sample=${i%.bam};
    java -jar '/home/jbogoin/AGeNT_3.0.4/agent3.0/lib/creak-1.0.5.jar' \
        $i -f -F -r --consensus-mode SINGLE --MBC-mismatch 1 \
        --bed-file $ng_target \
        --output-bam-file ${sample}.dedup.bam;
done

conda deactivate


## BAM TO FASTQ

echo "BAM to FASTQ"
echo ""

bash ~/SCRIPTS/Files_preparation/bam_to_fastq.sh

mkdir -p ../Fastq_UMI
mv *.fastq.gz ../Fastq_UMI
cd ../Fastq_UMI


## REALIGNEMENT AVEC STAR

echo "2eme Alignement avec STAR"
echo "" 

conda activate rnaseq

for R1 in *.R1.fastq.gz; 
    do R2=${R1/.R1./.R2.}; 
    SAMPLE=${R1%%.*}; 
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

mkdir -p ../BAM_STAR_avec_UMI
mv !(*.gz) ../BAM_STAR_avec_UMI
cd ../BAM_STAR_avec_UMI


## CREATION DES INDEXS

for i in *.sortedByCoord.out.bam; do samtools index -@ 16 $i; done


## QC

echo "QC"
echo ""

echo "CollectHsMetrics"
echo ""

conda activate gatk4

# Couverture totale a chaque position du bed
for i in *.sortedByCoord.out.bam; 
    do sample=${i%.sortedByCoord.out.bam}; 
    gatk CollectHsMetrics \
    -I $i \
    -O ${sample}.hsMetrics.txt \
    -R $ref \
    --BAIT_INTERVALS $ng_target_il \
    --TARGET_INTERVALS $ng_target_il \
    --PER_TARGET_COVERAGE ${sample}.hsMetrics_pertargetcoverage.txt;
done

conda deactivate
conda activate rnaseq

echo "RNASEQC"
echo ""

for i in *.sortedByCoord.out.bam; 
    do sample=${i%.sortedByCoord.out.bam}; 
    rnaseqc $gtf_gene $i ${sample}_RNA-SeQC --sample=${sample} --stranded=rf; 
done


mkdir ../QC_avec_UMI
mv *_RNA-SeQC ../QC_avec_UMI
mv *.hsMetrics.txt ../QC_avec_UMI
mv *_pertargetcoverage.txt ../QC_avec_UMI
cd ../QC_avec_UMI



conda deactivate
conda activate rnaseq

echo "MULTIQC"

multiqc -f .

conda deactivate


mkdir -p hsMetrics
mv *.hsMetrics.txt hsMetrics/
mkdir -p pertargetcoverage
mv *_pertargetcoverage.txt pertargetcoverage/
mkdir -p RNA-SeQC
mv *RNA-SeQC RNA-SeQC/


conda deactivate


echo ""
echo "rnaseq_with_umi.sh job done!"
echo ""
