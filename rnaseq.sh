#!/usr/bin/env bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

echo ""
echo "rnaseq.sh start"
echo ""


genome_dir='/media/jbogoin/Data1/References/RNA-seq/STAR'
ref='/media/jbogoin/Data1/References/fa_hg19/rna-seq/GRCh37.primary_assembly.genome.fa'
gtf_file='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.annotation.gtf'
refflat='/media/jbogoin/Data1/References/RNA-seq/refFlat_hg19.txt'


gtf_gene='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.genes.gtf'
# Obtenu en utilisant le script collapse_annotation.py sur gtf_annotation

# ng_target='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions.bed'
# ng_target_il='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions.interval_list'


ng_target='/media/jbogoin/Data2/Donnees_brutes/hg19/NG_TRS/bed/CEREBMDv1_hg19_13Mar2019_primary_targets.bed'
ng_target_il='/media/jbogoin/Data2/Donnees_brutes/hg19/NG_TRS/bed/CEREBMDv1_hg19_13Mar2019_primary_targets.interval_list'


# TRIMMING DES ADATATEURS

cd Fastq_cat
echo "TRIMMER"
echo ""
bash ~/SCRIPTS/RNA-Seq/agent_trimmer.sh
cd ..


## ALIGNEMENT
## STAR (Spliced Transcripts Alignment to a Reference)


echo "ALIGNEMENT"
echo ""


## GENERATING GENOME INDEXES
STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $genome_dir\
 --genomeFastaFiles $ref --sjdbGTFfile $gtf_file --sjdbOverhang 72


## RUNNING MAPPING JOB
cd Fastq_trimmed

for R1 in *_R1.fastq.gz; 
    do R2=${R1/_R1./_R2.}; 
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

mkdir -p ../BAM_STAR_sans_UMI
mv !(*.gz) ../BAM_STAR_sans_UMI
cd ../BAM_STAR_sans_UMI
mv *.properties ../Fastq_trimmed


## MarkDuplicates sans traitement des UMI

echo "MarkDuplicates"
echo ""

conda activate gatk4

for i in *Aligned.sortedByCoord.out.bam;
    do sample=${i%Aligned.sortedByCoord.out.bam};
    gatk MarkDuplicates \
        -I $i \
        -O ${sample}.marked_duplicates.bam \
        -M ${sample}.marked_dup_metrics.txt;
done

conda deactivate


## CREATION DES INDEXS

for i in *.marked_duplicates.bam; do samtools index -@ 16 $i; done


## QC
echo "RNASEQC"
echo ""

conda activate rnaseq

for i in *.marked_duplicates.bam; 
    do sample=${i%.marked_duplicates.bam}; 
    rnaseqc $gtf_gene $i ${sample}_RNA-SeQC --sample=${sample} --stranded=rf; 
done


conda deactivate
conda activate gatk4


echo "CollectHsMetrics"
echo ""


# Couverture totale a chaque position du bed
for i in *.marked_duplicates.bam; 
    do sample=${i%.marked_duplicates.bam}; 
    gatk CollectHsMetrics \
    -I $i \
    -O ${sample}.hsMetrics.txt \
    -R $ref \
    --BAIT_INTERVALS $ng_target_il \
    --TARGET_INTERVALS $ng_target_il \
    --PER_TARGET_COVERAGE ${sample}.hsMetrics_pertargetcoverage.txt;
done


mkdir ../QC_sans_UMI
mv *_RNA-SeQC ../QC_sans_UMI
mv *.hsMetrics.txt ../QC_sans_UMI
mv *_pertargetcoverage.txt ../QC_sans_UMI
cd ../QC_sans_UMI


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
echo "rnaseq.sh job done!"
echo ""
