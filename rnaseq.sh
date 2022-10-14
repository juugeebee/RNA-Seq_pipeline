#!/usr/bin/env bash

source ~/miniconda3/etc/profile.d/conda.sh


echo ""
echo "rnaseq.sh start"
echo ""


cd Fastq_cat
bash ~/SCRIPTS/RNA-Seq/agent_trimmer.sh
cd ..

###############################################################################
#### ALIGNEMENT
### STAR (Spliced Transcripts Alignment to a Reference)


conda deactivate
conda activate rnaseq

genome_dir='/media/jbogoin/Data1/References/RNA-seq/STAR'
ref='/media/jbogoin/Data1/References/fa_hg19/rna-seq/GRCh37.primary_assembly.genome.fa'
gtf_file='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.annotation.gtf'
refflat='/media/jbogoin/Data1/References/RNA-seq/refFlat_hg19.txt'

gtf_gene='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.genes.gtf'
# Obtenu en utilisant le script collapse_annotation.py sur gtf_annotation

ng_target='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions.bed'
ng_target_il='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions.interval_list'

## GENERATING GENOME INDEXES

STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $genome_dir\
 --genomeFastaFiles $ref --sjdbGTFfile $gtf_file --sjdbOverhang 71


## RUNNING MAPPING JOB
cd Fastq_trimmed

for R1 in *_R1.fastq.gz; 
    do R2=${R1/_R1_/_R2_}; 
    SAMPLE=${R1%%_*}; 
    FLOWCELL="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 3)"; 
    DEVICE="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 1 | cut -d "@" -f 2)"; 
    BARCODE="$(zcat $R1 | head -1 | awk '{print $2}' | cut -d ":" -f 4)"; 
    STAR --runThreadN 16 --genomeDir $genome_dir \
        --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
        --readFilesCommand zcat --readFilesIn $R1 $R2 \
        --outSAMattrRGline ID:${DEVICE}.${FLOWCELL}.${SAMPLE} PL:ILLUMINA PU:${FLOWCELL}.${BARCODE} LB:Illumina-StrandedmRNA-PrepLigation_${SAMPLE}_${BARCODE} SM:${SAMPLE} \
        --outFileNamePrefix ${SAMPLE} \
        --quantMode TranscriptomeSAM GeneCounts \
        --twopassMode Basic; 
done

for i in *Aligned.sortedByCoord.out.bam; do samtools index -@ 16 $i; done


mkdir -p ../BAM
mv !(*.gz|*.properties) ../BAM


## QC


for i in *Aligned.sortedByCoord.out.bam; 
    do sample=${i/Aligned.sortedByCoord.out.bam/}; 
    rnaseqc $gtf_gene $i ${sample}_RNA-SeQC --sample=${sample} --stranded=rf; 
done


conda deactivate
conda activate gatk4

for i in *Aligned.sortedByCoord.out.bam; 
    do sample=${i%Aligned.sortedByCoord.out.bam}; 
    gatk CollectRnaSeqMetrics -I $i -O ${sample}.RNAseqMetrics.txt \
    -REF_FLAT $refflat \
    -STRAND SECOND_READ_TRANSCRIPTION_STRAND; 
done


for i in *Aligned.sortedByCoord.out.bam; 
    do sample=${i%Aligned.sortedByCoord.out.bam}; 
    gatk CollectRnaSeqMetrics -I $i -O ${sample}.RNAseqMetrics.txt \
    -REF_FLAT $refflat \
    -STRAND SECOND_READ_TRANSCRIPTION_STRAND; 
done

conda deactivate
conda activate rnaseq

#Couverture totale a chaque position du bed
for i in *Aligned.sortedByCoord.out.bam; 
    do sample=${i%Aligned.sortedByCoord.out.bam}; 
    gatk CollectHsMetrics \
    -I $i \
    -O ${sample}.hsMetrics.txt \
    -R $ref \
    --BAIT_INTERVALS $ng_target_il \
    --TARGET_INTERVALS $ng_target_il;
done


mkdir ../QC
mv *.RNAseqMetrics.txt ../QC
mv *_RNA-SeQC ../QC
mv *.hsMetrics.txt ../QC
cd ../QC


multiqc .
conda deactivate


echo ""
echo "rnaseq.sh job done!"
echo ""