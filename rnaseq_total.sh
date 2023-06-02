#!/usr/bin/sudo bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq

echo ""
echo "rnaseq_total.sh start"
echo ""


genome_dir='/media/jbogoin/Data1/References/RNA-seq/hg38/STAR'
ref='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/GRCh38.v43.primary_assembly.genome.fa'
gtf_file='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.primary_assembly.basic.annotation.gtf'
refflat='/media/jbogoin/Data1/References/RNA-seq/hg38/refFlat_hg38.txt'


gtf_gene='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.primary_assembly.basic.annotation.genes.gtf'
# Obtenu en utilisant le script collapse_annotation.py sur gtf_annotation

glob='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/globines.bed'
glob_il='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/goblines.interval_list'


# ALIGNEMENT
# STAR (Spliced Transcripts Alignment to a Reference)

echo "ALIGNEMENT"
echo ""

# GENERATING GENOME INDEXES
STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $ref --sjdbGTFfile $gtf_file --sjdbOverhang 72


# RUNNING MAPPING JOB
cd Fastq

for R1 in *_R1_001.fastq.gz; 
do R2=${R1/_R1/_R2}; 
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


mkdir -p ../BAM
mv !(*.gz) ../BAM
mv `ls . | grep -v "\.gz$"` ../BAM
cd ../BAM


# CREATION DES INDEXS

for i in *Aligned.sortedByCoord.out.bam; do samtools index -@ 24 $i; done


## QC

echo "QC"
echo ""


conda activate rnaseq

echo "RNASEQC"
echo ""

for i in *Aligned.sortedByCoord.out.bam; 
   do sample=${i/Aligned.sortedByCoord.out.bam/}; 
   rnaseqc $gtf_gene $i ${sample}_RNA-SeQC --sample=${sample} --stranded='rf'; 
done


conda deactivate

conda activate gatk4


echo "CollectRnaSeqMetrics"
echo ""


for i in *Aligned.sortedByCoord.out.bam; 
   do sample=${i/Aligned.sortedByCoord.out.bam/}; 
   gatk CollectRnaSeqMetrics -I $i -O ${sample}.RNAseqMetrics.txt \
   --REF_FLAT $refflat \
   -R $ref -STRAND SECOND_READ_TRANSCRIPTION_STRAND \
   --RIBOSOMAL_INTERVALS /media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.rRNA.transcripts.interval_list; 
done


echo "pertargetcoverage"
echo ""


# Couverture totale a chaque position du bed
for i in *Aligned.sortedByCoord.out.bam; 
    do sample=${i/Aligned.sortedByCoord.out.bam/}; 
    gatk CollectHsMetrics \
    -I $i \
    -O ${sample}.hsMetrics.txt \
    -R $ref \
    --BAIT_INTERVALS $glob_il \
    --TARGET_INTERVALS $glob_il \
    --PER_TARGET_COVERAGE ${sample}.hsMetrics_pertargetcoverage.txt;
done


mkdir ../QC
mv *_RNA-SeQC ../QC
mv *.RNAseqMetrics.txt ../QC
mv *.hsMetrics.txt ../QC
mv *_pertargetcoverage.txt ../QC
cd ../QC


conda deactivate


conda activate rnaseq

echo "multiqc"
echo ""

multiqc -f .

conda deactivate


mkdir -p hsMetrics
mv *.hsMetrics.txt hsMetrics
mkdir -p pertargetcoverage
mv *_pertargetcoverage.txt pertargetcoverage
mkdir -p RNA-SeQC
mv *RNA-SeQC RNA-SeQC
mkdir -p RnaSeqMetrics
mv *.RNAseqMetrics.txt RnaSeqMetrics


echo ""
echo "rnaseq_total.sh job done!"
echo ""
