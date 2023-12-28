#!/usr/bin/sudo bash

source ~/miniconda3/etc/profile.d/conda.sh

echo ""
echo "rnaseq_total.sh start"
echo ""


###
mv Reports Reports qc/reports_ICM
rm -rf qc/multiqc*
###


genome_dir='/media/jbogoin/Data1/References/RNA-seq/hg38/STAR'

ref='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/GRCh38.v43.primary_assembly.genome.fa'
refflat='/media/jbogoin/Data1/References/RNA-seq/hg38/refFlat_hg38.txt'

gtf_file='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.primary_assembly.basic.annotation.gtf'

gtf_gene='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.primary_assembly.basic.annotation.genes.gtf'

gtf_transcript='/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v43.primary_assembly.basic.transcript.gtf'


# ALIGNEMENT
# STAR (Spliced Transcripts Alignment to a Reference)

echo ""
echo "ALIGNEMENT"
echo ""


conda activate rnaseq


cd fastq_merged_lanes


#***********************************************************************#
#RUNNING MAPPING JOB
for R1 in *_R1.fastq.gz; 
do R2=${R1/_R1/_R2}; 
   SAMPLE=${R1%%_*}; 
   FLOWCELL="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 3)"; 
   DEVICE="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 1 | cut -d "@" -f 2)"; 
   BARCODE="$(zcat $R1 | head -1 | awk '{print $2}' | cut -d ":" -f 4)"; 
   STAR --runThreadN 12 --genomeDir $genome_dir -n 10000 \
      --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
      --readFilesCommand zcat --readFilesIn $R1 $R2 \
      --outSAMattrRGline ID:${DEVICE}.${FLOWCELL}.${SAMPLE} PL:ILLUMINA PU:${FLOWCELL}.${BARCODE} LB:SureSelect-XT-HS2mRNA-Library_${SAMPLE}_${BARCODE} SM:${SAMPLE} \
      --outFileNamePrefix ${SAMPLE} --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic; 
done


#***********************************************************************#
# CREATION DES INDEXS
for i in *Aligned.sortedByCoord.out.bam; do samtools index -@ 24 $i; done


############# QC #############

echo "QC"
echo ""


***********************************************************************#
echo "RNASEQC"
echo ""

for i in *Aligned.sortedByCoord.out.bam; 
   do sample=${i/Aligned.sortedByCoord.out.bam/}; 
   rnaseqc $gtf_gene $i ${sample}_RNA-SeQC --sample=${sample} --stranded='rf'; 
done


conda deactivate

conda activate gatk4


#***********************************************************************#
echo "CollectRnaSeqMetrics"
echo ""

for i in *Aligned.sortedByCoord.out.bam; 
   do sample=${i/Aligned.sortedByCoord.out.bam/}; 
   gatk CollectRnaSeqMetrics -I $i -O ${sample}.RNAseqMetrics.txt \
   --REF_FLAT $refflat \
   -R $ref -STRAND SECOND_READ_TRANSCRIPTION_STRAND \
   --RIBOSOMAL_INTERVALS /media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.rRNA.transcripts.interval_list; 
done

conda deactivate


#***********************************************************************#
echo "salmon"
echo ""

conda activate salmon


# INDEX
#salmon index -t '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts.fa' \
#-i '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts-salmon.idx'



# COUNT
for R1 in *_R1.fastq.gz; 
   do R2=${R1/_R1/_R2};
   sample=${R1/_S**_R1_001.fastq.gz/};
   salmon quant -i '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts-salmon.idx' \
   -l ISR \
   -1 $R1 -2 $R2 \
   --validateMappings \
   -p 24 \
   -o ../qc/salmon/$sample;
done


conda deactivate 

mkdir -p ../BAM
mv `ls . | grep -v "\.gz$"` ../BAM


#***********************************************************************#
echo "fraser"
echo ""

bash ~/SCRIPTS/RNA-Seq/DROP/drop.sh



### CLEANING


mv *_RNA-SeQC ../qc
mv *.RNAseqMetrics.txt ../qc
mv drop ../qc

cd ../qc

mkdir -p RNA-SeQC
mv *RNA-SeQC RNA-SeQC
mkdir -p RnaSeqMetrics
mv *.RNAseqMetrics.txt RnaSeqMetrics


#***********************************************************************#
echo "multiqc"
echo ""

conda activate rnaseq

multiqc -f .

conda deactivate


echo ""
echo "rnaseq_total.sh job done!"
echo ""
