#!/usr/bin/sudo bash

source ~/miniconda3/etc/profile.d/conda.sh

echo ""
echo "rnaseq_total.sh start"
echo ""


#########################################
## A LANCER DANS LE DOSSIER RACINE DU RUN
#########################################


genome_dir='/media/jbogoin/Data1/References/RNA-seq/hg38/STAR'

ref='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/GRCh38.v43.primary_assembly.genome.fa'
refflat='/media/jbogoin/Data1/References/RNA-seq/hg38/refFlat_hg38.txt'

gtf_file='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.primary_assembly.basic.annotation.gtf'

gtf_gene='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.primary_assembly.basic.annotation.genes.gtf'
# Obtenu en utilisant le script collapse_annotation.py sur gtf_annotation

gtf_transcript='/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v43.primary_assembly.basic.transcript.gtf'


#***********************************************************************
#echo "FastQC"
#echo ""

cd Fastq

#conda activate fastqc

#mkdir -p ../QC/fastqc
#for R1 in *_R1_001.fastq.gz; do R2=${R1/_R1/_R2}; fastqc -o ../QC/fastqc -f fastq $R1 $R2; done

#conda deactivate


#***********************************************************************#
# ALIGNEMENT
# STAR (Spliced Transcripts Alignment to a Reference)


echo ""
echo "ALIGNEMENT"
echo ""


conda activate rnaseq


# ***********************************************************************#
# GENERATING GENOME INDEXES
# STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $ref --sjdbGTFfile $gtf_file --sjdbOverhang 72


#***********************************************************************#
#RUNNING MAPPING JOB
for R1 in *_R1_001.fastq.gz; 
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


#***********************************************************************#
echo "RNASEQC"
echo ""

for i in *Aligned.sortedByCoord.out.bam; 
   do sample=${i/Aligned.sortedByCoord.out.bam/}; 
   rnaseqc $gtf_gene $i ${sample}_RNA-SeQC --sample=${sample} --stranded='rf'; 
done


conda deactivate


#***********************************************************************#
echo "CollectRnaSeqMetrics"
echo ""


conda activate gatk4


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
for R1 in *_R1_001.fastq.gz; 
   do R2=${R1/_R1/_R2};
   sample=${R1/_S**_R1_001.fastq.gz/};
   salmon quant -i '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts-salmon.idx' \
   -l ISR \
   -1 $R1 -2 $R2 \
   --validateMappings \
   -p 24 \
   -o ../QC/salmon/$sample;
done


conda deactivate 


# #***********************************************************************#
# # echo "FeatureCounts"
# # echo ""

# # conda activate FeatureCounts

# # for i in *Aligned.sortedByCoord.out.bam;
# #    do sample=${i/Aligned.sortedByCoord.out.bam/};
# #    featureCounts -p -O -T 24 -s 2 \
# #    -t transcript \
# #    -a $gtf_transcript \
# #    -o ${sample}_featureCounts_output.txt \
# #    $i;
# # done

# # conda deactivate 


# #***********************************************************************#
# # echo "htseq-count"
# # echo ""


# # conda activate htseq

# # for i in *Aligned.sortedByCoord.out.bam;
# #    do sample=${i/Aligned.sortedByCoord.out.bam/}; 
# #    htseq-count -f bam -s reverse -t transcript \
# #    --secondary-alignments ignore --supplementary-alignments ignore \
# #    -c ../QC/htseq/${sample}_htseq_output.tsv \
# #    -p bam \
# #    -n 24 \
# #    $i \
# #    $gtf_transcript;
# # done

# # conda deactivate


# #***********************************************************************#
# # echo "rsem"
# # echo ""

# # conda activate rsem


# # # INDEX
# # #mkdir -p /media/jbogoin/Data1/References/RNA-seq/hg38/ref/human_gencode

# # #rsem-prepare-reference --gtf '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v43.primary_assembly.basic.annotation.gtf' \
# # #--star --star-path '/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/GRCh38.v43.primary_assembly.genome.fa' \
# # #/media/jbogoin/Data1/References/RNA-seq/hg38/ref/human_gencode


# # # COUNT
# # for R1 in *_R1_001.fastq.gz; 
# #    do R2=${R1/_R1/_R2};
# #    sample=${R1/_S**_R1_001.fastq.gz/};
# #    rsem-calculate-expression \
# #    --paired-end -p 24 --append-names --star --star-gzipped-read-file --no-bam-output\
# #    $R1 $R2 /media/jbogoin/Data1/References/RNA-seq/hg38/ref/human_gencode $sample;
# # done
 
# # conda deactivate


# #***********************************************************************
# # echo "kallisto"
# # echo ""

# # conda activate kallisto


# # # INDEX
# # #kallisto index '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts.fa' \
# # #-i '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts-kallisto.idx'


# # # COUNT
# # for R1 in *_R1_001.fastq.gz; 
# #    do R2=${R1/_R1/_R2};
# #    sample=${R1/_S**_R1_001.fastq.gz/};
# #    kallisto quant\
# #    -i '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts-kallisto.idx'\
# #    -o ../QC/kallisto/$sample \
# #    --rf \
# #    -t 24 \
# #    --genomebam --gtf $gtf_transcript \
# #    $R1 $R2 ;
# # done

# # conda deactivate


#***********************************************************************
### CLEANING


mv Logs ../QC
mv Reports ../QC
mv *_RNA-SeQC ../QC
mv *.RNAseqMetrics.txt ../QC
# mv *.hsMetrics.txt ../QC
# mv *_pertargetcoverage.txt ../QC
# mv *_featureCounts* ../QC
# mv *.stat ../QC
# mv *.results ../QC


cd ../QC


#***********************************************************************#
echo "multiqc"
echo ""

conda activate rnaseq

multiqc -f .

conda deactivate


# mkdir -p hsMetrics
# mv *.hsMetrics.txt hsMetrics
# mkdir -p pertargetcoverage
# mv *_pertargetcoverage.txt pertargetcoverage
mkdir -p RNA-SeQC
mv *RNA-SeQC RNA-SeQC
mkdir -p RnaSeqMetrics
mv *.RNAseqMetrics.txt RnaSeqMetrics
# mkdir -p FeatureCounts
# mv *_featureCounts* FeatureCounts
# mkdir -p rsem
# mv *.stat rsem
# mv *.results rsem 


mkdir -p ../BAM

cd ../Fastq
# mv !(*.gz) ../BAM
mv `ls . | grep -v "\.gz$"` ../BAM

cd ..


#***********************************************************************#
echo "fraser"
echo ""

bash ~/SCRIPTS/RNA-Seq/DROP/drop.sh


echo ""
echo "rnaseq_total.sh job done!"
echo ""
