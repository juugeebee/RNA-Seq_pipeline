#!/usr/bin/sudo bash

source ~/miniconda3/etc/profile.d/conda.sh

echo ""
echo "rnaseq_total.sh start"
echo ""


#########################################
## A LANCER DANS LE DOSSIER RACINE DU RUN
#########################################

genome_dir='/media/jbogoin/Data1/References/RNA-seq/hg38/STAR_71pb_v43'

# Ref gencode v43
ref='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/GRCh38.v43.primary_assembly.genome.fa'

refflat='/media/jbogoin/Data1/References/RNA-seq/hg38/refFlat_hg38.txt'

# GTF file gencode v43
gtf_file='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.primary_assembly.basic.annotation.gtf'

# GTF gene v43
gtf_gene='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.primary_assembly.basic.annotation.genes.gtf'
# Obtenu en utilisant le script collapse_annotation.py sur gtf_annotation

# GTF transcript v43
gtf_transcript='/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v43.primary_assembly.basic.transcript.gtf'


#***********************************************************************#
# ALIGNEMENT
# STAR (Spliced Transcripts Alignment to a Reference)


echo ""
echo "ALIGNEMENT"
echo ""


conda activate rnaseq


cd Fastq


# ***********************************************************************#
# GENERATING GENOME INDEXES
# STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $ref --sjdbGTFfile $gtf_file --sjdbOverhang 72


#***********************************************************************#
#RUNNING MAPPING JOB
for R1 in *_R1*.fastq.gz; 
do R2=${R1/_R1/_R2}; 
   SAMPLE=${R1%%_*};
   FLOWCELL="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 3)"; 
   DEVICE="$(zcat $R1 | head -1 | awk '{print $1}' | cut -d ":" -f 1 | cut -d "@" -f 2)"; 
   BARCODE="$(zcat $R1 | head -1 | awk '{print $2}' | cut -d ":" -f 4)"; 
   STAR --runThreadN 24 --genomeDir $genome_dir -n 10000 \
      --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within \
      --readFilesCommand zcat --readFilesIn $R1 $R2 \
      --outSAMattrRGline ID:${DEVICE}.${FLOWCELL}.${SAMPLE} PL:ILLUMINA PU:${FLOWCELL}.${BARCODE} LB:Il-str-mRNA-D SM:${SAMPLE} \
      --outFileNamePrefix ${SAMPLE}. \
      --twopassMode Basic; 
done


#***********************************************************************#
# CREATION DES INDEXS

# time parallel -j 12 "samtools index {}" ::: *Aligned.sortedByCoord.out.bam
for i in *Aligned.sortedByCoord.out.bam; do samtools index -@ 16 $i; done


############# QC #############

echo "QC"
echo ""


***********************************************************************
echo "FastQC"
echo ""

conda activate fastqc

mkdir -p ../QC/fastqc
time parallel -j 16 fastqc {} ::: *.fastq.gz
# for R1 in *_R1_001.fastq.gz; do R2=${R1/_R1/_R2}; fastqc -o ../QC/fastqc -f fastq $R1 $R2; done
# mv *fastqc* ../QC/fastqc
# for R1 in *_R1_001.fastq.gz; do R2=${R1/_R1/_R2}; fastqc -o ../QC/fastqc -f fastq $R1 $R2; done
# for R1 in *_R1.fastq.gz; do R2=${R1/_R1/_R2}; fastqc -o ../QC/fastqc -f fastq $R1 $R2; done

conda deactivate


#***********************************************************************#
echo "RNASEQC"
echo ""

conda activate rnaseq

# time parallel -j 12 "rnaseqc $gtf_gene {} {= s/Aligned.sortedByCoord.out.bam/_RNA-SeQC/; =} --sample={= s/Aligned.sortedByCoord.out.bam//; =} --stranded='rf' " ::: *Aligned.sortedByCoord.out.bam

for i in *Aligned.sortedByCoord.out.bam; 
   do sample=${i/Aligned.sortedByCoord.out.bam/}; 
   rnaseqc $gtf_gene $i ${sample}_RNA-SeQC --sample=${sample} --stranded='rf'; 
done


conda deactivate


#***********************************************************************#
# Couverture totale a chaque position du bed
echo "CollectHsMetrics"
echo ""

conda activate gatk4


#***********************************************************************#
echo "CollectRnaSeqMetrics"
echo ""

# time parallel -j 12 "gatk CollectRnaSeqMetrics -I {} -O {= s/Aligned.sortedByCoord.out.bam/.RNAseqMetrics.txt/; =} \
#  --REF_FLAT $refflat -R $ref -STRAND SECOND_READ_TRANSCRIPTION_STRAND \
#  --RIBOSOMAL_INTERVALS /media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.rRNA.transcripts.interval_list" ::: *Aligned.sortedByCoord.out.bam

for i in *Aligned.sortedByCoord.out.bam; 
   do sample=${i/Aligned.sortedByCoord.out.bam/}; 
   gatk CollectRnaSeqMetrics -I $i -O ${sample}.RNAseqMetrics.txt \
   --REF_FLAT $refflat \
   -R $ref -STRAND SECOND_READ_TRANSCRIPTION_STRAND \
   --RIBOSOMAL_INTERVALS /media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/gencode.v43.rRNA.transcripts.interval_list; 
done

conda deactivate


# #***********************************************************************#
echo "salmon"
echo ""

conda activate salmon


###### INDEX ######
# salmon index -t '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts.fa' \
# -i '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts-salmon.idx'

# salmon index -t '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts.fa' \
#    -i '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts-salmon-format.idx' --gencode

### Generating a decoy-aware transcriptome ###
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
# grep "^>" <(gunzip -c /media/jbogoin/Data1/References/RNA-seq/hg38/salmon/GRCh38.primary_assembly.genome.fa.gz) \
#    | cut -d " " -f 1 > /media/jbogoin/Data1/References/RNA-seq/hg38/salmon/decoys.txt
# sed -i.bak -e 's/>//g' /media/jbogoin/Data1/References/RNA-seq/hg38/salmon/decoys.txt
# cat '/media/jbogoin/Data1/References/RNA-seq/hg38/salmon/gencode.v47.transcripts.fa.gz' \
#    '/media/jbogoin/Data1/References/RNA-seq/hg38/salmon/GRCh38.primary_assembly.genome.fa.gz' \
#    > '/media/jbogoin/Data1/References/RNA-seq/hg38/salmon/gentrome.fa.gz'
# salmon index -t '/media/jbogoin/Data1/References/RNA-seq/hg38/salmon/gentrome.fa.gz' \
#    -d '/media/jbogoin/Data1/References/RNA-seq/hg38/salmon/decoys.txt' \
#    -p 12 i '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts-salmon-format-decoys.idx' --gencode


COUNT
for R1 in *_R1_001.fastq.gz; 
#for R1 in *_R1.fastq.gz; 
   do R2=${R1/_R1/_R2};
   sample=${R1%%_*};
   salmon quant -i '/media/jbogoin/Data1/References/RNA-seq/hg38/gencode.v38.transcripts-salmon-format.idx' \
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
mv Stats ../QC
mv *_RNA-SeQC ../QC
mv *.RNAseqMetrics.txt ../QC
mv *.hsMetrics.txt ../QC
mv *pertargetcoverage* ../QC
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
# mv *pertargetcoverage* pertargetcoverage
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
echo "DROP"
echo ""

bash ~/SCRIPTS/RNA-Seq/DROP/drop.sh


echo ""
echo "rnaseq_total.sh job done!"
echo ""
