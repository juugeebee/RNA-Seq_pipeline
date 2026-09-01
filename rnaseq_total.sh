#!/usr/bin/sudo bash

source ~/miniconda3/etc/profile.d/conda.sh

echo ""
echo "rnaseq_total.sh start"
echo ""


#########################################
## A LANCER DANS LE DOSSIER RACINE DU RUN
#########################################


#***********************************************************************#
echo "STAR"
echo ""

bash ~/SCRIPTS/RNA-Seq/STAR.sh


#***********************************************************************#
echo "DROP"
echo ""

bash ~/SCRIPTS/RNA-Seq/DROP/drop.sh


echo ""
echo "rnaseq_total.sh job done!"
echo ""


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