#!/usr/bin/env bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate fastq_bam_env

###############################################################################
##### BCL to FASTQ

# bcl2fastq --barcode-mismatches 0 --minimum-trimmed-read-length 35 --no-lane-splitting -R "." \
# --sample-sheet "SampleSheet.csv" -o ./Fastq_raw -r 4 -p 12 -w 4

conda deactivate

# cd Fastq_raw
# bash ~/SCRIPTS/RNA-seq/agent_trimmer.sh

###############################################################################
#### ALIGNEMENT
### STAR (Spliced Transcripts Alignment to a Reference)

## GENERATING GENOME INDEXES

conda deactivate
conda activate rnaseq

genome_dir='/media/jbogoin/Data1/References/RNA-seq/STAR'
ref='/media/jbogoin/Data1/References/fa_hg19/rna-seq/GRCh37.primary_assembly.genome.fa'
gtf_file='/media/jbogoin/Data1/References/fa_hg19/rna-seq/gencode.v41lift37.annotation.gtf'
refflat='/media/jbogoin/Data1/References/RNA-seq/refFlat_hg19.txt'


# STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $genome_dir\
#  --genomeFastaFiles $ref --sjdbGTFfile $gtf_file --sjdbOverhang 71

   

## RUNNING MAPPING JOB
# cd Fastq_trimmed

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

# for i in *Aligned.sortedByCoord.out.bam; do samtools index -@ 16 $i; done


# ## QC


# for i in *Aligned.sortedByCoord.out.bam; 
#     do sample=${i/Aligned.sortedByCoord.out.bam/}; 
#     rnaseqc $gtf_file $i ${sample}_RNA-SeQC --sample=${sample} --stranded=rf; 
# done


# conda deactivate
# conda activate gatk4

# for i in *Aligned.sortedByCoord.out.bam; 
#     do sample=${i%Aligned.sortedByCoord.out.bam}; 
#     gatk CollectRnaSeqMetrics -I $i -O ${sample}.RNAseqMetrics.txt \
#     --REF_FLAT $refflat \
#     -R /$REF \
#     -STRAND SECOND_READ_TRANSCRIPTION_STRAND; 
# done



# mkdir ../QC
# mv *.gct ../QC
# mv *_Metrics ../QC
# mv *.tsv ../QC

# cd ../QC
# multiqc .