conda activate rnaseq

###############################################################################
##### BCL to FASTQ

bcl2fastq --barcode-mismatches 0 --minimum-trimmed-read-length 35 --no-lane-splitting -R "." \
--sample-sheet "SampleSheet.csv" -o ./Fastq -r 4 -p 12 -w 4


###############################################################################
#### ALIGNEMENT
### STAR (Spliced Transcripts Alignment to a Reference)

## GENERATING GENOME INDEXES

genome_dir='/media/jbogoin/Data1/References/RNA-seq/STAR'
ref='/media/jbogoin/Data1/References/fa_hg38/hg38-rnaseq/GRCh38.primary_assembly.genome.fa'
gtf_file='/media/jbogoin/Data1/References/RNA-seq/STAR/UCSC_mRNA_EST_RefSeq.gtf'

STAR --runThreadN 24 --runMode genomeGenerate --genomeDir $genome_dir\
 --genomeFastaFiles $ref --sjdbGTFfile $gtf_file --sjdbOverhang 99
    

## RUNNING MAPPING JOB

start_end='/media/jbogoin/Data1/References/RNA-seq/STAR/sjdbList.fromGTF.out.tab'

cd Fastq

rm -Rf ../BAM


## BAM
for R1 in *_R1_001.fastq.gz; 
    do R2=${R1/_R1/_R2}; 
    
    SAMPLE=${R1%%_*};
    
    STAR --runThreadN 12 \
        --genomeDir $genome_dir \
        --readFilesIn $R1 $R2 \
        --readFilesCommand zcat \
        --outFileNamePrefix ../BAM/$SAMPLE \
        --outSAMtype BAM SortedByCoordinate\
        --sjdbGTFfile $gtf_file \
        --sjdbFileChrStartEnd $start_end;

done

cd ../BAM


## DUPLICATES
for BAM in *.bam;
    do SAMPLE=${BAM%A*};

    STAR --runThreadN 12 \
        --limitBAMsortRAM 30000000000 \
        --runMode inputAlignmentsFromBAM \
        --bamRemoveDuplicatesType UniqueIdentical \
        --inputBAMfile $BAM \
        --outFileNamePrefix ${SAMPLE}.dedup.bam;

done


## INDEX
parallel samtools index ::: *.dedup.bamProcessed.out.bam

conda deactivate


###############################################################################
#### QC


## GATK 
conda activate gatk4

# echo "Download refFlat.txt annotation file"
# wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz
# gunzip refFlat.txt.gz

for BAM in *.dedup.bamProcessed.out.bam;
    do SAMPLE=${BAM%.*};

gatk CollectRnaSeqMetrics \
      -I $BAM \
      -O $SAMPLE.RNA_Metrics \
      --REF_FLAT '/media/jbogoin/Data1/References/RNA-seq/refFlat.txt' \
      -STRAND FIRST_READ_TRANSCRIPTION_STRAND;

done

conda deactivate


## RNASeQC
conda activate rnaseq

# https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model
# cd ~/gtex-pipeline/gene_model
# python3 collapse_annotation.py \
#     '/media/jbogoin/Data1/References/RNA-seq/gencode.v41.basic.annotation.gtf' 
#     '/media/jbogoin/Data1/References/RNA-seq/gencode.v41.genes.gtf'


for BAM in *.dedup.bamProcessed.out.bam;
    do SAMPLE=${BAM%.*};

    rnaseqc '/media/jbogoin/Data1/References/RNA-seq/gencode.v41.genes.gtf' \
        $BAM . ;

done

mkdir ../QC
mv *.gct ../QC
mv *_Metrics ../QC
mv *.tsv ../QC

cd ../QC
multiqc .