#### SpliceLauncher ####
########################


conda activate SpliceLauncher

cd ~/SpliceLauncher
git pull


#Configure SpliceLauncher with INSTALL mode 
bash ./SpliceLauncher.sh --runMode INSTALL -O ./refSpliceLauncher \
    --STAR /home/jbogoin/miniconda3/envs/SpliceLauncher/bin/STAR \
    --samtools /home/jbogoin/miniconda3/envs/SpliceLauncher/bin/samtools \
    --bedtools /home/jbogoin/miniconda3/envs/SpliceLauncher/bin/bedtools \
    --gff /media/jbogoin/Data1/References/RNA-seq/hg38/SpliceLauncher/GRCh38_latest_genomic.gff \
    --threads 4 \
    --fasta /media/jbogoin/Data1/References/RNA-seq/hg38/SpliceLauncher/GRCh38_latest_genomic.fna


#Running the SpliceLauncher tests
bash ./SpliceLauncher.sh --runMode Align,Count,SpliceLauncher -F ./dataTest/fastq/ -O ./testSpliceLauncher/

