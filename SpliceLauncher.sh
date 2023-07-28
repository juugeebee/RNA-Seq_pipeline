#### SpliceLauncher ####
########################


#Installing SpliceLauncher

git clone https://github.com/raphaelleman/SpliceLauncher


#Conda environment

cd ./SpliceLauncher
conda env create -f environment.yml
conda activate spliceLauncher


#Reference files

mkdir -p '/media/jbogoin/Data11/References/RNA-seq/SpliceLauncher'
cd /media/jbogoin/Data1/References/RNA-seq/SpliceLauncher

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.fna.gz
gunzip ./GRCh37_latest_genomic.fna.gz

wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz
gunzip ./GRCh37_latest_genomic.gff.gz


#Configure SpliceLauncher with INSTALL mode

mkdir /media/jbogoin/Data1/References/RNA-seq/hg19/SpliceLauncher/refSpliceLauncher


bash ./SpliceLauncher.sh --runMode INSTALL \
-O '/media/jbogoin/Data1/References/RNA-seq/hg19/SpliceLauncher/refSpliceLauncher' \
--fasta '/media/jbogoin/Data1/References/RNA-seq/hg19/SpliceLauncher/GRCh37_latest_genomic.fna' \
--gff '/media/jbogoin/Data1/References/RNA-seq/hg19/SpliceLauncher/GRCh37_latest_genomic.gff' \
--STAR '/home/jbogoin/miniconda3/envs/SpliceLauncher/bin/STAR' \
--samtools '/home/jbogoin/miniconda3/envs/SpliceLauncher/bin/samtools' \
--bedtools '/home/jbogoin/miniconda3/envs/SpliceLauncher/bin/bedtools'


#Running the SpliceLauncher tests

bash '/home/jbogoin/SpliceLauncher/SpliceLauncher.sh' --runMode Align,Count,SpliceLauncher -F . -O ./testSpliceLauncher/