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
head ./GRCh37_latest_genomic.gff
grep ">" GRCh37_latest_genomic.fna | grep -v "unplaced genomic contig"| grep -v "unlocalized genomic contig" | grep -v "genomic patch"| grep -v "alternate locus" | sed 's/^>//' > chr_names
seqtk subseq GRCh37_latest_genomic.fna chr_names > GRCh37_latest_genomic.sub.fna
cut -f 1 -d ' ' chr_names > chr_names_id
head -n 9 GRCh37_latest_genomic.gff > GRCh37_latest_genomic.sub.gff
grep -f chr_names_id GRCh37_latest_genomic.gff >> GRCh37_latest_genomic.sub.gff


#Configure SpliceLauncher with INSTALL mode

mkdir /media/jbogoin/Data1/References/RNA-seq/hg19/SpliceLauncher/refSpliceLauncher


bash ./SpliceLauncher.sh --runMode INSTALL \
-O '/media/jbogoin/Data1/References/RNA-seq/hg19/SpliceLauncher/refSpliceLauncher' \
--fasta '/media/jbogoin/Data1/References/RNA-seq/hg19/SpliceLauncher/GRCh37_latest_genomic.fna' \
--gff '/media/jbogoin/Data1/References/RNA-seq/hg19/SpliceLauncher/GRCh37_latest_genomic.gff' \
--STAR /home/jbogoin/miniconda3/envs/SpliceLauncher/bin/STAR \
--samtools /home/jbogoin/miniconda3/envs/SpliceLauncher/bin/samtools \
--bedtools /home/jbogoin/miniconda3/envs/SpliceLauncher/bin/bedtools


#Running the SpliceLauncher tests

bash '/home/jbogoin/SpliceLauncher/SpliceLauncher.sh' --runMode Align,Count,SpliceLauncher -F . -O ./testSpliceLauncher/