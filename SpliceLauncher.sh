sudo apt-get update
sudo apt-get install g++
sudo apt-get install make


#STAR

mkdir -p STAR
wget https://github.com/alexdobin/STAR/archive/2.7.0c.tar.gz
tar -xzf 2.7.0c.tar.gz
cd STAR-2.7.0c
cd STAR/source
make STAR


#SAMTOOLS

mkir -p samtools
wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2
tar -xf samtools-1.16.1.tar.bz2
cd samtools-1.16.1
./configure --prefix=./configure --prefix='/home/jbogoin/samtools/samtools-1.16.1'
make
make install


#BEDTOOLS

mkdir -p bedtools
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
cd bedtools2
make


#SpliceLauncher

git clone https://github.com/raphaelleman/SpliceLauncher
cd ./SpliceLauncher
conda env create -f environment.yml
conda activate spliceLauncher


#Reference files

mkdir -p '/media/jbogoin/Data11/References/RNA-seq/SpliceLauncher'
cd /media/jbogoin/Data11/References/RNA-seq/SpliceLauncher
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

bash ./SpliceLauncher.sh --runMode INSTALL -O ./refSpliceLauncher \
    --fasta '/media/jbogoin/Data11/References/RNA-seq/SpliceLauncher/GRCh37_latest_genomic.fna' \
    --gff '/media/jbogoin/Data11/References/RNA-seq/SpliceLauncher/GRCh37_latest_genomic.gff' \
    --STAR '/home/jbogoin/STAR/STAR-2.7.0c/bin/Linux_x86_64_static/STAR' \
    --samtools '/home/jbogoin/samtools/samtools-1.16.1/samtools' \
    --bedtools '/home/jbogoin/bedtools/bedtools2/bin/bedtools'


#SpliceLauncher

bash '/home/jbogoin/SpliceLauncher/SpliceLauncher.sh' --runMode Align,Count,SpliceLauncher -F . -O ./testSpliceLauncher/