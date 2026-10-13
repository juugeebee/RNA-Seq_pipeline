#### SpliceLauncher ####
########################
source ~/miniconda3/etc/profile.d/conda.sh

echo ""
echo "SpliceLauncher start"
echo ""

conda activate SpliceLauncher


#################################################################################################
# #INSTALL
# cd ~/SpliceLauncher
# git pull


# bash ./SpliceLauncher.sh --runMode INSTALL -O ./refSpliceLauncher \
#     --STAR /home/jbogoin/miniconda3/envs/SpliceLauncher/bin/STAR \
#     --samtools /home/jbogoin/miniconda3/envs/SpliceLauncher/bin/samtools \
#     --bedtools /home/jbogoin/miniconda3/envs/SpliceLauncher/bin/bedtools \
#     --gff /media/jbogoin/Data1/References/RNA-seq/hg38/SpliceLauncher/GRCh38_latest_genomic.gff \
#     --fasta /media/jbogoin/Data1/References/RNA-seq/hg38/SpliceLauncher/GRCh38_latest_genomic.fna \
#     --assembly_report /media/jbogoin/Data1/References/RNA-seq/hg38/SpliceLauncher/GRCh38_latest_assembly_report.txt \
#     --threads 4 2>&1 | tee ~/install_log.txt
#################################################################################################


#Renommage des bam et bai
cd /media/jbogoin/Data4/Il-str-mRNA-D_17042024-TBR1/BAM

for f in *Aligned.sortedByCoord.out.bam*; do
    mv "$f" "${f/Aligned/.Aligned}"
done

ls *.bam *.bai 2>/dev/null


#Running SpliceLauncher
cd ~/SpliceLauncher/

bash ./SpliceLauncher.sh --runMode Count,SpliceLauncher \
    -B /media/jbogoin/Data4/Il-str-mRNA-D_17042024-TBR1/BAM \
    -O /media/jbogoin/Data4/Il-str-mRNA-D_17042024-TBR1/SpliceLauncher \
    -t 4 \
    --Graphics


conda deactivate SpliceLauncher


echo ""
echo "SpliceLauncher job done!"
echo ""