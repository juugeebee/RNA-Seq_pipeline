#!/usr/bin/sudo bash

source ~/miniconda3/etc/profile.d/conda.sh

echo ""
echo "pertargetcoverage.sh start"
echo ""


ref='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/GRCh38.v43.primary_assembly.genome.fa'
targets_il='/media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/BED/genes_panel_MB_V11_hg38.interval_list'


#***********************************************************************#
echo "pertargetcoverage"
echo ""


conda activate gatk4


# Couverture totale a chaque position du bed
for i in *Aligned.sortedByCoord.out.bam; 
    do sample=${i/Aligned.sortedByCoord.out.bam/}; 
    gatk CollectHsMetrics \
    -I $i \
    -O ${sample}.hsMetrics.txt \
    -R $ref \
    --BAIT_INTERVALS $mb_il \
    --TARGET_INTERVALS $mb_il \
    --PER_TARGET_COVERAGE ${sample}.hsMetrics_pertargetcoverage.txt;
done


conda deactivate


echo ""
echo "pertargetcoverage.sh job done!"
echo ""
