#!/usr/bin/env bash

source ~/miniconda3/etc/profile.d/conda.sh

echo ""
echo "rnaseq_cibles.sh start"
echo ""


echo "Choisissez l'UF :"

select ENV in NG OA; do
    case $ENV in
        NG)
            export ENV
            export target='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions_liftover_hg38_ucsc.bed'
            export target_il='/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions_liftover_hg38_ucsc.interval_list'
            break
            ;;
        OA)
            export ENV
            export target='/media/jbogoin/Data1/References/cibles_panel_OA/ONCO_BED_RNASEQ_GENE_DIAG_CODING_EXON_hg38.bed'
            export target_il='/media/jbogoin/Data1/References/cibles_panel_OA/ONCO_BED_RNASEQ_GENE_DIAG_CODING_EXON_hg38.interval_list'
            break
            ;;
        *)
            echo "Choix invalide."
            ;;
    esac
done


echo ""
echo "UF sélectionnée : $ENV"
echo "Fichier cible : $target"
echo ""



#***********************************************************************#
# #TRIMMING DES ADATATEURS
# echo "TRIMMER"
# echo ""

# conda activate rnaseq
# bash ~/SCRIPTS/RNA-Seq/agent_trimmer.sh
# cd ..


#***********************************************************************#
echo "STAR"
echo ""

bash ~/SCRIPTS/RNA-Seq/STAR.sh


#***********************************************************************#
echo "pertargetcoverage"
echo ""

conda activate gatk4

cd BAM


#CIBLES
for i in *Aligned.sortedByCoord.out.bam; 
    do sample=${i%Aligned.sortedByCoord.out.bam}; 
    gatk CollectHsMetrics \
    -I $i \
    -O ${sample}.hsMetrics.txt \
    -R /media/jbogoin/Data1/References/fa_hg38/hg38_rnaseq/GRCh38.v48.primary_assembly.chr25.fa \
    --BAIT_INTERVALS $target_il \
    --TARGET_INTERVALS $target_il \
    --PER_TARGET_COVERAGE ${sample}.pertargetcoverage_cibles.txt;
done


cd ..


# ***********************************************************************#
echo "DROP"
echo ""


bash ~/SCRIPTS/RNA-Seq/DROP/drop_cibles.sh


echo ""
echo "rnaseq_cibles.sh job done!"
echo ""