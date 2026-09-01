#!/usr/bin/env bash


source ~/miniconda3/etc/profile.d/conda.sh


echo ""
echo "drop_export_counts.sh start"
echo ""


# mkdir -p drop
# cd drop


# python ~/SCRIPTS/RNA-Seq/DROP/config_file.py
# python ~/SCRIPTS/RNA-Seq/DROP/sample_annotation.py


conda activate drop_env


if [ -d "output" ];then
    echo ''
    echo "Le dossier output existe !";
    echo ''
    rm -Rf output;
    snakemake aberrantSplicing --unlock
    snakemake aberrantExpression --unlock
    snakemake exportCounts --unlock
else :
    echo ''
    echo "Le dossier output n'existe pas !";
    echo ''
    drop init;
    drop update
fi


echo ""
echo "Lancement d'OUTRIDER"
echo ""
snakemake aberrantExpression --cores 1 --max-threads 24 --latency-wait 50 --resources mem_mb=100000 > drop_aberrantExpression.log


echo ""
echo "Lancement de FRASER2"
echo ""
snakemake aberrantSplicing --cores 1 --max-threads 24 --latency-wait 50 --resources mem_mb=100000 > drop_aberrantSplicing.log


echo ""
echo "Lancement de l'external counts"
echo ""
snakemake exportCounts --cores 1 --max-threads 24 --latency-wait 50 --resources mem_mb=100000 > drop_export_counts.log



conda deactivate


echo ""
echo "drop.sh job done!"
echo ""