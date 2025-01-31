#!/usr/bin/env bash


source ~/miniconda3/etc/profile.d/conda.sh


echo ""
echo "drop..sh start"
echo ""


mkdir -p drop
cd drop


python ~/SCRIPTS/RNA-Seq/DROP/config_file.py
python ~/SCRIPTS/RNA-Seq/DROP/sample_annotation.py


conda activate drop_env


if [ -d "output" ];then
    echo ''
    echo "Le dossier output existe !";
    echo ''
    rm -Rf output;
    snakemake aberrantSplicing --unlock
    snakemake aberrantExpression --unlock
else :
    echo ''
    echo "Le dossier output n'existe pas !";
    echo ''
    drop init;
    drop update
fi

echo ""
echo "Lancement de FRASER2"
echo ""
snakemake aberrantSplicing --cores 4 --max-threads 24 --latency-wait 50 --resources mem_mb=100 > drop_aberrantSplicing.log


echo ""
echo "Lancement d'OUTRIDER"
echo ""
snakemake aberrantExpression --cores 4 --max-threads 24 --latency-wait 50 --resources mem_mb=100 > drop_aberrantExpression.log


conda deactivate


echo ""
echo "Annotations des fichiers"
cd ..
python ~/SCRIPTS/RNA-Seq/DROP/prepare_annotation.py
python ~/SCRIPTS/RNA-Seq/DROP/gene_annotation.py


echo ""
echo "drop.sh job done!"
echo ""
