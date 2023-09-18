#!/usr/bin/env bash


source ~/miniconda3/etc/profile.d/conda.sh


echo ""
echo "drop.sh start"
echo ""


mkdir -p drop
cd drop


python ~/SCRIPTS/RNA-Seq/DROP/config_file.py
python ~/SCRIPTS/RNA-Seq/DROP/sample_annotation.py


conda activate drop_env


if [ -d "output" ];then
    echo "Le dossier output existe !";
    rm -Rf output;
    snakemake aberrantSplicing --unlock
else :
    echo "Le dossier output n'existe pas !";
    drop init;
    drop update
fi


snakemake aberrantSplicing --cores 16 --latency-wait 10


conda deactivate



echo ""
echo "drop.sh job done!"
echo ""