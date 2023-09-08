#!/usr/bin/env bash

source ~/miniconda3/etc/profile.d/conda.sh


echo ""
echo "drop.sh start"
echo ""


mkdir -p DROP
cd DROP


python ~/SCRIPTS/RNA-Seq/DROP/config_file.py
python ~/SCRIPTS/RNA-Seq/DROP/sample_annotation.py


conda activate drop_env


drop init
drop update
snakemake aberrantSplicing --cores 10


conda deactivate



echo ""
echo "drop.sh job done!"
echo ""