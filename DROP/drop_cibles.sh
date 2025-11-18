#!/usr/bin/env bash


source ~/miniconda3/etc/profile.d/conda.sh


echo ""
echo "drop_cibles.sh start"
echo ""


mkdir -p drop
cd drop


python ~/SCRIPTS/RNA-Seq/DROP/config_file.py
python ~/SCRIPTS/RNA-Seq/DROP/sample_annotation_OA.py
# python ~/SCRIPTS/RNA-Seq/DROP/sample_annotation_NG.py


conda activate drop_env


if [ -d "output" ];then
    echo ''
    echo "Le dossier output existe !";
    echo ''
    rm -Rf output;
    snakemake aberrantSplicing --unlock --default-resources
    snakemake aberrantExpression --unlock --default-resources
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
snakemake aberrantSplicing --cores 4 --max-threads 24 --latency-wait 50 --resources mem_mb=100 --rerun-triggers mtime -k > drop_aberrantSplicing.log
#snakemake aberrantSplicing --cores 6 --max-threads 24 --default-resources "tmpdir='/media/jbogoin/Data4/tmp'"

echo ""
echo "Lancement d'OUTRIDER"
echo ""
snakemake aberrantExpression --cores 4 --max-threads 24 --latency-wait 50 --resources mem_mb=100 --rerun-triggers mtime -k  > drop_aberrantExpression.log
#snakemake aberrantExpression --cores 6 --max-threads 24 --default-resources "tmpdir='/media/jbogoin/Data4/tmp'"

conda deactivate


echo ""
echo "Annotations des fichiers"
cd ..
python ~/SCRIPTS/RNA-Seq/DROP/prepare_annotation.py
python ~/SCRIPTS/RNA-Seq/DROP/gene_annotation.py


conda deactivate

echo ""
echo "drop_cibles.sh job done!"
echo ""
