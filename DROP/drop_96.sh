#!/usr/bin/env bash


source ~/miniconda3/etc/profile.d/conda.sh


echo ""
echo "drop.sh start"
echo ""


mkdir -p drop
cd drop


python ~/SCRIPTS/RNA-Seq/DROP/config_file_96.py
python ~/SCRIPTS/RNA-Seq/DROP/sample_annotation_96.py


conda activate drop_env


if [ -d "output" ];then
    echo ''
    echo "Le dossier output existe !";
    echo ''
    rm -Rf output;
    snakemake aberrantSplicing --unlock
    snakemake aberrantExpression --unlock
    # snakemake exportCounts --unlock
else :
    echo ''
    echo "Le dossier output n'existe pas !";
    echo ''
    drop init;
    drop update
fi


# echo ""
# echo "Lancement de l'external counts"
# echo ""
# snakemake exportCounts --cores 4 --max-threads 24 --latency-wait 50 --resources mem_mb=100 --rerun-incomplete > drop_export_counts.log


# echo ""
# echo "Lancement de FRASER2"
# echo ""
# snakemake aberrantSplicing --cores 4 --max-threads 24 --latency-wait 50 --resources mem_mb=100 > drop_aberrantExpression.log



echo ""
echo "Lancement d'OUTRIDER"
echo ""
snakemake aberrantExpression --cores 4 --max-threads 24 --latency-wait 50 --resources mem_mb=100 > drop_aberrantExpression.log


conda deactivate


echo ""
echo "Annotations des fichiers"
cd ..
python ~/SCRIPTS/RNA-Seq/DROP/prepare_annotation.py


cd ./drop/output/processed_results/aberrant_expression/v48/outrider/outrider

header=$(head -n 1 OUTRIDER_results.tsv)
total=$(($(wc -l < OUTRIDER_results.tsv) - 1))
half=$((total / 2))

# première moitié
{ echo "$header"; tail -n +2 OUTRIDER_results.tsv | head -n "$half"; } > OUTRIDER_results_partie_1.tsv

# deuxième moitié
{ echo "$header"; tail -n +$((half + 2)) OUTRIDER_results.tsv; } > OUTRIDER_results_partie_2.tsv

cd ../../../../../../..


python ~/SCRIPTS/RNA-Seq/DROP/gene_annotation_96.py

echo ""
echo "drop.sh job done!"
echo ""
