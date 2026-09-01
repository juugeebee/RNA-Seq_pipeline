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
# snakemake exportCounts --cores 4 --max-threads 24 --latency-wait 50 --resources mem_mb=100000 --rerun-incomplete > drop_export_counts.log


# echo ""
# echo "Lancement de FRASER2"
# echo ""
# snakemake aberrantSplicing --cores 1 --max-threads 24 --latency-wait 50 --resources mem_mb=100000 > drop_aberrantSplicing.log


echo ""
echo "Lancement d'OUTRIDER"
echo ""
# python ~/SCRIPTS/RNA-Seq/DROP/beforeOUTRIDER.py
# snakemake aberrantExpression --cores 4 --max-threads 24 --latency-wait 50 --resources mem_mb=100000 > drop_aberrantExpression.log


snakemake aberrantExpression --cores 4 --max-threads 24 \
  --latency-wait 50 --resources mem_mb=100 \
  2>&1 | tee drop_aberrantExpression.log


conda deactivate


echo ""
echo "Annotations des fichiers"


cd ..
python ~/SCRIPTS/RNA-Seq/DROP/prepare_annotation.py


### OURTIDER ###


OUTRIDER_DIR="./drop/output/processed_results/aberrant_expression/v48/outrider/outrider"


[[ -d "$OUTRIDER_DIR" ]] || {
	    echo "Erreur : le répertoire '$OUTRIDER_DIR' n'existe pas." >&2
    exit 1
}

cd "$OUTRIDER_DIR" || exit 1


header=$(head -n 1 OUTRIDER_results.tsv)
total=$(($(wc -l < OUTRIDER_results.tsv) - 1))
half=$((total / 2))

# première moitié
{ echo "$header"; tail -n +2 OUTRIDER_results.tsv | head -n "$half"; } > OUTRIDER_results_partie_1.tsv

# deuxième moitié
{ echo "$header"; tail -n +$((half + 2)) OUTRIDER_results.tsv; } > OUTRIDER_results_partie_2.tsv


cd ../../../../../../..


python ~/SCRIPTS/RNA-Seq/DROP/gene_annotation_96_outrider.py


python ~/SCRIPTS/RNA-Seq/DROP/outrider_volcano_plots.py


### FRASER2 ###


FRASER2_DIR="./drop/output/processed_results/aberrant_splicing/results/v48/fraser/fraser2/"


[[ -d "$FRASER2_DIR" ]] || {
	    echo "Erreur : le répertoire '$FRASER2_DIR' n'existe pas." >&2
    exit 1
}

cd "$FRASER2_DIR" || exit 1


header=$(head -n 1 results_gene_all.tsv)
total=$(($(wc -l < results_gene_all.tsv) - 1))
half=$((total / 2))

# première moitié
{ echo "$header"; tail -n +2 results_gene_all.tsv | head -n "$half"; } > FRASER2_results_partie_1.tsv

# deuxième moitié
{ echo "$header"; tail -n +$((half + 2)) results_gene_all.tsv; } > FRASER2_results_partie_2.tsv


cd ../../../../../../../..


python ~/SCRIPTS/RNA-Seq/DROP/gene_annotation_96_fraser2.py



echo ""
echo "drop.sh job done!"
echo ""
