#!/bin/bash

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh

echo ""
echo "rnaseq_tota_96.sh start"
echo ""


#########################################
## A LANCER DANS LE DOSSIER RACINE DU RUN
#########################################


#***********************************************************************#
echo "STAR"
echo ""

bash ~/SCRIPTS/RNA-Seq/STAR.sh


# ***********************************************************************#
echo "Tri des fichiers BAM par tissus"


BAM_DIR="./BAM"
SCRIPT_DIR="~/SCRIPTS/RNA-Seq/DROP"


cd "$BAM_DIR"

# Fichier de résumé demandé
SUMMARY="bam_par_tissus.txt"

declare -A patients_tissue

echo "Création des dossiers et liens symboliques..."

for bam in *Aligned.sortedByCoord.out.bam; do
    [[ -e "$bam" ]] || continue

    filename=$(basename "$bam")

    # Extraction du tissu
    tissue=$(echo "$filename" | sed -E 's/.*-([^_]+)Aligned.sortedByCoord.out.bam/\1/')

    # Extraction du patient
    patient=$(echo "$filename" | sed -E "s/-${tissue}_Aligned.*//")

    # Si le tissu n'est pas renseigné
    if [[ "$tissue" == "$filename" || -z "$tissue" ]]; then
        echo "ATTENTION : tissu absent pour $filename --> classement dans BAM-PAX"
        tissue="PAX"

        patient=$(echo "$filename" | sed -E 's/_Aligned.*//')
    fi

    outdir="BAM-${tissue}"

    mkdir -p "$outdir"

    # Création des liens symboliques
    ln -sf "../$bam" "$outdir/$bam"

    if [[ -f "${bam}.bai" ]]; then
        ln -sf "../${bam}.bai" "$outdir/${bam}.bai"
    fi

    # Mémorisation patient unique par tissu
    patients_tissue["${tissue},${patient}"]=1

done


# Comptage des patients uniques
declare -A count

for key in "${!patients_tissue[@]}"; do
    tissue=${key%%,*}

    if [[ -z "${count[$tissue]+x}" ]]; then
        count[$tissue]=1
    else
        ((count[$tissue]++))
    fi
done


# Création du fichier de synthèse
echo -e "Tissu\tNombre_patients" > "$SUMMARY"

for tissue in $(printf "%s\n" "${!count[@]}" | sort); do
    echo -e "${tissue}\t${count[$tissue]}" >> "$SUMMARY"
done


echo
echo "======================================"
echo "Résumé créé : $SUMMARY"
echo "======================================"

cat "$SUMMARY"


echo ""
echo "Lancement des scripts DROP..."


for dir in BAM-*; do

    tissue=${dir#BAM-}

    case "$tissue" in
        PAX)
            script="$SCRIPT_DIR/drop_96.sh"
            ;;
        EDTA)
            script="$SCRIPT_DIR/drop_96.sh"
            ;;
        FIBRO)
            script="$SCRIPT_DIR/drop_96.sh"
            ;;
        Emetine-FIBRO)
            script="$SCRIPT_DIR/drop_96.sh"
            ;;
        noEmetine-FIBRO)
            script="$SCRIPT_DIR/drop_96.sh"
            ;;
        Cyclo-FIBRO)
            script="$SCRIPT_DIR/drop_96.sh"
            ;;
        noCyclo-FIBRO)
            script="$SCRIPT_DIR/drop_96.sh"
            ;;
        ACD)
            script="$SCRIPT_DIR/drop_96.sh"
            ;;
        PEAU)
            script="$SCRIPT_DIR/drop_96.sh"
            ;;
        MUSCLE)
            script="$SCRIPT_DIR/drop_96.sh"
            ;; 
        TISSU)
            script="$SCRIPT_DIR/drop_96.sh"
            ;;                     
        *)
            echo "Aucun script défini pour le tissu '$tissue'."
            echo "Lancement du script PAX par défaut."
            script="$SCRIPT_DIR/drop_96.sh"
            ;;
    esac


    echo "=== $tissue ==="

    (
        cd "$dir"
        bash "$script"
    )

done


echo ""
echo "rnaseq_total_96.sh job done!"
echo ""
