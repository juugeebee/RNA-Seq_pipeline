#!/usr/bin/sudo bash

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
echo "DROP"
echo ""

bash ~/SCRIPTS/RNA-Seq/DROP/drop_96.sh


echo ""
echo "rnaseq_total_96.sh job done!"
echo ""
