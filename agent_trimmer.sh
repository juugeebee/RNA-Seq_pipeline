#!/usr/bin/env bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate agent_env


echo ""
echo "agent_trimmer.sh start"
echo ""


mkdir -p Trimmer

for R1 in *_R1_001.fastq.gz; 
    do R2=${R1/_R1/_R2}; 
    SAMPLE=${R1%%_*}; 
    java -jar '/home/jbogoin/AGeNT_3.0.4/agent3.0/lib/trimmer-3.0.3.jar' \
        -fq1 $R1 -fq2 $R2 -v2 -out ./Fastq_trimmed/$SAMPLE;

done


echo ""
echo "agent_trimmer.sh job done!"
echo ""