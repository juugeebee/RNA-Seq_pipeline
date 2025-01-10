rmats.py --b2 '/media/jbogoin/Data2/rMATS/Il-str-mRNA-D_21112024/6624CY001245/rmats_control_bams.txt' \
    --b1 '/media/jbogoin/Data2/rMATS/Il-str-mRNA-D_21112024/6624CY001245/RNU5_bams.txt' \
    --gtf '/media/jbogoin/Data1/References/RNA-seq/hg38/rmats/gencode.v47.annotation.gtf' \
    -t paired --readLength 101 --od '/media/jbogoin/Data2/rMATS/Il-str-mRNA-D_21112024/6624CY001245/rmats_out' \
    --tmp '/media/jbogoin/Data2/rMATS/Il-str-mRNA-D_21112024/6624CY001245/tmp' \
    --anchorLength 1 --nthread 12 --libType fr-firststrand --task both --novelSS --variable-read-length --allow-clipping


rmats.py --b2 '/media/jbogoin/Data2/rMATS/Il-str-mRNA-D_21112024/6624CY001359/rmats_control_bams.txt' \
    --b1 '/media/jbogoin/Data2/rMATS/Il-str-mRNA-D_21112024/6624CY001359/RNU5_bams.txt' \
    --gtf '/media/jbogoin/Data1/References/RNA-seq/hg38/rmats/gencode.v47.annotation.gtf' \
    -t paired --readLength 101 --od '/media/jbogoin/Data2/rMATS/Il-str-mRNA-D_21112024/6624CY001359/rmats_out' \
    --tmp '/media/jbogoin/Data2/rMATS/Il-str-mRNA-D_21112024/6624CY001359/tmp' \
    --anchorLength 1 --nthread 12 --libType fr-firststrand --task both --novelSS --variable-read-length --allow-clipping


# Benjamin
python /data1/benjamin/miniconda3/envs/rmats/bin/rmats.py \
    --gtf /data2/Exome_analysis_BC/Server_S_exome/bank/annotation/ensembl/currentHomo_sapiens.GRCh38.106.gtf \
    --tmp /data1/tmp --readLength 76 --b1 b1.txt --b2 b2.txt --od results_21062024_13RNU -t paired \
    --anchorLength 1 --nthread 5 --libType fr-firststrand --task both --novelSS --variable-read-length --allow-clipping