import pandas
from os import listdir
import os


run_l = os.getcwd()
run = run_l.split('/')
runid = run[-1]


rnaid = []
rna_bam_file = []


for f in listdir("../BAM") :
    if ('Aligned.sortedByCoord.out.bam' in f) and ('bai' not in f) :
        rna_bam_file.append("../BAM/" + f)
        id_l = f.split('A')
        sample = id_l[0]
        rnaid.append(sample)


df = pandas.DataFrame(columns = ['RNA_ID', 'RNA_BAM_FILE', 'DROP_GROUP', 'PAIRED_END', \
    'STRAND', 'COUNT_MODE', 'COUNT_OVERLAPS', 'INDIVIDAL_ID', 'DNA_VCF_FILE', 'DNA_ID'])


df['RNA_ID'] = rnaid
df['RNA_BAM_FILE'] = rna_bam_file
df['DROP_GROUP'] = 'outrider,outrider_external,fraser,fraser_external,batch_0,group1,blood,WES'
df['PAIRED_END'] = 'FALSE'
df['STRAND'] = 'REVERSE'
df['COUNT_MODE'] = 'Union'
df['COUNT_OVERLAPS'] = 'FALSE'
df['INDIVIDAL_ID'] = rnaid
df['DNA_VCF_FILE'] = '/media/jbogoin/Data1/References/RNA-seq/hg38/DROP/qc_vcf_1000G_hg38.vcf.gz'
df['DNA_ID'] = rnaid


df.to_csv('sample_annotation.tsv', sep='\t', index=False)