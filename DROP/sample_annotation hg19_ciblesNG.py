import pandas
from os import listdir
import os

chemin = os.getcwd()
chemin_minus_drop = chemin.strip("/drop")
chemin_l = chemin.split('/')
runid = chemin_l[-1]


rnaid = []
rna_bam_file = []


for f in listdir("../BAM") :
    if ('Aligned.sortedByCoord.out.bam' in f) and ('bai' not in f) :
        rna_bam_file.append('/' + chemin_minus_drop + '/BAM/' + f)
        id_l = f.split('Aligned.sortedByCoord')
        sample = id_l[0]
        rnaid.append(sample)


df = pandas.DataFrame(columns = ['RNA_ID', 'RNA_BAM_FILE', 'DROP_GROUP', 'PAIRED_END', \
    'STRAND', 'COUNT_MODE', 'COUNT_OVERLAPS', 'INDIVIDAL_ID', 'DNA_VCF_FILE', 'DNA_ID'])


rnaid.append('6622NG002672_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6622NG002672Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG001598_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623NG001598Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000227_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000227Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000235_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000235Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000236_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000236Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000237_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000237Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000244_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000244Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000245_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000245Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000253_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000253Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000254_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000254Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000259_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000259Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000324_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000324Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000325_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000325Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000326_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000326Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000328_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000328Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000329_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000329Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000352_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000352Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000353_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000353Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000370_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000370Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000393_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000393Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000394_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000394Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000395_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000395Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000396_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000396Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000397_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000397Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000426_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000426Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000445_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000445Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000446_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000446Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000490_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000490Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000491_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000491Aligned.sortedByCoord.out.bam')
rnaid.append('6623CY000537_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg38/PAXgene/6623CY000537Aligned.sortedByCoord.out.bam')


df['RNA_ID'] = rnaid
df['RNA_BAM_FILE'] = rna_bam_file
#df['DROP_GROUP'] = 'fraser,fraser_external,batch_0,group1,PAXgene'
df['DROP_GROUP'] = 'fraser,batch_0,group1,PAXgene'
df['PAIRED_END'] = 'TRUE'
df['STRAND'] = 'REVERSE'
df['COUNT_MODE'] = 'Union'
df['COUNT_OVERLAPS'] = 'FALSE'
df['INDIVIDAL_ID'] = rnaid
df['DNA_VCF_FILE'] = '/media/jbogoin/Data1/References/RNA-seq/hg38/DROP/qc_vcf_1000G_hg38.vcf.gz'
df['DNA_ID'] = rnaid


df.to_csv('sample_annotation.tsv', sep='\t', index=False)


print('\n#####\nFichier sample_annotation.tsv genere !\n')
