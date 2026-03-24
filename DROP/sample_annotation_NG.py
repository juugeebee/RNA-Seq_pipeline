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
    if ('_Aligned.sortedByCoord.out.bam' in f) and ('bai' not in f) :
        rna_bam_file.append('/' + chemin_minus_drop + '/BAM/' + f)
        id_l = f.split('Aligned.sortedByCoord')
        sample = id_l[0]
        rnaid.append(sample)


df = pandas.DataFrame(columns = ['RNA_ID', 'RNA_BAM_FILE', 'DROP_GROUP', 'PAIRED_END', \
    'STRAND', 'COUNT_MODE', 'COUNT_OVERLAPS', 'INDIVIDAL_ID', 'DNA_VCF_FILE', 'DNA_ID'])


# Get id and path for witness BAM file
rnaid.append('17NG002384-FSP_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/17NG002384-FSP_Aligned.sortedByCoord.out.bam')
rnaid.append('18NG001626_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/18NG001626_Aligned.sortedByCoord.out.bam')
rnaid.append('18NG001723-PAR_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/18NG001723-PAR_Aligned.sortedByCoord.out.bam')
rnaid.append('18NG002044-ALD_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/18NG002044-ALD_Aligned.sortedByCoord.out.bam')
rnaid.append('21NG002225-EEE_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/21NG002225-EEE_Aligned.sortedByCoord.out.bam')
rnaid.append('6619NG001842-EEE_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6619NG001842-EEE_Aligned.sortedByCoord.out.bam')
rnaid.append('6620NG001457-FSP_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6620NG001457-FSP_Aligned.sortedByCoord.out.bam')
rnaid.append('6622NG000396-ITD_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6622NG000396-ITD_Aligned.sortedByCoord.out.bam')
rnaid.append('6622NG001445-FSP_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6622NG001445-FSP_Aligned.sortedByCoord.out.bam')
rnaid.append('6622NG001886-FHM_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6622NG001886-FHM_Aligned.sortedByCoord.out.bam')
rnaid.append('6622NG003246_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6622NG003246_Aligned.sortedByCoord.out.bam')
rnaid.append('6622NG003716_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6622NG003716_Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000155_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG000155_Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000221_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG000221_Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000604_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG000604_Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000823_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG000823_Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000990_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG000990_Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG001056_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG001056_Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG001599_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG001599_Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG002822-F_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG002822-FAligned.sortedByCoord.out.bam')
rnaid.append('6623NG002917-F_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG002917-FAligned.sortedByCoord.out.bam')
rnaid.append('6623NG003378-F_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG003378-FAligned.sortedByCoord.out.bam')
rnaid.append('6623NG003412-F_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6623NG003412-FAligned.sortedByCoord.out.bam')
rnaid.append('6624NG000450-F_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6624NG000450-FAligned.sortedByCoord.out.bam')
rnaid.append('6624NG000698-M_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6624NG000698-MAligned.sortedByCoord.out.bam')
rnaid.append('6624NG000891-F_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6624NG000891-FAligned.sortedByCoord.out.bam')
rnaid.append('6624NG001409-F_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6624NG001409-FAligned.sortedByCoord.out.bam')
rnaid.append('6624NG001620-F_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6624NG001620-FAligned.sortedByCoord.out.bam')
rnaid.append('6624NG002428-M_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_NG_hg38/PAXgene/6624NG002428-MAligned.sortedByCoord.out.bam')


df['RNA_ID'] = rnaid
df['RNA_BAM_FILE'] = rna_bam_file
df['DROP_GROUP'] = 'patient,fraser2,outrider,mae,variant_calling,PAXgene'
df['PAIRED_END'] = 'TRUE'
df['STRAND'] = 'REVERSE'
df['COUNT_MODE'] = 'Union'
df['COUNT_OVERLAPS'] = 'FALSE'
df['INDIVIDAL_ID'] = rnaid
df['DNA_VCF_FILE'] = '/media/jbogoin/Data1/References/RNA-seq/hg38/DROP/qc_vcf_1000G_hg38.vcf.gz'
df['DNA_ID'] = rnaid


df.to_csv('sample_annotation.tsv', sep='\t', index=False)


print('\n#####\nFichier sample_annotation_NG.tsv genere !\n')
