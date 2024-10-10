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


rnaid.append('617NG002384-FSP_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/17NG002384-FSPAligned.sortedByCoord.out.bam')
rnaid.append('18NG001723-PAR_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/18NG001723-PARAligned.sortedByCoord.out.bam')
rnaid.append('18NG002044-ALD_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/18NG002044-ALDAligned.sortedByCoord.out.bam')

#double
rnaid.append('21NG002225-EEE_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/21NG002225-EEEAligned.sortedByCoord.out.bam')
rnaid.append('21NG002225b-EEE_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/21NG002225-EEEAligned.sortedByCoord.out.bam')
#double

#double
rnaid.append('6619NG001842-EEE_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6619NG001842-EEEAligned.sortedByCoord.out.bam')
rnaid.append('6619NG001842b-EEE_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6619NG001842-EEEAligned.sortedByCoord.out.bam')
#double

#double
rnaid.append('6620NG001457-FSP_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6620NG001457-FSPAligned.sortedByCoord.out.bam')
rnaid.append('6620NG001457b-FSP_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6620NG001457-FSPAligned.sortedByCoord.out.bam')
#double

#double
rnaid.append('6622NG000396-ITD_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6620NG001457-FSPAligned.sortedByCoord.out.bam')
rnaid.append('6622NG000396b-ITD_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6620NG001457-FSPAligned.sortedByCoord.out.bam')
#double

#double
rnaid.append('6622NG001445-FSP_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6622NG001445-FSPAligned.sortedByCoord.out.bam')
rnaid.append('6622NG001445b-FSP_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6622NG001445-FSPAligned.sortedByCoord.out.bam')
# double

#double
rnaid.append('6622NG001886-FHM_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6622NG001886-FHMAligned.sortedByCoord.out.bam')
rnaid.append('6622NG001886b-FHM_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6622NG001886-FHMAligned.sortedByCoord.out.bam')
#double

rnaid.append('6622NG003246_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6622NG003246Aligned.sortedByCoord.out.bam')
rnaid.append('6622NG003716_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6622NG003716Aligned.sortedByCoord.out.bam')

#double
rnaid.append('6623NG000155_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000155Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000155b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000155Aligned.sortedByCoord.out.bam')
#double

#double
rnaid.append('6623NG000221_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000221Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000221b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000221Aligned.sortedByCoord.out.bam')
#double

#double
rnaid.append('6623NG000604_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000604Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000604b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000604Aligned.sortedByCoord.out.bam')
#double

#double
rnaid.append('6623NG000823_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000823Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000823b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000823Aligned.sortedByCoord.out.bam')
#double

#double
rnaid.append('6623NG000990_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000990Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000990b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000990Aligned.sortedByCoord.out.bam')
#double

rnaid.append('6623NG001056_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG001056Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG001599_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG001599Aligned.sortedByCoord.out.bam')
rnaid.append('18NG001626_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/18NG001626Aligned.sortedByCoord.out.bam')



df['RNA_ID'] = rnaid
df['RNA_BAM_FILE'] = rna_bam_file
#df['DROP_GROUP'] = 'fraser,fraser_external,batch_0,group1,PAXgene'
df['DROP_GROUP'] = 'fraser2,outrider,batch_0,group1,PAXgene'
df['PAIRED_END'] = 'TRUE'
df['STRAND'] = 'REVERSE'
df['COUNT_MODE'] = 'Union'
df['COUNT_OVERLAPS'] = 'FALSE'
df['INDIVIDAL_ID'] = rnaid
df['DNA_VCF_FILE'] = '/media/jbogoin/Data1/References/RNA-seq/hg19/DROP/qc_vcf_1000G_hg19.vcf.gz'
df['DNA_ID'] = rnaid


df.to_csv('sample_annotation.tsv', sep='\t', index=False)


print('\n#####\nFichier sample_annotation.tsv genere !\n')
