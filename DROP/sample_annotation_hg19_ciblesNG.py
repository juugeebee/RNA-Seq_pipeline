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
rnaid.append('21NG002225-EEE_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/21NG002225-EEEAligned.sortedByCoord.out.bam')
rnaid.append('6619NG001842-EEE_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6619NG001842-EEEAligned.sortedByCoord.out.bam')
rnaid.append('6620NG001457-FSP_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6620NG001457-FSPAligned.sortedByCoord.out.bam')
rnaid.append('6622NG000396-ITD_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6620NG001457-FSPAligned.sortedByCoord.out.bam')
rnaid.append('6622NG001445-FSP_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6622NG001445-FSPAligned.sortedByCoord.out.bam')
rnaid.append('6622NG001886-FHM_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6622NG001886-FHMAligned.sortedByCoord.out.bam')
rnaid.append('DFT-AMI-GAU-2086-008-emneg_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/DFT-AMI-GAU-2086-008-emnegAligned.sortedByCoord.out.bam')
rnaid.append('DFT-AMI-GAU-2086-008-empos_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/DFT-AMI-GAU-2086-008-emposAligned.sortedByCoord.out.bam')
rnaid.append('DFT-SAL-REB-787-002_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/DFT-SAL-REB-787-002Aligned.sortedByCoord.out.bam')
rnaid.append('DFT-SAL-WAT-792-001_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/DFT-SAL-WAT-792-001Aligned.sortedByCoord.out.bam')
rnaid.append('ITD-GRE-ROU-1293-002-emneg_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/ITD-GRE-ROU-1293-002-emnegAligned.sortedByCoord.out.bam')
rnaid.append('ITD-GRE-ROU-1293-002-empos_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/ITD-GRE-ROU-1293-002-emposAligned.sortedByCoord.out.bam')
rnaid.append('PED-SAL-KIC-1414-003_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/PED-SAL-KIC-1414-003Aligned.sortedByCoord.out.bam')
rnaid.append('6622NG003246_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6622NG003246Aligned.sortedByCoord.out.bam')
rnaid.append('6622NG003716_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6622NG003716Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000155_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000155Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000221_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000221Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000604_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000604Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000823_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000823Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG000990_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG000990Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG001056_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG001056Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG001599_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG001599Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG001605_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG001605Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG001606_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG001606Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG001883_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG001883Aligned.sortedByCoord.out.bam')
rnaid.append('6623NG001884_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/6623NG001884Aligned.sortedByCoord.out.bam')
rnaid.append('KIF1Aex35c-3749-28_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/KIF1Aex35c-3749-28Aligned.sortedByCoord.out.bam')
rnaid.append('KIF1Aex35WT_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesNG/PAXgene/KIF1Aex35WTAligned.sortedByCoord.out.bam')


df['RNA_ID'] = rnaid
df['RNA_BAM_FILE'] = rna_bam_file
#df['DROP_GROUP'] = 'fraser,fraser_external,batch_0,group1,PAXgene'
df['DROP_GROUP'] = 'fraser,batch_0,group1,PAXgene'
df['PAIRED_END'] = 'TRUE'
df['STRAND'] = 'REVERSE'
df['COUNT_MODE'] = 'Union'
df['COUNT_OVERLAPS'] = 'FALSE'
df['INDIVIDAL_ID'] = rnaid
df['DNA_VCF_FILE'] = '/media/jbogoin/Data1/References/RNA-seq/hg19/DROP/qc_vcf_1000G_hg19.vcf.gz'
df['DNA_ID'] = rnaid


df.to_csv('sample_annotation.tsv', sep='\t', index=False)


print('\n#####\nFichier sample_annotation.tsv genere !\n')
