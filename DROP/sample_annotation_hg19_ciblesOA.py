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


rnaid.append('14OA000513_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/14OA000513Aligned.sortedByCoord.out.bam')
rnaid.append('20OA000508_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/20OA000508Aligned.sortedByCoord.out.bam')


#double
rnaid.append('20OA000549_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/20OA000549Aligned.sortedByCoord.out.bam')
rnaid.append('20OA000549b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/20OA000549Aligned.sortedByCoord.out.bam')
rnaid.append('20OA001880_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/20OA001880Aligned.sortedByCoord.out.bam')
rnaid.append('20OA001880b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/20OA001880Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002094_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6619OA002094Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002094b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6619OA002094Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002863_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6619OA002863Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002863b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6619OA002863Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002865_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6619OA002865Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002865b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6619OA002865Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002899_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6619OA002899Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002899b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6619OA002899Aligned.sortedByCoord.out.bam')
rnaid.append('6620OA003962_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6620OA003962Aligned.sortedByCoord.out.bam')
rnaid.append('6620OA003962b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6620OA003962Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA002468_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6621OA002468Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA002468b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6621OA002468Aligned.sortedByCoord.out.bam')
rnaid.append('66621OA003244_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6621OA003244Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA003244b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6621OA003244Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA003399_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6621OA003399Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA003399b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6621OA003399Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA003582_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6621OA003582Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA003582b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6621OA003582Aligned.sortedByCoord.out.bam')
rnaid.append('6622OA004101_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6622OA004101Aligned.sortedByCoord.out.bam')
rnaid.append('6622OA004101b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6622OA004101Aligned.sortedByCoord.out.bam')
rnaid.append('6623OA001275_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6623OA001275Aligned.sortedByCoord.out.bam')
rnaid.append('6623OA001275b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6623OA001275Aligned.sortedByCoord.out.bam')
rnaid.append('6623OA004806_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6623OA004806Aligned.sortedByCoord.out.bam')
rnaid.append('6623OA004806b_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_hg19_ciblesOA/PAXgene/6623OA004806Aligned.sortedByCoord.out.bam')
#double


df['RNA_ID'] = rnaid
df['RNA_BAM_FILE'] = rna_bam_file
#df['DROP_GROUP'] = 'fraser,fraser_external,batch_0,group1,PAXgene'
df['DROP_GROUP'] = 'fraser,outrider,batch_0,group1,PAXgene'
df['PAIRED_END'] = 'TRUE'
df['STRAND'] = 'REVERSE'
df['COUNT_MODE'] = 'Union'
df['COUNT_OVERLAPS'] = 'FALSE'
df['INDIVIDAL_ID'] = rnaid
df['DNA_VCF_FILE'] = '/media/jbogoin/Data1/References/RNA-seq/hg19/DROP/qc_vcf_1000G_hg19.vcf.gz'
df['DNA_ID'] = rnaid


df.to_csv('sample_annotation.tsv', sep='\t', index=False)


print('\n#####\nFichier sample_annotation.tsv genere !\n')
