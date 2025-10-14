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



# Get id and path for witness BAM file
rnaid.append('14OA000513_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/14OA000513Aligned.sortedByCoord.out.bam')
rnaid.append('20OA000508_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/20OA000508Aligned.sortedByCoord.out.bam')
rnaid.append('20OA000549_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/20OA000549Aligned.sortedByCoord.out.bam')
rnaid.append('20OA001880_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/20OA001880Aligned.sortedByCoord.out.bam')
rnaid.append('6618OA004357_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6618OA004357Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002094_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6619OA002094Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002863_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6619OA002863Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA002865_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6619OA002865Aligned.sortedByCoord.out.bam')
rnaid.append('6619OA003864_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6619OA003864Aligned.sortedByCoord.out.bam')
rnaid.append('6620OA002899_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6620OA002899Aligned.sortedByCoord.out.bam')
rnaid.append('6620OA003962_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6620OA003962Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA002468_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6621OA002468Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA003244_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6621OA003244Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA003399_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6621OA003399Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA003582_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6621OA003582Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA005066_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6621OA005066Aligned.sortedByCoord.out.bam')
rnaid.append('6621OA006190_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6621OA006190Aligned.sortedByCoord.out.bam')
rnaid.append('6622OA004101_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6622OA004101Aligned.sortedByCoord.out.bam')
rnaid.append('6623OA001275_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6623OA001275Aligned.sortedByCoord.out.bam')
rnaid.append('6623OA004806_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6623OA004806Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000124_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000124Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000127_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000127Aligned.sortedByCoord.out.bam ')
rnaid.append('6624CY000262_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000262Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000263_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000263Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000264_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000264Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000265_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000265Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000266_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000266Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000267_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000267Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000268_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000268Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000714_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000714Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000730_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000730Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000731_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000731Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000742_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000742Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000799_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000799Aligned.sortedByCoord.out.bam')
rnaid.append('6624CY000812_ref')
rna_bam_file.append('/media/jbogoin/Data1/DROP_BAM_OA_hg38/PAXgene/6624CY000812Aligned.sortedByCoord.out.bam')


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


print('\n#####\nFichier sample_annotation_OA.tsv genere !\n')
