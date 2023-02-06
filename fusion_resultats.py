import pandas
from os import listdir
from os.path import isfile, join

print('\n*** Fusion resulats pertargetcoverage** \n')


fichiers = []

for f in listdir('.'):
    if '_pertargetcoverage' in f :
        fichiers.append(f)


tissus_dic = {"6619NG001842": "Paxgene", "6622NG001886": "Paxgene", "DFT-AMI-GAU-2086-008" : "lymphocyte",
    "ITD-GRE-ROU-1293-002": "lymphocyte", "PED-SAL-KIC-1414-003": "lymphocyte", "DFT-SAL-REB-787-002": "fibroblaste",
    "DFT-SAL-WAT-792-001": "fibroblaste", "6622NG001445": "Paxgene", "6620NG001457": "Paxgene", "6622NG000396": "Paxgene",
    "21NG002225": "Paxgene", "18NG001723": "Paxgene", "18NG002044": "Paxgene", "17NG002384": "Paxgene"}

tissus = []

for fichier in fichiers:
    for cle, valeur in tissus_dic.items():
        if cle in fichier:
            tissus.append(valeur)



print('Nbre de fichiers : {}'.format(len(fichiers)))
print('Nbre de tissus : {}'.format(len(tissus)))


pax = []
fibro = []
lympho = []


for i in range(len(tissus)):

    sp = fichiers[i].split('.')
    ech = sp[0]
    
    if tissus[i] == 'Paxgene':

        df = pandas.read_csv(fichiers[i], header=[0], sep='\t')
        del df['length']
        del df['%gc']
        del df['read_count'] 
        del df['normalized_coverage'] 
        del df['min_normalized_coverage'] 
        del df['max_normalized_coverage']
        del df['min_coverage']
        del df['max_coverage'] 
        del df['pct_0x']
        df.rename(columns = {'mean_coverage': 'mean_coverage_' + ech}, inplace = True)
        
        if not pax :
            pax.append(df)
        else :
            pax.append(df['mean_coverage_' + ech])
    
    if tissus[i] == 'fibroblaste':
        df = pandas.read_csv(fichiers[i],header = [0], sep='\t')
        del df['length']
        del df['%gc']
        del df['read_count'] 
        del df['normalized_coverage'] 
        del df['min_normalized_coverage'] 
        del df['max_normalized_coverage']
        del df['min_coverage']
        del df['max_coverage'] 
        del df['pct_0x']
        df.rename(columns = {'mean_coverage': 'mean_coverage_' + ech}, inplace = True)

        if not fibro :
            fibro.append(df)
        else : 
            fibro.append(df['mean_coverage_' + ech])

    if tissus[i] == 'lymphocyte':
        df = pandas.read_csv(fichiers[i],header = [0], sep='\t')
        del df['length']
        del df['%gc']
        del df['read_count'] 
        del df['normalized_coverage'] 
        del df['min_normalized_coverage'] 
        del df['max_normalized_coverage']
        del df['min_coverage']
        del df['max_coverage'] 
        del df['pct_0x']
        df.rename(columns = {'mean_coverage': 'mean_coverage_' + ech}, inplace = True)

        if not lympho :
            lympho.append(df)
        else :
            lympho.append(df['mean_coverage_' + ech])


df_pax = pandas.concat(pax, axis=1)
df_pax.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first', inplace=True, ignore_index=False)

df_fibro = pandas.concat(fibro, axis=1)
df_fibro.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first', inplace=True, ignore_index=False)

df_lympho = pandas.concat(lympho, axis=1)
df_lympho.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first', inplace=True, ignore_index=False)


## Ajouter les panels dans une nouvelle colonne

df_pax['Panel'] = None
df_fibro['Panel'] = None
df_lympho['Panel'] = None


df_dips = pandas.read_csv('/media/jbogoin/Data1/References/cibles_panels_NG/DIPS_v4_cibles_20210728.bed',header = None, sep='\t')
df_dips['Panel'] = 'DIPS'
del df_dips[0]
del df_dips[1]
del df_dips[2]
del df_dips[4]
del df_dips[5]


df_ee = pandas.read_csv('/media/jbogoin/Data1/References/cibles_panels_NG/EE_v5_cibles_20220301.bed',header = None, sep='\t')
df_ee['Panel'] = 'EE'
del df_ee[0]
del df_ee[1]
del df_ee[2]


df_spatax = pandas.read_csv('/media/jbogoin/Data1/References/cibles_panels_NG/SPATAX_CMT_v1_cibles_20210329.bed',header = None, sep='\t')
df_spatax['Panel'] = 'SPATAX-CMT'
del df_spatax[0]
del df_spatax[1]
del df_spatax[2]


df_cibles = pandas.concat([df_dips, df_ee, df_spatax])


df_rnaseq = pandas.read_csv('/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions.bed',header = None, sep='\t')
del df_rnaseq[0]
del df_rnaseq[1]
del df_rnaseq[2]


df_panel = df_rnaseq.merge(df_cibles, left_on=3, right_on=3, how='left')
df_panel.to_csv('cibles.csv', sep='\t', index=False)



df_pax_f = df_pax.merge(df_panel, left_on='name', right_on=3, how='outer')
del df_pax_f['Panel_x']
del df_pax_f[3]
df_pax_f.rename(columns = {'Panel_y': 'Panel'}, inplace = True)
df_pax_f.sort_values(by=['name', 'start'], inplace=True)
df_pax_f.drop_duplicates(keep='first', inplace=True, ignore_index=True)


df_fibro_f = df_fibro.merge(df_panel, left_on='name', right_on=3, how='outer')
del df_fibro_f['Panel_x']
del df_fibro_f[3]
df_fibro_f.rename(columns = {'Panel_y': 'Panel'}, inplace = True)
df_fibro_f.sort_values(by=['name', 'start'], inplace=True)
df_fibro_f.drop_duplicates(keep='first', inplace=True, ignore_index=True)

df_lympho_f = df_lympho.merge(df_panel, left_on='name', right_on=3, how='outer')
del df_lympho_f['Panel_x']
del df_lympho_f[3]
df_lympho_f.rename(columns = {'Panel_y': 'Panel'}, inplace = True)
df_lympho_f.sort_values(by=['name', 'start'], inplace=True)
df_lympho_f.drop_duplicates(keep='first', inplace=True, ignore_index=True)


## Fichiers finaux

df_pax_f.to_csv('synthese_paxgene.csv', sep='\t', index=False)
df_fibro_f.to_csv('synthese_fibroblaste.csv', sep='\t', index=False)
df_lympho_f.to_csv('synthese_lymphocyte.csv', sep='\t', index=False)


print('\nJob done!\n')