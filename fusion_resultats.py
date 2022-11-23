import pandas
from os import listdir
from os.path import isfile, join

print('\n*** Fusion resulats pertargetcoverage** \n')


fichiers = []

for f in listdir('.'):
    if '_pertargetcoverage' in f :
        fichiers.append(f)


tissus = ['Paxgene','Paxgene','fibroblaste','Paxgene','Paxgene','lymphocyte','Paxgene','Paxgene','lymphocyte','Paxgene','lymphocyte',
'Paxgene','lymphocyte','lymphocyte','Paxgene','fibroblaste']

print('Nbre de fichiers : {}'.format(len(fichiers)))
print('Nbre de tissus : {}'.format(len(tissus)))


pax = []
fibro = []
lympho = []



for i in range(len(tissus)):

    sp = fichiers[i].split('.')
    ech = sp[0]
    
    if tissus[i] == 'Paxgene':

        df = pandas.read_csv(fichiers[i],header = [0], sep='\t')
        del df['length']
        del df['name']
        del df['%gc']
        del df['mean_coverage'] 
        del df['normalized_coverage'] 
        del df['min_normalized_coverage'] 
        del df['max_normalized_coverage']
        del df['min_coverage']
        del df['max_coverage'] 
        del df['pct_0x']
        df.rename(columns = {'read_count': 'read_count_' + ech}, inplace = True)
        
        if not pax :
            pax.append(df)
        else :
            pax.append(df['read_count_' + ech])
    
    if tissus[i] == 'fibroblaste':
        df = pandas.read_csv(fichiers[i],header = [0], sep='\t')
        del df['length']
        del df['name']
        del df['%gc']
        del df['mean_coverage'] 
        del df['normalized_coverage'] 
        del df['min_normalized_coverage'] 
        del df['max_normalized_coverage']
        del df['min_coverage']
        del df['max_coverage'] 
        del df['pct_0x']
        df.rename(columns = {'read_count': 'read_count_' + ech}, inplace = True)

        if not fibro :
            fibro.append(df)
        else : 
            fibro.append(df['read_count_' + ech])

    if tissus[i] == 'lymphocyte':
        df = pandas.read_csv(fichiers[i],header = [0], sep='\t')
        del df['length']
        del df['name']
        del df['%gc']
        del df['mean_coverage'] 
        del df['normalized_coverage'] 
        del df['min_normalized_coverage'] 
        del df['max_normalized_coverage']
        del df['min_coverage']
        del df['max_coverage'] 
        del df['pct_0x']
        df.rename(columns = {'read_count': 'read_count_' + ech}, inplace = True)

        if not lympho :
            lympho.append(df)
        else :
            lympho.append(df['read_count_' + ech])


df_pax = pandas.concat(pax, axis=1)
df_pax.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first', inplace=True, ignore_index=False)
df_pax.to_csv('synthese_paxgene.csv', sep='\t', index=False)

df_fibro = pandas.concat(fibro, axis=1)
df_fibro.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first', inplace=True, ignore_index=False)
df_fibro.to_csv('synthese_fibroblaste.csv', sep='\t', index=False)

df_lympho = pandas.concat(lympho, axis=1)
df_lympho.drop_duplicates(subset=['chrom', 'start', 'end'], keep='first', inplace=True, ignore_index=False)
df_lympho.to_csv('synthese_lymphocyte.csv', sep='\t', index=False)


print('\nJob done!\n')