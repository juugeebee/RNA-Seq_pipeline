import pandas
import numpy
import os


print('\n*** Outsingle *** \n')

dossier = './QC/salmon'

i = 0

final = pandas.DataFrame()

for patient in os.listdir(dossier) :

    print(patient)
    fichier = dossier + '/' + patient + '/quant.sf'


    # Importer le csv dans un dataframe
    brut = pandas.read_csv(fichier, delimiter ="\t", dtype=str, header=[0])

    if i == 0 :

        final['name_of_genes'] = brut['Name']
        final[patient] = brut['NumReads']

        i = i + 1

    else :
        
        final[patient] = brut['NumReads']


ensg_full = final['name_of_genes'].str.split(pat='ENSG', n=-1, expand=True, regex=None)

ensg = ensg_full[1].str.split(pat='|', n=-1, expand=True, regex=None)
final['name_of_genes'] = ensg[0]


final.to_csv(('./outsingle_dataset.csv'), sep='\t', index=False)

print('\n*** outsingle-dataset, job done ! *** \n')

print('\n Lancement Outsingle')

python ~/outsingle/inject_outliers_fzse_pysvdcc.py ./outsingle_dataset.csv
python ~/outsingle/fast_zscore_estimation.py ./outsingle_dataset-wo-f1-b-z6.00.txt
python ~/outsingle/optht_svd_zs.py ./outsingle_dataset-wo-f1-b-z6.00-fzse-zs.txt

print('\n outsingle_dataset-wo-f1-b-z6.00-fzse-zs-svd-optht-zs.csv OutSingle score file genere!'\n)