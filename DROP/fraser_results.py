import pandas
import os


print('\n*** Mise en forme et annotations fraser. *** \n')


# Importer le tsv dans un dataframe
brut = pandas.read_csv('results.tsv', delimiter ="\t", dtype=str)


# Supprimer les colonnes inutiles
del brut['padjustGene']
del brut ['pValueGene']
del brut['PAIRED_END']
del brut['DROP_GROUP']
del brut['INDIVIDAL_ID']
del brut['DNA_ID']
del brut ['isExternal']


# Supprimer les patients _ref
indexNames_ref = brut[brut["sampleID"].str.contains("_ref") == True].index
brut.drop(indexNames_ref, inplace=True)


# Supprimer les hgncSymbol HLA
indexNames_hla = brut[brut["hgncSymbol"].str.contains("HLA") == True].index
brut.drop(indexNames_hla, inplace=True)


# Creer une colonne IGV_coordinate
brut['IGV_coordinate'] = brut['seqnames'] + ':' + brut['start'] + '-' + brut['end']


# Recuperer le nom du run 
run_path = os.getcwd()
run_path_l = run_path.split('/')
run_name = run_path_l[8]


####################################################################################
## ANNOTATIONS


panelapp_f = '/media/jbogoin/Data1/Annotations_RNA-Seq/panelapp/6sep2023/full_panelapp.tsv'
omim_f = '/media/jbogoin/Data1/Annotations_RNA-Seq/OMIM/6sep2023/genemap2.txt'
loeuf_f = '/media/jbogoin/Data1/Annotations_RNA-Seq/LOEUF/supplement/loeuf_dataset_grch38.tsv.gz'


panelapp_df = pandas.read_csv(panelapp_f, delimiter ="\t", dtype=str)
print(panelapp_df.columns)
# 'Gene Symbol'


omim_df = pandas.read_csv(omim_f, delimiter ="\t", dtype=str)
print(omim_df.columns)
# 'Gene Symbols'


loeuf_df = pandas.read_csv(loeuf_f, delimiter ="\t", dtype=str, compression='gzip')
print(loeuf_df.columns)
# '#gene'


print(brut.columns)
# 'hgncSymbol'


####################################################################################


# Trier sur la colonne pvalue
brut = brut.sort_values(by=["pValue"])


# Creer le fichier de sortie
file_name = 'fraser_results_' + run_name + '.csv'
brut.to_csv(file_name, sep ="\t", header=True, index=False)


print('\n*** Mise en forme et annotations fraser. Job Done ! *** \n')
