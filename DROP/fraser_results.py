import pandas
import os


print('\n*** Mise en forme et annotations fraser. *** \n')


# Importer le tsv dans un dataframe
brut = pandas.read_csv('results_gene_all.tsv', delimiter ="\t", dtype=str)


# Supprimer les colonnes inutiles
del brut['padjustGene']
del brut ['pValueGene']
# del brut['PAIRED_END']
# del brut['DROP_GROUP']
# del brut['INDIVIDAL_ID']
# del brut['DNA_ID']
# del brut ['isExternal']


# Supprimer les patients _ref
indexNames_ref = brut[brut["sampleID"].str.contains("_ref") == True].index
brut.drop(indexNames_ref, inplace=True)


# Supprimer les hgncSymbol HLA
indexNames_hla = brut[brut["hgncSymbol"].str.contains("HLA") == True].index
brut.drop(indexNames_hla, inplace=True)


# Creer une colonne IGV_coordinate
brut['IGV_coordinates'] = brut['seqnames'] + ':' + brut['start'] + '-' + brut['end']


# Supprimer les lignes avec une pValue à 1
indexNames_ref = brut[brut["pValue"] == 1].index
brut.drop(indexNames_ref, inplace=True)


# Recuperer le nom du run 
run_path = os.getcwd()
run_path_l = run_path.split('/')
run_name = run_path_l[8]


# ####################################################################################
# ## ANNOTATIONS


# panelapp_f = '/media/jbogoin/Data1/Annotations_RNA-Seq/panelapp/6sep2023/full_panelapp.tsv'
# omim_f = '/media/jbogoin/Data1/Annotations_RNA-Seq/OMIM/6sep2023/genemap2.txt'
# loeuf_f = '/media/jbogoin/Data1/Annotations_RNA-Seq/LOEUF/supplement/loeuf_dataset_grch38.tsv.gz'


# ##### PANELAPP
# panelapp_df = pandas.read_csv(panelapp_f, delimiter ="\t", dtype=str)
# # Supprimer les infos hg19
# del panelapp_df['Position GRCh37 Start']
# del panelapp_df['Position GRCh37 End']
# del panelapp_df['EnsemblId(GRch37)']
# # print(panelapp_df.columns)


# ##### OMIM
# omim_df = pandas.read_csv(omim_f, delimiter ="\t", dtype=str)
# # Creer une colonne IGV_coordinate
# omim_df['IGV_coordinates'] = omim_df['Chromosome'] + ':' + omim_df['Genomic Position Start'] + '-' + omim_df['Genomic Position End']
# #print('\n')
# # print(omim_df.columns)


# ##### LOEUF
# loeuf_df = pandas.read_csv(loeuf_f, delimiter ="\t", dtype=str, compression='gzip')
# loeuf_df['IGV_coordinates'] = loeuf_df['chromosome'] + ':' + loeuf_df['start_position'] + '-' + loeuf_df['end_position']
# #print('\n')
# # print(loeuf_df.columns)


# #print('\n')
# # print(brut.columns)


# ####################################################################################
# ### MERGING ###

# #### OMIM
# # brut = brut.sort_values(by=['IGV_coordinates'])
# # print(brut['IGV_coordinates'])
# # omim_df = omim_df.sort_values(by=['IGV_coordinates'])
# # print(omim_df['IGV_coordinates'])

# # merge_omim = brut.merge(omim_df, how='inner',\
# #  left_on='IGV_coordinates', right_on='IGV_coordinates', \
# #  suffixes=('_fraser', '_omim'))

# # print(merge_omim)


# #### LOEUF
# brut = brut.sort_values(by=['hgncSymbol'])
# # data types to string 
# brut = brut.astype(str)
# print(brut['hgncSymbol']) 

# loeuf_df = loeuf_df.sort_values(by=['#gene'])
# loeuf_df = loeuf_df.astype(str) 
# print(loeuf_df['#gene'])

# merge_loeuf = brut.merge(loeuf_df, how='inner',\
#   left_on='hgncSymbol', right_on='#gene', \
#   suffixes=('_fraser', '_loeuf'))

# print(merge_loeuf)





####################################################################################
# Trier sur la colonne pvalue
brut = brut.sort_values(by=["pValue"])


####################################################################################
# Creer le fichier de sortie
file_name = 'fraser2_resultats_filtres_' + run_name + '.tsv'
brut.to_csv(file_name, sep ="\t", header=True, index=False, float_format=str, )


print('\n*** Mise en forme et annotations fraser. Job Done ! *** \n')
