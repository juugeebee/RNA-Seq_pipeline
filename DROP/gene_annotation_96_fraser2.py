#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys, os, pprint, re, gzip, subprocess, xlsxwriter, glob, argparse
from datetime import date
from subprocess import Popen, PIPE
import pandas


pp = pprint.PrettyPrinter(indent=4)


print('\n***GENE ANNOTATIONS FRASER2***\n ')


hgnc_f = "/media/jbogoin/Data1/Annotations_WES_pipeline.v2/HGNC/1mar2024/hgnc.tsv"
panelapp_eng_f = "/media/jbogoin/Data1/Annotations_WES_pipeline.v2/panelapp/5aug2025/panelapp_england.tsv"
panelapp_aus_f = '/media/jbogoin/Data1/Annotations_WES_pipeline.v2/panelapp/5aug2025/panelapp_australia.tsv'
clinvar_f = "/media/jbogoin/Data1/Annotations_WES_pipeline.v2/clinvar/1mar2024/gene_condition_source_id"
omim_f = "/media/jbogoin/Data1/Annotations_WES_pipeline.v2/OMIM/1mar2024/genemap2.txt"
hpo_f = "/media/jbogoin/Data1/Annotations_WES_pipeline.v2/HPO/1mar2024/genes_to_phenotype.txt"
loeuf_f = "/media/jbogoin/Data1/Annotations_WES_pipeline.v2/LOEUF/supplement/loeuf_dataset_grch38.tsv.gz"
gencodeGene_f = "./Fichiers_annotes/gencode.basic.gene.prot_coding.bed"


# Gencode v48
#fraser2_f = './drop/output/processed_results/aberrant_splicing/results/v48/fraser/fraser2/results.tsv'

#fraser2_f = './drop/output/processed_results/aberrant_splicing/results/v48/fraser/fraser2/results_gene_all.tsv'
fraser2_f = './drop/output/processed_results/aberrant_splicing/results/v48/fraser/fraser2/FRASER2_results_partie_1.tsv'
fraser2_f2 = './drop/output/processed_results/aberrant_splicing/results/v48/fraser/fraser2/FRASER2_results_partie_2.tsv'


def tabix_query(file_fn, chrom_fn, start_fn, end_fn):
    """Call tabix (1-based) and generate an array of strings for each line it returns."""
    query = '{}:{}-{}'.format(chrom_fn, start_fn, end_fn)
    process = Popen(['tabix', '-f', file_fn, query], stdout=PIPE)
    for line in process.stdout:
        yield line.strip().split()


def hgnc_parser(ensemblD, hgnc_f):
	with open(hgnc_f, 'r') as hgnc:
		for line in hgnc:
			line = line.strip('\n')
			if not line.startswith("HGNC ID\tApproved symbol"):
				hgnc_id, gene_symbol, alias_symbols, refseq_id, entrez_id, ensembl_id = "", "", "", "", "", ""
				items = line.split("\t")
				hgnc_id = items[0]
				gene_symbol = items[1]
				alias_symbols = items[2]
				entrez_id = items[5]
				omim_geneID = items[6]
				refseq_id = items[4]
				ensembl_id = items[8]

				""" ENSEMBL ID DICTIONARY """
				if ensembl_id != '':
					if ensembl_id not in ensemblD:
						ensemblD[ensembl_id] = {}
						ensemblD[ensembl_id]['hgnc_id'] = hgnc_id
						ensemblD[ensembl_id]['gene_symbol'] = gene_symbol
						ensemblD[ensembl_id]['entrez_id'] = entrez_id
						ensemblD[ensembl_id]['omim_geneID'] = omim_geneID
					else:
						if ensemblD[ensembl_id]['hgnc_id'] == '':
							ensemblD[ensembl_id]['hgnc_id'] = hgnc_id
						if ensemblD[ensembl_id]['gene_symbol'] == '':
							ensemblD[ensembl_id]['gene_symbol'] = gene_symbol
						if ensemblD[ensembl_id]['entrez_id'] == '':
							ensemblD[ensembl_id]['entrez_id'] = entrez_id
						if ensemblD[ensembl_id]['omim_geneID'] == '':
							ensemblD[ensembl_id]['omim_geneID'] = omim_geneID			
	return ensemblD


def panelapp_parser(ensembl_d, panelapp_f):
	with open(panelapp_f, 'r') as panelapp:
		if "panelapp_australia.tsv" in panelapp_f:
			which_panelapp = "panelapp_aus"
		elif "panelapp_england.tsv" in panelapp_f:
			which_panelapp = "panelapp_eng"
		for line in panelapp:
			line = line.strip('\n')
			if "Mode of Inheritance\tPanel Name\tConfidence" not in line:
				gene_symbol, ensembl_id, moi, level, color, panel, abv = "","","","","","",""

				gene_symbol = line.split("\t")[1]
				ensembl_id = line.split("\t")[2]
				hgnc_id = line.split("\t")[3]
				moi = line.split("\t")[4]
				panel = line.split("\t")[5]
				if panel == "COVID-19 research":
					continue
				level = line.split("\t")[6]
				evidence = line.split("\t")[7]

				if moi.startswith('MONOALLELIC, autosomal or pseudoautosomal'):
					if 'maternally imprinted' in moi:
						abv = 'IMP_patexp'
					elif 'paternally imprinted' in moi:
						abv = 'IMP_matexp'
					else:
						abv = 'AD'
				elif moi == 'BIALLELIC, autosomal or pseudoautosomal' or moi == 'BOTH monoallelic and biallelic, autosomal or pseudoautosomal':
					abv = 'AR'
				elif moi.startswith('BOTH monoallelic and biallelic'):
					abv = 'AD/AR'
				elif moi == 'X-LINKED: hemizygous mutation in males, biallelic mutations in females':
					abv = 'XLR'
				elif moi == 'X-LINKED: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)' or moi == 'X linked: hemizygous mutation in males, monoallelic mutations in females may cause disease (may be less severe, later onset than males)':
					abv = 'XL'
				elif moi == 'MITOCHONDRIAL':
					abv = 'mito'
				elif moi == '' or moi == 'Unknown' or moi.startswith('Othe'):
					abv = 'unk'
				else:
					abv = moi
			
				cat_hits = "Lvl" + level + "|" + abv + "|" + panel
				# print(gene_symbol,ensembl_id,cat_hits)

				""" ENSEMBL ID DICTIONARY """
				if ensembl_id != '':
					if ensembl_id not in ensembl_d:
						ensembl_d[ensembl_id] = {}
						ensembl_d[ensembl_id]['hgnc_id'] = hgnc_id
						ensembl_d[ensembl_id]['gene_symbol'] = gene_symbol
						ensembl_d[ensembl_id][which_panelapp] = cat_hits
					else:
						if ensembl_d[ensembl_id]['hgnc_id'] == '':
							ensembl_d[ensembl_id]['hgnc_id'] = hgnc_id
						if ensembl_d[ensembl_id]['gene_symbol'] == '':
							ensembl_d[ensembl_id]['gene_symbol'] = gene_symbol
						if ensembl_d[ensembl_id][which_panelapp] == '':
							ensembl_d[ensembl_id][which_panelapp] = cat_hits
						if ensembl_d[ensembl_id][which_panelapp] != '':
							if cat_hits not in ensembl_d[ensembl_id][which_panelapp]:
								ensembl_d[ensembl_id][which_panelapp] += "," + cat_hits
				else:
					# pas d'ensembl_id
					continue

	return ensemblD


def omim_parser(ensemblD, omim_f):
	with open(omim_f, 'r') as omim:
		for line in omim:
			line = line.strip('\n')
			if not line.startswith("#"):
				mimNumber, gene_symbol, entrez_id, ensembl_id, phenotypeString = "", "", "", "", "" 
				items = line.split('\t')
				mimNumber = items[5]
				gene_symbol = items[8]
				entrez_id = items[9]
				ensembl_id = items[10]
				phenotypeString = items[12]

				if not phenotypeString:
					continue

				for phenotype in phenotypeString.split(';'):
					phenotype = phenotype.strip()
					inherit = ''
					# print(">>> ",gene_symbol, phenotype)
					# Long phenotype
					matcher = re.match(r'^(.*),\s(\d{6})\s\((\d)\)(|, (.*))$', phenotype)
					if matcher:
						phenotype = matcher.group(1)
						phenotypeMimNumber = matcher.group(2)
						phenotypeMappingKey = matcher.group(3)
						inheritances = matcher.group(5)
						if inheritances:
							inherit = inheritances.replace("Autosomal dominant", "AD" )
							inherit = inherit.replace("Autosomal recessive", "AR" )
							inherit = inherit.replace("Somatic mutation", "SOM Mut" )
							inherit = inherit.replace("Somatic mosaicism", "SOM Mos" )
							inherit = inherit.replace("Digenic dominant", "DD" )
							inherit = inherit.replace("Digenic recessive", "DR" )
							inherit = inherit.replace("X-linked dominant", "XLD" )
							inherit = inherit.replace("X-linked recessive", "XLR" )
							inherit = inherit.replace("X-linked", "XL" )
							inherit = inherit.replace("Y-linked", "YL" )
							inherit = inherit.replace("Pseudoautosomal dominant", "PAD" )
							inherit = inherit.replace("Pseudoautosomal recessive", "PAR" )
							inherit = inherit.replace("Multifactorial", "MULTI" )

					else:
						# Short phenotype
						matcher = re.match(r'^(.*)\((\d)\)(|, (.*))$', phenotype)
						if matcher:
							phenotype = matcher.group(1)
							phenotypeMimNumber = ""
							phenotypeMappingKey = matcher.group(2)
							inheritances = matcher.group(3)

							if inheritances:
								inherit = inheritances.replace("Autosomal dominant", "AD" )
								inherit = inherit.replace("Autosomal recessive", "AR" )
								inherit = inherit.replace("Somatic mutation", "SOM Mut" )
								inherit = inherit.replace("Somatic mosaicism", "SOM Mos" )
								inherit = inherit.replace("Digenic dominant", "DD" )
								inherit = inherit.replace("Digenic recessive", "DR" )
								inherit = inherit.replace("X-linked dominant", "XLD" )
								inherit = inherit.replace("X-linked recessive", "XLR" )
								inherit = inherit.replace("X-linked", "XL" )
								inherit = inherit.replace("Y-linked", "YL" )
								inherit = inherit.replace("Pseudoautosomal dominant", "PAD" )
								inherit = inherit.replace("Pseudoautosomal recessive", "PAR" )
								inherit = inherit.replace("Multifactorial", "MULTI" )
					phenotype = phenotype.replace(", Autosomal dominant", "")
					phenotype = phenotype.replace(", Autosomal recessive", "")
					phenotype = phenotype.replace(", Somatic mutation", "")
					phenotype = phenotype.replace(", Somatic mosaicism", "")
					phenotype = phenotype.replace(", Digenic dominant", "")
					phenotype = phenotype.replace(", Digenic recessive", "")
					phenotype = phenotype.replace(", X-linked dominant", "")
					phenotype = phenotype.replace(", X-linked recessive", "")
					phenotype = phenotype.replace(", X-linked", "")
					phenotype = phenotype.replace(", Y-linked", "")
					phenotype = phenotype.replace(", Pseudoautosomal dominant", "")
					phenotype = phenotype.replace(", Pseudoautosomal recessive", "")
					phenotype = phenotype.replace(", Multifactorial", "")

					# omim = phenotype + " MOI:" + inherit + " OMIM:" + mimNumber
					if phenotypeMimNumber != "":
						omim = phenotype + " OMIM:" + phenotypeMimNumber
					else:
						omim = phenotype + " INHERIT:"

					""" ENSEMBL ID DICTIONARY """
					if ensembl_id != '':
						if ensembl_id not in ensemblD:
							ensemblD[ensembl_id] = {}
							ensemblD[ensembl_id]['omim_geneID'] = mimNumber
							ensemblD[ensembl_id]['omim_disease'] = omim
							ensemblD[ensembl_id]['omim_inheritance'] = inherit
							ensemblD[ensembl_id]['entrez_id'] = entrez_id
							ensemblD[ensembl_id]['gene_symbol'] = gene_symbol
						else:
							if ensemblD[ensembl_id]['gene_symbol'] == '':
								ensemblD[ensembl_id]['gene_symbol'] = gene_symbol
							if ensemblD[ensembl_id]['entrez_id'] == '':
								ensemblD[ensembl_id]['entrez_id'] = entrez_id
							if ensemblD[ensembl_id]['omim_geneID'] == '':
								ensemblD[ensembl_id]['omim_geneID'] = mimNumber
							if ensemblD[ensembl_id]['omim_disease'] == '':
								ensemblD[ensembl_id]['omim_disease'] = omim
							if ensemblD[ensembl_id]['omim_inheritance'] == '':
								ensemblD[ensembl_id]['omim_inheritance'] = inherit
	return ensemblD


def clinvar_parser(ensemblD, clinvar_f):
	with open(clinvar_f, 'r') as clinvar:
		curr_id = '00000'
		curr_gene = '00000'
		clinvar_l = []
		for line in clinvar:
			line = line.strip('\n')
			if not line.startswith("#GeneID"):
				entrez_id, gene_symbol, disease = "", "", ""
				items = line.split("\t")
				entrez_id = items[0]
				gene_symbol = items[1]
				if curr_id == '00000':
					curr_id = entrez_id
				if curr_gene == '00000':
					curr_gene = gene_symbol
				disease = items[4].replace(';','')
				omim_id = items[7]
				if omim_id != '':
					clinvar = disease + " OMIM:" + omim_id
				else:
					clinvar = disease
				# print(gene_symbol, entrez_id, clinvar)

				if entrez_id != curr_id:
					# print(">>> DIFFERENT GENE, write")
					### ensemblD ###
					for k in ensemblD:
						if curr_id == ensemblD[k]['entrez_id']:
							if ensemblD[k]['clinvar'] != '' and clinvar not in ensemblD[k]['clinvar']:
								ensemblD[k]['clinvar'] += ',' + ','.join(clinvar_l)
							else:
								ensemblD[k]['clinvar'] = ','.join(clinvar_l)
							break
					
					# print("<<< new gene and reset clinvar_l")
					clinvar_l = []
					curr_id = entrez_id
					curr_gene = gene_symbol
					clinvar_l.append(clinvar)
				else:
					if clinvar not in clinvar_l:
						# print(" --- add info to clinvar_l")
						clinvar_l.append(clinvar)
	return ensemblD


def hpo_parser(ensemblD, hpo_f):
	with open(hpo_f, 'r') as hpo:
		curr_id = '00000'
		curr_gene = '00000'
		hpo_l = []
		for line in hpo:
			line = line.strip('\n')
			if not line.startswith("#Format:"):
				entrez_id, gene_symbol, hpo_id, hpo_name, = "", "", "", ""
				items = line.split("\t")
				entrez_id = items[0]
				gene_symbol = items[1]
				if curr_id == '00000':
					curr_id = entrez_id
				if curr_gene == '00000':
					curr_gene = gene_symbol
				hpo_id = items[2]
				hpo_name = items[3]
				hpo = hpo_name + " " + hpo_id
				# print(gene_symbol, entrez_id, hpo)

				if entrez_id != curr_id:
					# print(">>> DIFFERENT GENE, write")
					### ensemblD ###
					for k in ensemblD:
						if curr_id == ensemblD[k]['entrez_id']:
							if ensemblD[k]['hpo'] != '' and (','.join(hpo_l) not in ensemblD[k]['hpo']):
								ensemblD[k]['hpo'] += "," + ','.join(hpo_l)
							else:
								ensemblD[k]['hpo'] = ','.join(hpo_l)
							break

					# print("<<< new gene and reset hpo_l")
					hpo_l = []
					curr_id = entrez_id
					curr_gene = gene_symbol
					hpo_l.append(hpo)
				else:
					if hpo not in hpo_l:
						# print(" --- add info to hpo_l")
						hpo_l.append(hpo)
	return ensemblD


def loeuf_parser(ensemblD, loeuf_f):
	with gzip.open(loeuf_f,'rb') as loeuf:
		for line in loeuf:
			line = line.decode("utf-8").strip("\n")
			if not line.startswith("#gene\ttranscript"):
				gene_symbol, ensembl_id, loeuf_score = "", "", ""
				gene_symbol = line.split("\t")[0]
				loeuf_score = line.split("\t")[30] # oe_lof_upper
				ensembl_id = line.split("\t")[64]

				if loeuf_score == "NA":
					loeuf_score = ""
					continue
				else:
					loeuf_score = format(float(loeuf_score), '.3f')

				""" ENSEMBL ID DICTIONARY """
				if ensembl_id not in ensemblD:
					ensemblD[ensembl_id] = {}
					ensemblD[ensembl_id]['loeuf'] = loeuf_score
					ensemblD[ensembl_id]['gene_symbol'] = gene_symbol
				else:
					if ensemblD[ensembl_id]['gene_symbol'] == '':
						ensemblD[ensembl_id]['gene_symbol'] = gene_symbol
					if ensemblD[ensembl_id]['loeuf'] == '':
						ensemblD[ensembl_id]['loeuf'] = loeuf_score
	return ensemblD


def complete_missing_keys(ensemblD):
	for i in ensemblD:
		if "gene_symbol" not in ensemblD[i]:
			ensemblD[i]['gene_symbol'] = ""
		if "hgnc_id" not in ensemblD[i]:
			ensemblD[i]['hgnc_id'] = ""
		if "entrez_id" not in ensemblD[i]:
			ensemblD[i]['entrez_id'] = ""
		if "panelapp_aus" not in ensemblD[i]:
			ensemblD[i]['panelapp_aus'] = ""
		if "panelapp_eng" not in ensemblD[i]:
			ensemblD[i]['panelapp_eng'] = ""		
		if "omim_geneID" not in ensemblD[i]:
			ensemblD[i]['omim_geneID'] = ""		
		if "omim_disease" not in ensemblD[i]:
			ensemblD[i]['omim_disease'] = ""
		if "omim_inheritance" not in ensemblD[i]:
			ensemblD[i]['omim_inheritance'] = ""
		if "clinvar" not in ensemblD[i]:
			ensemblD[i]['clinvar'] = ""
		if "hpo" not in ensemblD[i]:
			ensemblD[i]['hpo'] = ""
		if "loeuf" not in ensemblD[i]:
			ensemblD[i]['loeuf'] = ""
	return ensemblD


def replace_chars(string):
	new_string = string.replace(' - ', '-').replace(', ',',').replace(' ', '_')
	return new_string


def genes_coord(gencodeGene_f_fn):
	genesD_fn = {}
	with open(gencodeGene_f_fn, "r") as genes:
		for line in genes:
			line = line.rstrip()
			if ";gene" in line:
				chrom = line.split('\t')[0]
				start_0b = line.split('\t')[1] # 0-based
				end = line.split('\t')[2]
				info = line.split('\t')[3]
				geneID = info.split(';')[0]
				gene_symbol = info.split(';')[1]
				genesD_fn[geneID] = chrom + ":" + start_0b + "-" + end
	return genesD_fn 


def get_annotation_fraser2(fraser2_df, ensemblD_fn, genesD_fn):	
	for k in fraser2_df:
		if k.startswith('chr'):
			if "CDS" in fraser2_df[k]['feature_l']:
				fraser2_df[k]['mainFeature'] = "CDS"
			elif "UTR" in fraser2_df[k]['feature_l']:
				fraser2_df[k]['mainFeature'] =  "UTR"
			elif "intron" in fraser2_df[k]['feature_l']:
				fraser2_df[k]['mainFeature'] = "intron"
			else:
				fraser2_df[k]['mainFeature'] = "other"

			fraser2_df[k]['chrom'] = fraser2_df[k]['key'].split(':')[0]
			fraser2_df[k]['start_0b'] = int(fraser2_df[k]['key'].split(':')[1].split('-')[0]) # 0-based
			fraser2_df[k]['start_1b'] = int(fraser2_df[k]['key'].split(':')[1].split('-')[0]) + 1 # 1-based
			fraser2_df[k]['end'] = int(fraser2_df[k]['key'].split(':')[1].split('-')[1])
			fraser2_df[k]['size'] = (fraser2_df[k]['end'] - fraser2_df[k]['start_0b'])

			if len(fraser2_df[k]['geneID_l']) == 1:
				gene_chrom = genesD_fn[fraser2_df[k]['geneID_l'][0]].split(':')[0]
				gene_start_0b = int(genesD_fn[fraser2_df[k]['geneID_l'][0]].split(':')[1].split('-')[0]) # 0-based
				gene_end = int(genesD_fn[fraser2_df[k]['geneID_l'][0]].split(':')[1].split('-')[1])
				if fraser2_df[k]['chrom'] == gene_chrom and fraser2_df[k]['start_0b'] >= gene_start_0b and fraser2_df[k]['end'] <= gene_end:
					fraser2_df[k]['inGene'] = 1
				else:
					fraser2_df[k]['inGene'] = 0
			else:
				fraser2_df[k]['inGene'] = 0

			for i in range(0, len(fraser2_df[k]['geneID_l'])): # for each gene
				panelapp, hpo, clinvar, omim_inh, omim_dis, loeuf = '', '', '', '', '', ''
				if fraser2_df[k]['geneID_l'][i] in ensemblD_fn:
					panelapp_aus = ensemblD_fn[cnvD_fn['geneID_l'][i]]['panelapp_aus']
					panelapp_eng = ensemblD_fn[cnvD_fn['geneID_l'][i]]['panelapp_eng']
					hpo = ensemblD_fn[fraser2_df[k]['geneID_l'][i]]['hpo']
					clinvar = ensemblD_fn[fraser2_df[k]['geneID_l'][i]]['clinvar']
					omim_dis = ensemblD_fn[fraser2_df[k]['geneID_l'][i]]['omim_disease']
					omim_inh = ensemblD_fn[fraser2_df[k]['geneID_l'][i]]['omim_inheritance']
					loeuf = ensemblD_fn[fraser2_df[k]['geneID_l'][i]]['loeuf']
				fraser2_df[k]['panelapp_eng_l'].append(panelapp_eng)
				fraser2_df[k]['panelapp_aus_l'].append(panelapp_aus)
				fraser2_df[k]['hpo_l'].append(hpo)
				fraser2_df[k]['clinvar_l'].append(clinvar)
				fraser2_df[k]['omim_dis_l'].append(omim_dis)
				fraser2_df[k]['omim_inh_l'].append(omim_inh)
				fraser2_df[k]['loeuf_l'].append(loeuf)
	return fraser2_df


def fraser21TSV2BED():
	# importer le fichier brut
	fraser21_df = pandas.read_csv(fraser2_f, header=[0], sep='\t')

	# supprimer tous les Undetermined
	fraser21_df.drop(fraser21_df[fraser21_df['sampleID'] == 'Undetermined'].index, inplace = True)

	# Supprimer les patients _ref (normal)
	ik = fraser21_df[fraser21_df["sampleID"].str.contains("_ref") == True]
	print('Nombres de temoins Fraser2:')
	lik =ik['sampleID'].unique()
	print(len(lik))
	ik['sampleID'].to_csv('./Fichiers_annotes/temoins_fraser2.csv', header=False, index=False, sep='\t')

	indexNames_ref = fraser21_df[fraser21_df["sampleID"].str.contains("_ref") == True].index
	fraser21_df.drop(indexNames_ref, inplace=True)

	# Supprimer les hgncSymbol HLA
	indexNames_hla = fraser21_df[fraser21_df["hgncSymbol"].str.contains("HLA") == True].index
	fraser21_df.drop(indexNames_hla, inplace=True)

	# Supprimer les lignes avec une pValue à 1
	indexNames_pv = fraser21_df[fraser21_df["pValue"] == 1].index
	fraser21_df.drop(indexNames_pv, inplace=True)

	# Supprimer les colonnes inutiles
	del fraser21_df['pValueGene']
	###
	# del fraser21_df['PAIRED_END']
	# del fraser21_df['DROP_GROUP']
	# del fraser21_df['INDIVIDAL_ID']
	# del fraser21_df['DNA_ID']
	# del fraser21_df['isExternal']
	###

	# trier le df
	fraser21_df.sort_values(['sampleID','seqnames', 'start'], ascending=True, inplace=True)
	
	# creer un fichier bed 
	fraser2_bed = pandas.DataFrame(columns=['seqnames', 'start', 'end', 'hgncSymbol', 'sampleID', 'pValue', 'type'])
	fraser2_bed['seqnames'] = fraser21_df['seqnames']
	fraser2_bed['start'] = fraser21_df['start']
	fraser2_bed['end'] = fraser21_df['end']
	fraser2_bed['hgncSymbol'] = fraser21_df['hgncSymbol']
	fraser2_bed['sampleID'] = fraser21_df['sampleID']
	fraser2_bed['pValue'] = fraser21_df['pValue']
	fraser2_bed['type'] = fraser21_df['type']
	fraser2_bed.to_csv('./Fichiers_annotes/fraser2.bed', header=False, index=False, sep='\t')

	subprocess.call("bedtools intersect -a ./Fichiers_annotes/gencode.basic.CDS.UTR.intron.prot_coding.canonical.bed \
				 -b ./Fichiers_annotes/fraser2.bed -wb | sort -V > ./Fichiers_annotes/fraser2_bedtools.bed", shell="/bin/bash")
	return fraser21_df


def fraser22TSV2BED():
	# importer le fichier brut
	fraser22_df = pandas.read_csv(fraser2_f2, header=[0], sep='\t')

	# supprimer tous les Undetermined
	fraser22_df.drop(fraser22_df[fraser22_df['sampleID'] == 'Undetermined'].index, inplace = True)

	# Supprimer les patients _ref (normal)
	ik = fraser22_df[fraser22_df["sampleID"].str.contains("_ref") == True]
	print('Nombres de temoins Fraser2:')
	lik =ik['sampleID'].unique()
	print(len(lik))
	ik['sampleID'].to_csv('./Fichiers_annotes/temoins_fraser2.csv', header=False, index=False, sep='\t')

	indexNames_ref = fraser22_df[fraser22_df["sampleID"].str.contains("_ref") == True].index
	fraser22_df.drop(indexNames_ref, inplace=True)

	# Supprimer les hgncSymbol HLA
	indexNames_hla = fraser22_df[fraser22_df["hgncSymbol"].str.contains("HLA") == True].index
	fraser22_df.drop(indexNames_hla, inplace=True)

	# Supprimer les lignes avec une pValue à 1
	indexNames_pv = fraser22_df[fraser22_df["pValue"] == 1].index
	fraser22_df.drop(indexNames_pv, inplace=True)

	# Supprimer les colonnes inutiles
	del fraser22_df['pValueGene']
	###
	# del fraser22_df['PAIRED_END']
	# del fraser22_df['DROP_GROUP']
	# del fraser22_df['INDIVIDAL_ID']
	# del fraser22_df['DNA_ID']
	# del fraser22_df['isExternal']
	###

	# trier le df
	fraser22_df.sort_values(['sampleID','seqnames', 'start'], ascending=True, inplace=True)
	
	# creer un fichier bed 
	fraser2_bed = pandas.DataFrame(columns=['seqnames', 'start', 'end', 'hgncSymbol', 'sampleID', 'pValue', 'type'])
	fraser2_bed['seqnames'] = fraser22_df['seqnames']
	fraser2_bed['start'] = fraser22_df['start']
	fraser2_bed['end'] = fraser22_df['end']
	fraser2_bed['hgncSymbol'] = fraser22_df['hgncSymbol']
	fraser2_bed['sampleID'] = fraser22_df['sampleID']
	fraser2_bed['pValue'] = fraser22_df['pValue']
	fraser2_bed['type'] = fraser22_df['type']
	fraser2_bed.to_csv('./Fichiers_annotes/fraser2.bed', header=False, index=False, sep='\t')

	subprocess.call("bedtools intersect -a ./Fichiers_annotes/gencode.basic.CDS.UTR.intron.prot_coding.canonical.bed \
				 -b ./Fichiers_annotes/fraser2.bed -wb | sort -V > ./Fichiers_annotes/fraser2_bedtools.bed", shell="/bin/bash")
	return fraser22_df


def fraser2Annot(fraser2_df):
	with open('./Fichiers_annotes/fraser2_bedtools.bed','r') as bedtoolsSmoove:
		for line in bedtoolsSmoove:
			line = line.rstrip('\n')
			chrom = line.split('\t')[6]
			start_0b = line.split('\t')[7] # 0-based
			end = line.split('\t')[8]
			info = line.split('\t')[3]
			geneID = info.split(';')[0]
			gene_symbol = info.split(';')[1]
			feature = info.split(';')[2]
			strand = line.split('\t')[8]
			key_0b = chrom + ":" + start_0b + "-" + end # 0-based
			if key_0b in fraser2_df:
				if 'geneID_l' not in fraser2_df[key_0b]:
					fraser2_df[key_0b]['geneID_l'] = []
					fraser2_df[key_0b]['geneID_l'].append(geneID)
				else:
					if geneID not in fraser2_df[key_0b]['geneID_l']:
						fraser2_df[key_0b]['geneID_l'].append(geneID)
				if 'gene_symbol_l' not in fraser2_df[key_0b]:
					fraser2_df[key_0b]['gene_symbol_l'] = []
					fraser2_df[key_0b]['gene_symbol_l'].append(gene_symbol)
				else:
					if gene_symbol not in fraser2_df[key_0b]['gene_symbol_l']:
						fraser2_df[key_0b]['gene_symbol_l'].append(gene_symbol)
				if 'feature_l' not in fraser2_df[key_0b]:
					fraser2_df[key_0b]['feature_l'] = []
					fraser2_df[key_0b]['feature_l'].append(feature)
				else:
					if feature not in fraser2_df[key_0b]['feature_l']:
						fraser2_df[key_0b]['feature_l'].append(feature)

	fraser2_df = get_annotation_fraser2(fraser2_df, ensemblD, genesD)
	return fraser2_df


###### MAIN ######


print(">>> Loading annotations...\n")

print("> Genes coordinates...")
genesD = genes_coord(gencodeGene_f)
df_genes = pandas.DataFrame.from_dict(genesD, orient='index', columns=['coord'])
# df_genes['ENSG'] = df_genes.index.values.tolist()
df_genes.reset_index(inplace=True)
df_genes = df_genes.rename(columns={'index': 'ensg'})
# print(df_genes.columns.values.tolist())


print("\n>> ENSEMBL")
print("	> HGNC...")
ensemblD = {}
ensemblD = hgnc_parser(ensemblD, hgnc_f)
ensemblD = complete_missing_keys(ensemblD)

print("	> PANELAPP ENGLAND...")
ensemblD = panelapp_parser(ensemblD, panelapp_eng_f)
ensemblD = complete_missing_keys(ensemblD)

print("	> PANELAPP AUSTRALIA...")
ensemblD = panelapp_parser(ensemblD, panelapp_aus_f)
ensemblD = complete_missing_keys(ensemblD)

print("	> OMIM...")
ensemblD = omim_parser(ensemblD, omim_f)
ensemblD = complete_missing_keys(ensemblD)

print("	> CLINVAR...")
ensemblD = clinvar_parser(ensemblD, clinvar_f)
ensemblD = complete_missing_keys(ensemblD)

print("	> HPO...")
ensemblD = hpo_parser(ensemblD, hpo_f)
ensemblD = complete_missing_keys(ensemblD)

print("	> LOEUF...")
ensemblD = loeuf_parser(ensemblD, loeuf_f)
ensemblD = complete_missing_keys(ensemblD)

if '' in ensemblD: 
    del ensemblD['']


df_ensembl = pandas.DataFrame.from_dict(ensemblD, orient='index')
df_ensembl.reset_index(inplace=True)
df_ensembl = df_ensembl.rename(columns={'index': 'ensg'})
# print(df_ensembl.columns.values.tolist())
df_ensembl.to_csv('./Fichiers_annotes/ensembl.csv', sep='\t', index=False)


#### FRASER2 ####
print("\n> FRASER2 TSV TO DF...")
fraser2D1 = fraser21TSV2BED()
fraser2D2 = fraser22TSV2BED()


print("> FRASER2 Annotation...")
fraser2D1 = fraser2Annot(fraser2D1)
fraser2D2 = fraser2Annot(fraser2D2)

fraser2D1['coord'] = fraser2D1['seqnames'].astype(str) + ':' + fraser2D1['start'].astype(str) + '-' + fraser2D1['end'].astype(str)
# print(fraser2D1.columns.values.tolist())

fraser2D2['coord'] = fraser2D2['seqnames'].astype(str) + ':' + fraser2D2['start'].astype(str) + '-' + fraser2D2['end'].astype(str)
# print(fraser2D2.columns.values.tolist())


## MERGE ##
df_db = df_ensembl.merge(df_genes, left_on='ensg', right_on='ensg', suffixes=('_ensembl', '_genes'), how='outer')
df_db.to_csv('./Fichiers_annotes/db.csv', sep='\t', index=False)


df_final_fraser21 = fraser2D1.merge(df_db, left_on='hgncSymbol', right_on='gene_symbol', suffixes=('_fraser2', '_db'), how='left')
df_final_fraser22 = fraser2D2.merge(df_db, left_on='hgncSymbol', right_on='gene_symbol', suffixes=('_fraser2', '_db'), how='left')

## EXCEL ##
path = os.getcwd()
path_l = path.split('/')
run_name = path_l[-1]


# trier le df_final_fraser2
df_final_fraser21.sort_values(['pValue'], ascending=True, inplace=True)
del df_final_fraser21['coord_db']
# del df_final_fraser21['distNearestGene']

df_final_fraser22.sort_values(['pValue'], ascending=True, inplace=True)
del df_final_fraser22['coord_db']
# del df_final_fraser22['distNearestGene']


# ordonner les colonnes fraser2
cols_fraser2 = ['sampleID', 'panelapp_eng', 'panelapp_aus', 'omim_disease', 'omim_inheritance', 'hgncSymbol', 'pValue', 'padjustGene', 'psiValue',  
 'deltaPsi', 'hpo',  'clinvar', 'loeuf', 'seqnames', 'start', 'end', 'coord_fraser2', 'width', 'strand', 'gene_symbol', 'ensg', 'omim_geneID', 
 'hgnc_id', 'entrez_id', 'UTR_overlap', 'blacklist', 'type', 'potentialImpact', 'annotatedJunction', 'causesFrameshift', 
 'counts', 'totalCounts', 'meanCounts', 'meanTotalCounts', 'nonsplitCounts', 'nonsplitProportion', 'nonsplitProportion_99quantile']
df_final_fraser21 = df_final_fraser21.reindex(cols_fraser2, axis=1)
df_final_fraser22 = df_final_fraser22.reindex(cols_fraser2, axis=1)


# FICHIER FINAL FRASER2 ##
output_file = './Fichiers_annotes/FRASER2_' + run_name + '_annote.xlsx'

with pandas.ExcelWriter(output_file, engine='xlsxwriter') as writer2:
    max_rows = 1048575  # Excel limit - header
    workbook2 = writer2.book
    # formats (créés une seule fois)
    greenFormat2 = workbook2.add_format({'bg_color': 'lime'})
    redFormat2 = workbook2.add_format({'bg_color': 'red'})
    yellowFormat2 = workbook2.add_format({'bg_color': 'yellow'})

    dataframes = [df_final_fraser21, df_final_fraser22]  # à adapter aux noms réels de tes 2 df

    for i, df in enumerate(dataframes, start=1):
        sheet_name = f"FRASER2_{i}"
        print(f"Writing {sheet_name}...")
        df.to_excel(writer2, sheet_name=sheet_name, index=False)
        worksheet2 = writer2.sheets[sheet_name]
        longueur = len(df) + 1  # +1 header
        for col in ['B', 'C']:
            for lvl, fmt in [('Lvl3', greenFormat2), ('Lvl2', redFormat2), ('Lvl1', yellowFormat2)]:
                worksheet2.conditional_format(
                    f'{col}2:{col}{longueur}',
                    {'type': 'text', 'criteria': 'containing', 'value': lvl, 'format': fmt}
                )


print('\n***GENE ANNOTATIONS FRASER2 DONE !***\n ')


### EOF ###