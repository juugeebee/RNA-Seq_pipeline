#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys, os, argparse, pprint, subprocess


pp = pprint.PrettyPrinter(indent=4)

print('\n***PREPARATION DES ANNOTATIONS***\n ')


gencodeBasic = '/media/jbogoin/Data1/Annotations_RNA-Seq/gencode.v47.chr_patch_hapl_scaff.basic.annotation.gff3.gz'


subprocess.call("zcat " + gencodeBasic + " | awk '{if ($3==\"CDS\" || $3==\"three_prime_UTR\" || $3==\"five_prime_UTR\") print $0}' \
    | grep \"gene_type=protein_coding;\" > gencode.basic.CDS.UTR.prot_coding.gff3", shell="/bin/bash")

subprocess.call("zcat " + gencodeBasic + " | awk '{if ($3==\"gene\") print $0}' | grep \"gene_type=protein_coding;\" \
    > gencode.basic.gene.prot_coding.gff3", shell="/bin/bash")

subprocess.call("zcat " + gencodeBasic + " | awk '{if ($3==\"exon\") print $0}' | grep \"gene_type=protein_coding;\" \
    | grep \"Ensembl_canonical\" > gencode.basic.exon.prot_coding.canonical.gff3", shell="/bin/bash")


gencodeCDS_UTR_f = "gencode.basic.CDS.UTR.prot_coding.gff3"
gencodeGene_f = "gencode.basic.gene.prot_coding.gff3"
gencodeExons_f = "gencode.basic.exon.prot_coding.canonical.gff3"


prev_exon = 0
prev_start = 0
prev_end = 0
prev_geneID = "0"
    
    
with open(gencodeExons_f, 'r') as gff, open("gencode.basic.exon.intron.prot_coding.canonical.bed", 'w') as out1 :
    for line in gff:
        line = line.rstrip()
        if not line.startswith('#'):
            seqid = line.split('\t')[0]
            source = line.split('\t')[1]
            feature = line.split('\t')[2]
            exonStart = int(line.split('\t')[3])-1 # 1-based to 0-based for conversion to bed
            exonEnd = int(line.split('\t')[4])
            score = line.split('\t')[5]
            strand = line.split('\t')[6]
            phase = line.split('\t')[7]
            attributes = line.split('\t')[8]
            key = seqid + '-' + str(exonStart) + '-' + str(exonEnd)
            geneID, geneName, exonNumber = None, None, None
            for i in attributes.split(';'):
                if i.startswith('gene_id'): geneID = i.split('=')[1].split('.')[0]
                if i.startswith('gene_name'): geneName = i.split('=')[1]
                if i.startswith('exon_number'): exonNumber = int(i.split('=')[1])
            if geneID != prev_geneID:
                # print(">>>>>>>>>>>>>> NEW GENE <<<<<<<<<<<<<<<")
                intron = {}
                prev_geneID = geneID
                prev_exon = 0

            if exonNumber != prev_exon and geneID == prev_geneID:
                prev_exon = exonNumber
                """ STRAND + """
                if strand == "+":
                    if exonNumber == 1:
                        intron['start'] = exonEnd
                        intron['end'] = None
                    else:
                        intron['end'] = exonStart
                        intron['number'] = exonNumber
                        # print(geneID,geneName,strand,"intron",intron['number'],seqid, intron['start'],intron['end'])
                        out1.write(seqid + '\t' + str(intron['start']) + '\t' + str(intron['end']) + '\t' + geneID + ';' + geneName + ';intron' + str(intron['number']) + '\t' '1000' + '\t' + strand + '\n')
                        intron['start'] = exonEnd

                    # print(geneID,geneName,strand,"exon",exonNumber,key.split('-')[0],key.split('-')[1],key.split('-')[2])
                    out1.write(key.split('-')[0] + '\t' + key.split('-')[1] + '\t' + key.split('-')[2] + '\t' + geneID + ';' + geneName + ';exon' + str(exonNumber) + '\t' + '1000' + '\t' + strand + '\n')

                """ STRAND - """
                if strand == "-":
                    if exonNumber == 1:
                        intron['end'] = exonStart
                        intron['start'] = None
                    else:
                        intron['start'] = exonEnd
                        intron['number'] = exonNumber - 1
                        # print(geneID,geneName,strand,"intron",intron['number'],seqid, intron['start'],intron['end'])
                        out1.write(seqid + '\t' + str(intron['start']) + '\t' + str(intron['end']) + '\t' + geneID + ';' + geneName + ';intron' + str(intron['number']) + '\t' '1000' + '\t' + strand + '\n')
                        intron['end'] = exonStart

                    # print(geneID,geneName,strand,"exon",exonNumber,key.split('-')[0],key.split('-')[1],key.split('-')[2])
                    out1.write(key.split('-')[0] + '\t' + key.split('-')[1] + '\t' + key.split('-')[2] + '\t' + geneID + ';' + geneName + ';exon' + str(exonNumber) + '\t' + '1000' + '\t' + strand + '\n')


subprocess.call("grep \"intron\" gencode.basic.exon.intron.prot_coding.canonical.bed | sed 's/intron[0-9]*/intron/' > gencode.basic.intron.prot_coding.canonical.bed", shell="/bin/bash")


with open(gencodeGene_f, 'r') as gff, open("gencode.basic.gene.prot_coding.bed", 'w') as out2:
    for line in gff:
        line = line.rstrip()
        if not line.startswith('#'):
            seqid = line.split('\t')[0]
            source = line.split('\t')[1]
            feature = line.split('\t')[2]
            start = int(line.split('\t')[3])-1 # 1-based to 0-based for conversion to bed
            end = int(line.split('\t')[4])
            score = line.split('\t')[5]
            strand = line.split('\t')[6]
            phase = line.split('\t')[7]
            attributes = line.split('\t')[8]
            key = seqid + '-' + str(start) + '-' + str(end)
            geneID, geneName = None, None
            for i in attributes.split(';'):
                if i.startswith('gene_id'): geneID = i.split('=')[1].split('.')[0]
                if i.startswith('gene_name'): geneName = i.split('=')[1]
            
            out2.write(seqid + '\t' + str(start) + '\t' + str(end) + '\t' + geneID + ';' + geneName + ';' + feature + '\t' + '1000' + '\t' + strand + '\n')


with open(gencodeCDS_UTR_f, 'r') as gff, open("gencode.basic.CDS.UTR.prot_coding.canonical.bed", 'w') as out2:
    for line in gff:
        line = line.rstrip()
        if not line.startswith('#'):
            seqid = line.split('\t')[0]
            source = line.split('\t')[1]
            feature = line.split('\t')[2]
            start = int(line.split('\t')[3])-1 # 1-based to 0-based for conversion to bed
            end = int(line.split('\t')[4])
            score = line.split('\t')[5]
            strand = line.split('\t')[6]
            phase = line.split('\t')[7]
            attributes = line.split('\t')[8]
            key = seqid + '-' + str(start) + '-' + str(end)
            geneID, geneName = None, None
            for i in attributes.split(';'):
                if i.startswith('gene_id'): geneID = i.split('=')[1].split('.')[0]
                if i.startswith('gene_name'): geneName = i.split('=')[1]


            out2.write(seqid + '\t' + str(start) + '\t' + str(end) + '\t' + geneID + ';' + geneName + ';' + feature + '\t' + '1000' + '\t' + strand + '\n')


subprocess.call("cat gencode.basic.intron.prot_coding.canonical.bed gencode.basic.CDS.UTR.prot_coding.canonical.bed | sort -V > gencode.basic.CDS.UTR.intron.prot_coding.canonical.bed", shell="/bin/bash")


os.remove(gencodeCDS_UTR_f)
os.remove(gencodeGene_f)
os.remove(gencodeExons_f)
os.remove("gencode.basic.intron.prot_coding.canonical.bed")
os.remove("gencode.basic.CDS.UTR.prot_coding.canonical.bed")


subprocess.call("mkdir -p Fichiers_annotes", shell="/bin/bash")
subprocess.call("mv gencode.basic.CDS.UTR.intron.prot_coding.canonical.bed Fichiers_annotes", shell="/bin/bash")
subprocess.call("mv gencode.basic.exon.intron.prot_coding.canonical.bed Fichiers_annotes", shell="/bin/bash")
subprocess.call("mv gencode.basic.gene.prot_coding.bed Fichiers_annotes", shell="/bin/bash")


print('***PREPARATION DES ANNOTATIONS DONE !***\n ')

### EOF ###