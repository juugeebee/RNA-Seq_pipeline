import pysam
import os
import array


print('\n#### COVER PER TARGET ###\n')


### Fichier CIBLES ###
######################

cibles = '/media/jbogoin/Data1/References/cibles_panels_NG/RNAseq_UFNeuro_v1_Regions.bed'


chrom_l = []
start_l = []
stop_l = []
gene_l = []


with open(cibles, 'r') as filin:
    lignes = filin.readlines()
    for ligne in lignes:
        n = ligne.split("\n")
        t = n[0].split( "\t")
        chrom_l.append(t[0])
        start_l.append(t[1])
        stop_l.append(t[2])
        gene_l.append(t[3])

print("\nNombre de regions cibles = {}.\n".format(len(start_l)))



### Fichiers BAM ###
####################

files = os.listdir(".")

for filename in files:


    if ('.marked_duplicates' in filename) and ('.bai' not in filename):
        
        name = filename.split('.marked_duplicates')
        sample = name[0]


        samfile = pysam.AlignmentFile(filename, 'rb')


        print('\nEchantillon : {}.'.format(sample))


        with open(sample + '_cover_per_target.txt', 'w') as f:

            f.write("Chromosome\tStart\tStop\tCouverture\tGene\n")

            nbre_region = 0

            # for region in bed
            for i in range(len(start_l)):

                nbre_region = nbre_region + 1

                region_coverage = 0
                nbre_base = 0

                # for base in region 
                for base in range((int(start_l[i])-1),int(stop_l[i])):

                    coverage_per_base = samfile.count_coverage(chrom_l[i],int(base),int(base+1), read_callback='all', quality_threshold=0)
                    total_coverage_per_base = coverage_per_base[0][0] + coverage_per_base[1][0] + coverage_per_base[2][0] + coverage_per_base[3][0]

                    region_coverage = region_coverage + total_coverage_per_base

                    nbre_base = nbre_base + 1

                mean_coverage = int(region_coverage/nbre_base)

                f.write(str(chrom_l[i]) + "\t" + str(start_l[i]) + "\t" + str(stop_l[i]) + "\t" + str(mean_coverage) + "\t" + str(gene_l[i]) + "\n")
                

                # print("Region : {0}:{1}-{2}, {3}.".format(chrom_l[i],start_l[i],stop_l[i],gene_l[i]))
                # print("Taille de la region : {}.".format(int(stop_l[i])-int(start_l[i])+1))
                # print("Nombre de bases = {}.\n".format(nbre_base))


            print("Nombre de regions : {}.".format(nbre_region))


        print('Echantillon : {} : OK.\n'.format(sample))


        samfile.close()


print('\n#### JOB DONE ! ###\n')




# print(dir(samfile))
# ['__class__', '__delattr__', '__dir__', '__doc__', '__enter__', '__eq__', '__exit__', '__format__', '__ge__', '__getattribute__', \
# '__gt__', '__hash__', '__init__', '__init_subclass__', '__iter__', '__le__', '__lt__', '__ne__', '__new__', '__next__', '__pyx_vtable__',\
#  '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__', '_open', \
#  'add_hts_options', 'category', 'check_index', 'check_truncation', 'close', 'closed', 'compression', 'count', \
#  'count_coverage', 'description', 'duplicate_filehandle', 'fetch', 'filename', 'find_introns', 'find_introns_slow', \
#  'format', 'get_index_statistics', 'get_reference_length', 'get_reference_name', 'get_tid', 'getrname', 'gettid', 'has_index', 'head', \
#  'header', 'index_filename', 'is_bam', 'is_bcf', 'is_closed', 'is_cram', 'is_open', 'is_read', 'is_remote', 'is_sam', 'is_stream', \
#  'is_valid_reference_name', 'is_valid_tid', 'is_vcf', 'is_write', 'lengths', 'mapped', 'mate', 'mode', 'nocoordinate', 'nreferences', \
#  'parse_region', 'pileup', 'reference_filename', 'references', 'reset', 'seek', 'tell', 'text', 'threads', 'unmapped', 'version', 'write']


# print(dir(read))
# ['__class__', '__copy__', '__deepcopy__', '__delattr__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__',\
#  '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__', '__pyx_vtable__', '__reduce__', \
#  '__reduce_ex__', '__repr__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__', 'aend', 'alen', 'aligned_pairs',\
#   'bin', 'blocks', 'cigar', 'cigarstring', 'cigartuples', 'compare', 'flag', 'from_dict', 'fromstring', 'get_aligned_pairs', 'get_blocks',\
#    'get_cigar_stats', 'get_forward_qualities', 'get_forward_sequence', 'get_overlap', 'get_reference_positions', 'get_reference_sequence', \
#    'get_tag', 'get_tags', 'has_tag', 'header', 'infer_query_length', 'infer_read_length', 'inferred_length', 'is_duplicate', 'is_paired',\
#     'is_proper_pair', 'is_qcfail', 'is_read1', 'is_read2', 'is_reverse', 'is_secondary', 'is_supplementary', 'is_unmapped', 'isize', \
#     'mapping_quality', 'mapq', 'mate_is_reverse', 'mate_is_unmapped', 'mpos', 'mrnm', 'next_reference_id', 'next_reference_name', \
#     'next_reference_start', 'opt', 'overlap', 'pnext', 'pos', 'positions', 'qend', 'qlen', 'qname', 'qqual', 'qstart', 'qual', 'query', \
#     'query_alignment_end', 'query_alignment_length', 'query_alignment_qualities', 'query_alignment_sequence', 'query_alignment_start', \
#     'query_length', 'query_name', 'query_qualities', 'query_sequence', 'reference_end', 'reference_id', 'reference_length', 'reference_name',\
#      'reference_start', 'rlen', 'rname', 'rnext', 'seq', 'setTag', 'set_tag', 'set_tags', 'tags', 'template_length', 'tid', 'tlen', 'to_dict',\
#       'to_string', 'tostring']