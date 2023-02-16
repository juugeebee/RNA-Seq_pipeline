import pysam


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


samfile = pysam.AlignmentFile('6620NG001457-FSP.marked_duplicates.bam', 'rb')

sample = '6620NG001457-FSP'


print('Echantillon {}.\n'.format(sample))


with open(sample + '_cover_per_base.txt', 'w') as f:

    f.write("Chromosome\tPosition\tCouverture\tGene\n")

    nbre_region = 0 

    # for region in bed
    for i in range(len(start_l)+1):

        nbre_region = nbre_region + 1
        
        nbre_base = 0

        # for base in region 
        for base in range((int(start_l[i])-1),int(stop_l[i])):

            coverage_per_base = samfile.count_coverage(chrom_l[i], int(base), int(base+1), read_callback='all', quality_threshold=0)
            total_coverage_per_base = coverage_per_base[0][0] + coverage_per_base[1][0] + coverage_per_base[2][0] + coverage_per_base[3][0]

            f.write(str(chrom_l[i]) + ":" + str(base+1) + "\t" + str(total_coverage_per_base) + "\t" + str(gene_l[i]) + "\n")

            nbre_base = nbre_base + 1
        
        print("Region : {0}:{1}-{2}, {3}.".format(chrom_l[i],start_l[i],stop_l[i],gene_l[i]))
        print("Taille de la region : {}.".format(int(stop_l[i])-int(start_l[i])+1))
        print("Nombre de bases = {}.\n".format(nbre_base))

    
    print("Nombre de regions : {}.".format(nbre_region))


print('Echantillon {}: OK.'.format(sample))    


samfile.close()