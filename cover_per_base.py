import pysam


cibles = '/media/jbogoin/Data1/References/cibles_panels_NG/cibles_neuro_01032022.bed'

chrom_l = []
start_l = []
stop_l = []

with open(cibles, 'r') as filin:
    lignes = filin.readlines()
    for ligne in lignes:
        n = ligne.split("\n")
        t = n[0].split( "\t")
        chrom_l.append(t[0])
        start_l.append(t[1])
        stop_l.append(t[2])


print("\nnumber of intervals = {}".format(len(start_l)))


samfile = pysam.AlignmentFile('21NG002225-EEE.marked_duplicates.bam', 'rb')


pos_l = []
n_l = []
nbre_read = 0 


with open('cover_per_base.txt', 'w') as f:

    f.write("Positions\tCouvertures\n")

    for interval in range(len(start_l)):

        compteur = 0
        pos = start_l[interval]

        for read in samfile.fetch(chrom_l[interval], int(start_l[interval]), int(stop_l[interval])+1):

            nbre_read = nbre_read + 1 
        
            if read.pos == pos :
                compteur = compteur + 1

            else :
                pos_l.append(read.pos)
                n_l.append(compteur)
                f.write(str(read.pos) + "\t" + str(compteur) + "\n")
                compteur = 0
                pos = read.pos


print("number of reads = {}".format(nbre_read))
print("number of positions = {}".format(len(pos_l)))
print("number of cover counts = {}\n".format(len(n_l)))


samfile.close()