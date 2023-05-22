
from os import listdir


id_list = []
nb_list = []

for f in listdir(".") :

    if "_STARpass1" in f :

        id = f.split("_")
        id_list.append(id[0])

        for g in listdir("./" + f) :

            if ".final.out" in g :

                with open("./" + f + "/" + g, "r") as log:

                    for ligne in log.readlines() :

                        if "Uniquely mapped reads number" in ligne :

                            read = ligne.split("|	")
                            nb_list.append(read[1])

print(len(id_list))
print(len(nb_list))


with open("uniquely_mapped_reads.txt", "w") as fichier:
    
    fichier.write("ID\tUniquely mapped reads number\n")
    
    for i in range(len(id_list)) :

        fichier.write(id_list[i] + "\t" + nb_list[i])




