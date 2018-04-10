from sys import argv

script, infile = argv

inl =open(infile,"rU")
oup = open("indel_by_gene.txt","w")

ID_list = []
temp = []
line = inl.readline()


while line != "":
	line = line.split("\t")
	ID = line[2].replace("\n","")
	indel = line[0]
	ID_list.append(ID)
	temp.append((ID,indel))
	line = inl.readline()
	
setIDs = set(ID_list)

oup.write("gene\tnum_indels\tnum_dels\tnum_triplets\n")

for gene in setIDs:
    indel_no3 = []
    indel3 = []
    
    for x in temp:
        if x[0] == gene:
            if int(x[1])%3 == 0:
                indel3.append(x[1])
            elif int(x[1])%3 != 0:
                indel_no3.append(x[1])
            else:
                pass
    num_indels = len(indel_no3) + len(indel3)
    num_triplets = len(indel3)
    num_dels = len([i for i in (indel_no3 + indel3) if int(i) < 0])
    
    oup.write("%s\t%d\t%d\t%d\n"%(gene,num_indels,num_dels,num_triplets))
    
oup.close()
            
            
