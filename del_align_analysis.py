from Bio import AlignIO
from itertools import groupby
from operator import itemgetter
from sys import argv
import os
from Bio.SeqUtils import GC

script, align_folder = argv

folder = os.listdir(align_folder)

oup_del = open("del_statistics.txt","w")
oup_ins = open("ins_statistics.txt","w")
oup_GC = open("GCstatistics.txt","w")

del_list = []
ins_list = []
num_indels = []
GClist_pseudo = []
GClist_funct = []
identities = []

def percentID(seq1,seq2):
	x = 0
	for i,j in enumerate(seq1):
		if seq2[i] == j:
			x = x+1
		else:
			pass
	percID = 100*float(x)/len(seq1)
	
	return percID


for inl in folder:
	print inl
	
	os.system("trimal -in %s/%s -out temp.aln -gt 0.9"%(align_folder,inl)) #generate a trimmed alignment
	# to calculate GC content on conserved sections
	
	align_temp = AlignIO.read("temp.aln","clustal")
	

	align = AlignIO.read("%s/%s"%(align_folder,inl),"clustal")
	seq_IDs = []
	
	for seq in align:
		seq_IDs.append(seq.id)
	
	seqID = seq_IDs[-1]
	refID = seq_IDs[-2]
	
	dels = []
	ins = []
	
	for seq2 in align:
		if seq2.id == seqID:
			
			del_groups = [i for i,j in enumerate(seq2) if j == "-"]
			for k,g in groupby(enumerate(del_groups), lambda (index,item):(index - item)):
				group = map(itemgetter(1),g)
				dels.append(group)
			
			
			
		elif seq2.id == refID:
			
			
			ins_groups = [i for i,j in enumerate(seq2) if j == "-"]
			for k,g in groupby(enumerate(ins_groups), lambda (index,item):(index - item)):
				group = map(itemgetter(1),g)
				ins.append(group)
			#GCfunct = str(GC(seq.seq.replace("-","")))
			
			
		else:
			pass
	total_indels = len(dels) + len(ins)
	
	try:
		os.listdir("alignments_%d_indels"%total_indels)
		os.system("cp %s/%s alignments_%d_indels"%(align_folder,inl,total_indels))
	
	except OSError:
		os.system("mkdir alignments_%d_indels"%total_indels)
		os.system("cp %s/%s alignments_%d_indels"%(align_folder,inl,total_indels))
	
	num_indels.append(total_indels)
	
	
	
	
	seqA = ""
	seqB = ""
	for seq3 in align_temp:
		
		
		if seq3.id == seqID:
			seqA = seq3
			GClist_pseudo.append(str(GC(str(seq3.seq).replace("-",""))))
			
			
		elif seq3.id == refID:
			seqB = seq3
			GClist_funct.append(str(GC(str(seq3.seq).replace("-",""))))
			
		else:
			pass

	percID = percentID(str(seqA.seq),str(seqB.seq))
	identities.append(percID)
			
		
	print inl, num_indels, percID
		
	for x in dels:
		del_list.append(len(x))
	
	for i in ins:
		ins_list.append(len(i))
		
for x in del_list:
	oup_del.write("%d\n"%x)

for y in ins_list:
	oup_ins.write("%d\n"%y)

z = 0
oup_GC.write("pseudo\tfunctional\tindels\tperc id\n")

while z < len(GClist_pseudo):

	oup_GC.write("%s\t%s\t%s\t%s\n" %(GClist_pseudo[z],GClist_funct[z],num_indels[z],identities[z]))
	z = z +1


oup_ins.close()
oup_del.close()
oup_GC.close()

		
