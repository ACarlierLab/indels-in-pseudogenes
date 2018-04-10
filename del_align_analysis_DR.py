from Bio import AlignIO
from itertools import groupby
from operator import itemgetter
from sys import argv
import os
from Bio.SeqUtils import GC

script, align_folder = argv

folder = os.listdir(align_folder)

oup_del = open("del_statistics_DR.txt","w")
#oup_ins = open("ins_statistics.txt","w")
#oup_GC = open("GCstatistics.txt","w")
oup_DR = open("DR_statistics.txt","w")

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

def distance(a, b):
    #This function returns a score of 0 if 2 sequences have no matching letters, or a score equal to the matching letters (overlap length). 
    #If two sequences are completely identical, this score is equal to the length of the sequences
    return sum(map(lambda (x, y): 1 if x == y else 0, zip(a, b)))

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
    num_indels = 0
    
    for seq2 in align:
        if seq2.id == seqID:
            pseudo_seq = seq2
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

            print "ins: ",ins
            
            
            
        else:
            pass
        
    num_indels = int(len(dels) + len(ins))
    print num_indels
    
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
            
        
    #print inl, num_indels, percID
        
    for x in dels:
        print num_indels
        del_list.append((len(x),num_indels,inl))
    
    for i in ins:
        ins_list.append((len(i),num_indels,inl))
        #extracting the sequence of the insertion
        ins_seq = ""
        start = min(i)
        end = max(i)
        for k in i:
                   
            ins_seq = ins_seq + pseudo_seq[k]
        #extracting flanking regions, allowing for 1 mismatch (if n>3) and +/- 1 bp shift
        flankseqs=[]
        if (start > 0) and (end<len(seq2)):

            flankseqs.append(seq2[start-len(ins_seq):start].seq)
            flankseqs.append(seq2[start-len(ins_seq)+1:start+1].seq)
            flankseqs.append(seq2[start-1-len(ins_seq):start-1].seq)
            flankseqs.append(seq2[end:end+len(ins_seq)].seq)
            flankseqs.append(seq2[end+1:end+len(ins_seq) + 1].seq)
            flankseqs.append(seq2[end+2:end+len(ins_seq)+2].seq)
        
        else:
            pass
        #comparing insertion sequence to flanking regions of same length
        match = ("","")
        for flank in flankseqs:
            print "flank",flank
            
            if 1 < len(ins_seq) <= 3: #if insertion is smaller than 3, no mismatches are allowed
                if distance(str(flank),ins_seq) == len(ins_seq):
                    match = (ins_seq,flank)
                    print match
                else:
                    pass
            elif len(ins_seq) > 3:
                if distance(ins_seq,str(flank)) >= (len(ins_seq) -1): 
                    print distance, distance(ins_seq,str(flank)), len(ins_seq)
                    match = (ins_seq,str(flank))
                    print  "match",match
                else:
                    pass
            else:
                pass
            print "ins_seq",ins_seq,len(ins_seq)
        oup_DR.write("%s\t%s\t%s\t%d\n"%(inl,match[0],match[1],len(match[0])))
                
        
for x in del_list:
    oup_del.write("-%d\t%d\t%s\n"%(x[0],x[1],x[2]))

for y in ins_list:
    #oup_ins.write("%d\n"%y)
    oup_del.write("%d\t%d\t%s\n"%(y[0],y[1],y[2]))

z = 0
#oup_GC.write("pseudo\tfunctional\tindels\tperc id\n")

#while z < len(GClist_pseudo):

    #oup_GC.write("%s\t%s\t%s\t%s\n" %(GClist_pseudo[z],GClist_funct[z],num_indels[z],identities[z]))
    #z = z +1


#oup_ins.close()
oup_del.close()
#oup_GC.close()

        
