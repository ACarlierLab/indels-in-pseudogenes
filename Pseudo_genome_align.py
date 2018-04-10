from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
import StringIO,subprocess,os
from Bio.Blast import NCBIXML,NCBIWWW
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio import Data
import datetime
import os
import subprocess
from sys import argv
from Bio import AlignIO
from Bio.SeqUtils import GC

script,genome_folder = argv

class TranslationError(ValueError):
    pass


def extractprots(inlgenome):
    
    genome = SeqIO.read(open(inlgenome,"rU"),"genbank")

    protlist = []

    for feat in genome.features:
        try:
            pseudo = feat.qualifiers["pseudo"]
            pass
        except KeyError:
            try:
                seq = Seq(feat.qualifiers["translation"][0],generic_protein)
                name = feat.qualifiers["locus_tag"][0]
                try:
                    protID = feat.qualifiers["protein_id"][0]
                except KeyError:
                    protID = ""
                try:
                    product = feat.qualifiers["product"][0]
                except KeyError:
                    product = ""
                seqrec = SeqRecord(seq,id=name,description=product)
                protlist.append(seqrec)
                print seqrec
                
            except KeyError:
                pass
    return protlist


def intergenic(genome):
    loclist = []
    records = []
    for feat in genome.features:
        try: 
            feat.qualifiers["pseudo"]
            pass
        except KeyError:
            if feat.type in ["CDS","rRNA","tRNA"]:
                start = feat.location.nofuzzy_start
                end = feat.location.nofuzzy_end
                loclist.append((start,end))
            else:
                pass
        
    for i,item in enumerate(loclist[1:]):
        start = int(item[0])
        before_end = int(loclist[i][1])
        if start - before_end > 100:
            iter_seq = genome[before_end +1:start-1].seq
            name = "%d..%d"%(before_end +1, start)
            #rejecting completely masked intergenic regions. The 50nt threshold is completely arbitrary.
            #Necessary measure so blast doesn't crash on "all N" queries
            if iter_seq.count("N") < (len(iter_seq) - 50): 
                records.append(SeqRecord(iter_seq,id = name))
            else:
                pass
            
        else:
            pass
    
    return records

def run_blast(seqrecords,prot_fasta):

    #returns a list of tuples with query_ID, subject_ID, hit_start, hit_end
    #create blast database with prot_fasta
    
    if os.path.isfile("%s_db.pin"%prot_fasta) == False:
        print "formatting blast database\n"

        os.system("makeblastdb -in %s -dbtype prot -out %s_db" %(prot_fasta,prot_fasta))
    else:
        pass
    
    
    
    #run blastx 
    prot_list = []
    print "running blastx with %d sequences\n" %(len(seqrecords))
    i = 0
    pair_IDs = []
    for seq in seqrecords:

        if i in range(100,10000,100):
            print "blasting sequence #%d" %(i)
        else:
            pass
        i = i + 1
        query_length = len(seq.seq) 
        SeqIO.write(seq, "temp_blastx_query.fasta", "fasta")
        
        os.system("/home/aurelien/blast-2.2.25/bin/blastall -p blastx -F \"m S\"\
        -i temp_blastx_query.fasta -d %s_db -m 7 -a 8 -b 50 -e 0.000001 -o blastx_out.xml"%prot_fasta)
        blast_records = NCBIXML.parse(open("blastx_out.xml","rU"))
        
        for blast_record in blast_records:
            if blast_record.alignments != []:
                
                hit_coord = []
                for alignment in blast_record.alignments:

                    hsp_length = 0
                    hsp_coord = [] #contains a list of tuples with hit start and end coordinates
                    for hsp in alignment.hsps:
                        identities = 0
                        hsp_length = float((hsp.query_end - hsp.query_start)) / 3 #length of HSP in AA
                        
                        if hsp_length > 29: #keep only hsps larger than 30AA
                            identities = (float(hsp.identities)/hsp_length)
                            if hsp_coord != []:
                                if identities > 0.4: #keep only alignments with over 40% id
                                    if abs(hsp.query_start - max(hsp_coord)) < 50:
                                        hsp_coord.append(hsp.query_start)
                                        hsp_coord.append(hsp.query_end)
                                    elif abs(min(hsp_coord) - hsp.query_end) < 50:
                                        hsp_coord.append(hsp.query_start)
                                        hsp_coord.append(hsp.query_end)
                                    else:
                                        pass
                                else:
                                    
                                    pass

                            else:
                                if identities > 0.4:
                                    hsp_coord.append(hsp.query_start)
                                    hsp_coord.append(hsp.query_end)
                                else:
                                    pass
                        else:
                            pass
                        
                    #print "%s hsp_coord: %r"%(alignment.title.split(" ")[1],hsp_coord)    
                    if (hsp_coord != []) and (hit_coord != []):
                        overlaps = []
                        for coord in hit_coord:
                            # check for overlapping hits
                            boundary1 = int(coord[0]-100)
                            boundary2 = int(coord[1] +100)
                            overlap_size = 0
                            hit_start = int(min(hsp_coord))
                            hit_end = int(max(hsp_coord))
                            
                            #print "coordinates %d %d %d %d"% (hit_start, hit_end, boundary1, boundary2)
                            
                            #Calculating size of the overlap between hits
                            if boundary1  < hit_end < boundary2:
                                overlap_size = hit_end - max(boundary1,hit_start)
                            
                            elif boundary1 < hit_start < boundary2:
                                overlap_size = min(boundary2,hit_end) - hit_start
                                
                            else:
                                overlap_size = 0
                            
                                
                            if overlap_size < 0.1*(hit_end-hit_start): #True if overlap is  smaller than 10% of hit sequence
                                overlaps.append(False)
                            
                            else:
                                overlaps.append(True)
                                
                            #if (coord[0]-10) < min(hsp_coord) < (coord[1] +10):
                                #x = False
                                ##print "start conflict %r"%str(coord)
                            #elif (coord[0] -10) < max(hsp_coord) < (coord[1] +10):
                                #x = False
                                ##print "end conflict %r"%str(coord)
                            #else:
                                #pass
                        if True not in overlaps: 
                            hit_coord.append((min(hsp_coord),max(hsp_coord),alignment.title.split(" ")[1],seq.id))
                        else:
                            pass

                    elif (hsp_coord != []) and (hit_coord == []):
                        hit_coord.append((min(hsp_coord),max(hsp_coord),alignment.title.split(" ")[1],seq.id))
                        
                    else:
                        pass
                
                if hit_coord != []:

                    for item in hit_coord:
                        #print item[2],item[3]
                        start = end = 0
                        start = item[0]
                        end = item[1]
                        
                        #if item[0] - 200 <0:
                            #start = 0
                        #else:
                            #start = item[0] -200
                            
                        #if item[1] + 200 > len(seq):
                            #end = len(seq)
                        #else:
                            #end = item[1] + 200
                            
                        pair_IDs.append((item[2],item[3],start,end))

                else:
                    pass
            else:
                pass
            

    return pair_IDs

def parse_fasta(fasta_out):
    nlist = [] #nlist contains [query_id,query_len,subject_id,alignment length,alignment start (subject),alignment end (subject),query_align(sequence),subject_align(sequence) ,strand]
    query_id = query_len = subject_id = query_align = subject_align = query_done = strand = ""
    align_start = align_stop = 0
    for line in fasta_out:
        if line[:3] == ">>>" and "library" in line :
            query_id = line.lstrip(">>>").split(",")[0]
            nlist.append(query_id)
            query_len = float(line.split(" ")[1])
            nlist.append(query_len)
            
        elif ">>" in line[:3] and line[:3] != ">>>":
            subject_id = line.lstrip(">>").split(" ")[0]
            nlist.append(subject_id)
        
        elif line[0] != ">" and line[0] != ";" and subject_id != "" and query_id != "" and query_done == "":
            query_align = "%s%s" %(query_align,line.rstrip("\n"))
            
                    
        elif query_align != "" and line[0] == ">" and ">>>" not in line:
            query_done = "query found"
                
        elif subject_id != "" and query_id != "" and query_done == "query found" and line[0] != ">" and line[0] != ";":
            subject_align = "%s%s" %(subject_align,line.rstrip("\n"))
            
            
        elif query_done == "query found" and "al_start" in line:
            align_start = float(line.rstrip("\n").split(" ")[-1])
        
        elif query_done == "query found" and "al_stop" in line:
            align_stop = float(line.rstrip("\n").split(" ")[-1])
            nlist.append(abs(align_stop - align_start)/3)
            nlist.append(align_start)
            nlist.append(align_stop)
            
        elif "tfy_frame" in line:
            strand = line.split(" ")[2].rsplit("\n")[0]
        else:
            pass
            
        if ">--" in line: #not considering secondary hits
            break
        else:
            pass
        
        if ">>><<<" in line: #reached the end of the alignment
            break
        else:
            pass
    
    nlist.append(query_align)
    nlist.append(subject_align)
    nlist.append(strand)

    return nlist

def analyze_fasta(nlist): # returns a tuple with (geneID, pseudogene_call)
    
    try:
        
        if nlist[7] != "":
            nlist[7] = nlist[7].replace("\\","|") # replace escape character in subject alignment \ by |
            #replace all characters representing deleterious mutations by | 
            nlist[7] = nlist[7].replace("/","|")
            #nlist[7] = nlist[7].replace("*","|")
            fixed_seq = nlist[7].translate(None,"*/\|-")
            
            seq = list(nlist[7])
            mutpos = [i for i,j in enumerate(nlist[7]) if j in ["|","/","\\"]]
            #mutpos = [i for i,j in enumerate(nlist[7]) if j in ["|","*","/","\\"]]

            
            #try:
                #i = seq.index("|")

                #while seq.index("|"):
                    #try:
                        #i = seq.index("|")
                        #seq[i] = " "
                        #mutpos.append((i,"frameshift"))            
                    #except ValueError:
                        #pass
            #except ValueError:
                #pass
            
            #try:
                #i = seq.index("*")
                #while seq.index("*"):
                    #try:
                        #i = seq.index("*")
                        #seq[i] = " "
                        #mutpos.append((i,"stop"))            
                    #except ValueError:
                        #pass
            #except ValueError:
                #pass
            
            #stop = nlist[5].find("*")
            #frame1 = nlist[5].find("/")
            #frame2 = nlist[5].find("|")
            
            len_align = float(nlist[3])
            mut = False
            delmut = []
            
            if mutpos != []:
                for i in mutpos:
                    if (0.1*len_align < i < 0.8*len_align):
                        mut = True
                        delmut.append(i)
                        break
                    else:
                        pass
            
            else:
                pass
                
            #mut = (stop != -1 and 0.1*len_align < stop < 0.8*len_align) or \
            #(frame1 != -1 and 0.1*len_align < frame1 < 0.8*len_align)\
            #or (frame2 != -1 and 0.1*len_align< frame2 < 0.8*len_align)
            
            if len_align >= 0.8*float(nlist[1]) and mut: #checks if alignment length is at least 80% of query length and if a mutation is present
                #print "pseudogene because mut = %r"%mut
                pseudo = "pseudo"
            
            elif len_align >= 0.8*float(nlist[1]) and not mut:
                #print "functional"
                pseudo = "functional"
            
            elif len_align < 0.8*float(nlist[1]):
                #print "pseudogene because len_align < %r"% (0.5*float(nlist[1]))
                pseudo = "truncated"
                
            elif len_align < 0.3*float(nlist[1]):
                pseudo = "poor alignment"
                
            else:
                pseudo = "ambiguous"
            
            return (nlist[2],nlist[0],pseudo,nlist[4],nlist[5],nlist[8],fixed_seq,len(mutpos))
        
        else:
            return [] 
    except IndexError:
        return []
        
    
        
#nlist contains [query_id,query_len,subject_id,alignment length,alignment start (subject),alignment end (subject),query_align(sequence),subject_align(sequence), strand]
    #print (nlist[2],nlist[0],pseudo,nlist[4],nlist[5],nlist[8],fixed_seq,len(mutpos))
     


def ortho_intergenic(infile,protdbfile):
    
    
    fileID = infile.split(".")[0]
    oup = "%s_annot_pseudo.gbk"%fileID
    genome = SeqIO.read(infile,"genbank")

    records = intergenic(genome)
    SeqIO.write(records,"intergenic_seq.fa","fasta")
    IDlist = run_blast(records,protdbfile) 
    intergenic_dict = {}
    protdb_dict = SeqIO.index(protdbfile,"fasta")

    for record in records:
        intergenic_dict[record.id] = record

    calls = []
    print "running tfasty on %d pairs of sequences"%len(IDlist)
    for ID in IDlist:
        subject_seq = intergenic_dict[ID[1]]
        slice_start = min(ID[2]-200,0)
        slice_end = max(ID[3] +200 , len(subject_seq))
        subject_seq = subject_seq[slice_start:slice_end]
        query_seq = protdb_dict[ID[0]]
        SeqIO.write(subject_seq, "temp_subj.fasta","fasta")
        SeqIO.write(query_seq, "temp_query.fasta","fasta")
        process = subprocess.Popen("tfasty36 -s BL50 -d 1 -m 10 temp_query.fasta temp_subj.fasta", \
        shell = True,stdout=subprocess.PIPE)
        fasta_out = process.stdout.readlines()
        fasta_list = parse_fasta(fasta_out)
        call = analyze_fasta(fasta_list)
        calls.append(call)

    for x in calls:

        if x != []:
            start_intergenic = int(x[0].split("..")[0])
            end_intergenic = int(x[0].split("..")[1])
            start_hit = start_intergenic + int(min(x[3],x[4]))
            end_hit = start_intergenic + int(max(x[3],x[4]))
            locat = FeatureLocation(start_hit, end_hit)
            
            if x[5] == "f":
                strnd = 1
            elif x[5] == "r":
                strnd = -1
            else:
                pass
            
            NewFeature = SeqFeature(locat, type="CDS", strand = strnd)
            NewFeature.qualifiers["pseudo"] = ['']
            NewFeature.qualifiers["translation"] = x[6]
            NewFeature.qualifiers["note"] = "%d deleterious mutations"%x[7]
            genome.features.append(NewFeature)
        else:
            pass
        
    return(genome)
    
def translate_CDS(genome_file):
    protrecordlist = []
#nuclrecordlist = []

    genome = SeqIO.parse(open(genome_file,"rU"),"genbank")


    for contig in genome:
        
        for feat in contig.features:
            mut = pseudo = ""
            feat_start = feat_end = 0
            if feat.type == "CDS":
                
                try:
                    feat.qualifiers["pseudo"]
                    mut = "_" + feat.qualifiers["note"][0].split(" ")[0]
                    pseudo = "_p" + mut

                except KeyError:
                    pseudo = ""
                

                try:
                    prot = Seq(feat.qualifiers["translation"][0],"generic_protein")
                    
                except KeyError:                    
                    try:
                        #prot = feat.extract(contig.seq).translate(table=11, cds=True)
                        prot = feat.extract(contig.seq).translate(table=11)
                    
                    except TranslationError:
                        prot = ""
                        pass
                    
                    
                
                feat_strand = feat.location.strand
                
                if feat.strand == 1:
                    feat_start = feat.location.start +1 
                    feat_end = feat.location.end
                elif feat.strand == -1:
                    feat_start = feat.location.end 
                    feat_end = feat.location.start +1
                else:
                    pass
                
                feat_id = str(feat_start) + ".." + str(feat_end) + pseudo
                
                
                if prot != "":
                    #nucl = feat.extract(contig.seq)
                    protrecord = SeqRecord(prot,id = feat_id)
                    #nuclrecord = SeqRecord(nucl,id = feat.qualifiers["locus_tag"][0])
                    protrecordlist.append(protrecord)
                    #nuclrecordlist.append(nuclrecord)
                else:
                    pass

                
            else:
                pass
    return(protrecordlist)


def parse_line(line):
    
#example of line:
#ORTHOMCL678(3 genes,3 taxa):     2123579..2124223(Bschum_prot_pseudo) 2500635..2501279(Bpun_prot_pseudo) 439384..440028(Bkir_prot_pseudo)
#ORTHOMCL1552(3 genes,3 taxa):     111201..106373_p_5(Bpun_prot_pseudo) 1184283..1178249_p_3(Bschum_prot_pseudo) 163332..168625_p_13(Bkir_prot_pseudo)
    line = line.rstrip("\n")
    #line = line.split("\t")
    ortho_group = line.split("\t")[0].split("(")[0]
    line = line.split("\t")[1].split(" ")
    #genes is a list of tuples for each orthologous group
    genes = []
    genes.append(ortho_group)
    
    while len(line) > 1:

        gene = line.pop()
        species = gene.split("(")[1].split(")")[0]
        name = gene.split("(")[0]
        location = name.split("_")[0]
        
        start = location.split("..")[0]
        end = location.split("..")[1]
        
        start = start.replace("<","")
        start = int(start.replace(">",""))
        
        end = end.replace("<","")
        end = int(end.replace(">",""))
        
        
        try:
            pseudo = name.split("_")[1] + name.split("_")[2]
        except IndexError:
            pseudo = ""
        
        if start < end:
            strand = 1
        elif end < start:
            strand = -1
        else:
            pass
        
        geneinfo = (species,name)
        genes.append(geneinfo)
    
    return genes

def check_synteny(genome_IDs,gene1,gene2,direction):
    # this function checks if gene1 and gene2 are contiguous in genome
    # gene1 and gene2 names follow the format ddd..ddd_p_d
    # genome_IDs is an ordered list of IDs following the format ddd..ddd_p_d
    # direction is an integer taking the values 1 or -1, telling where to look for the contiguous gene
    # +1 will look in increasing value of position, -1 decreasing
    
    pos_gene1 = 0
    pos_gene2 = 0
    for pos,gene in enumerate(genome_IDs):
        
        if gene == gene1:
            pos_gene1 = pos
            
        elif gene == gene2:
            pos_gene2 = pos
            
        else:
            pass
    
    #choosing a cut-off value of 3, allowing for detection of larger deletions
    #including a whole gene
    
    if (direction == 1) and (1 <= pos_gene2 - pos_gene1 <= 2): 
        return True
    
    elif (direction == -1) and (-2 <= pos_gene2 - pos_gene1 <= -1): 
        return True
    
    else:
        return False
    
def is_pseudo1(gene_name):
    try:
        #print gene_name
        mut = int(gene_name.split("_p_")[1])
        if mut >= 0:
            pseudo = True
        else:
            pseudo = False
    except IndexError:
        pseudo = False
    return pseudo


                

def synteny_blocks(orthomcl_out,genome_folder):
    
    folder = os.listdir(genome_folder)
    inl = open(orthomcl_out,"rU")
    oup = []
    #oup = open("syntenic_blocks.txt","w")

    line = inl.readline()

    groups = []
    #group_dict = {}

    while line != "":
        group = parse_line(line)
        #group_dict = [orthomcl_group,(species1,gene1),(species2,gene2),...]
        groups.append(group)
        #groups = [[orthomcl_group1,(species1,gene1),(species2,gene2)...],..] all genes in a group are orthologs
        #group_dict[group[0]] = group[1:]
        line = inl.readline()
        
    ortho_dict = {}
    #ortho_dict is a nested dictionary containing orthologs for each gene for each genome
    species_list = []
    gene_dict = {}

    for fasta_file in folder:
        in_species = fasta_file.split(".fa")[0]
        species_list.append(in_species)
        sub_ortho_dict = {}
        
        genome_in = SeqIO.parse("%s/%s"%(genome_folder,fasta_file),"fasta")
        genome_IDs = [gene.id for gene in genome_in]
        genome_IDs.sort(key= lambda gene : int(gene.split("..")[0].replace(">","").replace("<","")))
        
        gene_dict[in_species] = genome_IDs
        
        for group in groups:
            print group
            #finding item of group from the reference species
            #each group is a list of tuples
            ref_gene = ""
            sub_genes = []
            for gene in group[1:]:
                #print "%s %s"%(gene[0],gene[1])

                if gene[0] == in_species:
                    #ref_gene = ref_species + "_" + gene[1] + "_" + gene[2] + "_" + gene[3] + "_" + gene[4]
                    ref_gene = gene[1]
                else:
                    #print gene
                    sub_genes.append(gene)
            #print sub_genes
            sub_ortho_dict[ref_gene] = sub_genes
        ortho_dict[in_species] = sub_ortho_dict
        

    align_blocks = []
    # align_blocks is a list of tuple with alignment coordinates like this
    # ((species1, coordinates),(species2,coordinates_segment,coordinates_funct_gene),(species3,coordinates_segment,coordinates_funct_gene))

    for fasta_file in folder:
        print fasta_file
        species = fasta_file.split(".fa")[0]
        ref_genome_IDs = gene_dict[species]
        
        #genome_fasta = SeqIO.parse("%s/%s"%(genome_folder,fasta_file),"fasta")
        #genome_IDs = [gene.id for gene in ref_genome_fasta]
        ##sorting sequences by start position in the genome
        #genome_IDs.sort(key= lambda gene : int(gene.split("..")[0]))
        
        for pos,gene in enumerate(ref_genome_IDs):
            print species,gene
            spe_align_blocks = []
            pseudo_r = True
            pseudo_l = True
            gene_r = ""
            gene_l = ""
            rel_gene_l = ""
            rel_gene_r = ""
            i = j = 1
        
            if is_pseudo1(gene) == True:
                #Looking for functional gene right-flanking gene
                while pseudo_r == True:
                    try:
                        gene_r = ref_genome_IDs[pos+i]
                        pseudo_r = is_pseudo1(gene_r)
                        i = i+1
                    except KeyError:
                        gene_r == ""
                        pseudo_r = False
                    except IndexError:
                        pass
                
                #doing the same on the left
                
                while pseudo_l == True: 
                    try:
                        gene_l = ref_genome_IDs[pos-j]
                        pseudo_l = is_pseudo1(gene_l)
                        j = j+1
                    except KeyError:
                        gene_l == ""
                        pseudo_l = False
                    except IndexError:
                        pass
                
            
                #find orthologs in other genomes
                if (gene_l != "") and (gene_r != ""):
                    #print species,gene_l,gene,gene_r
                    other_species = [i for i in species_list if i != species]
                    
                    for rel_spe in other_species:

                        try:
                            
                            rel_l = ortho_dict[species][gene_l]
                            rel_r = ortho_dict[species][gene_r]
                            
                                                    
                            for x in rel_l:
                                if x[0] == rel_spe:
                                    rel_gene_l = x[1]

                                else:
                                    pass
                                
                            for j in rel_r:
                                if j[0] == rel_spe:
                                    rel_gene_r = j[1]
                                else:
                                    pass
                                
                            
                        except KeyError:
                            rel_gene_l = rel_gene_r = ""
                            pass
                        
                        
                        if rel_gene_l != "" and rel_gene_r != "":
                            
                            print rel_gene_l, rel_gene_r
                            
                            if (is_pseudo1(rel_gene_l) == False) and (is_pseudo1(rel_gene_r) == False):
                                                        
                                if int(rel_gene_l.split("..")[0].replace(">","").replace("<","")) < int(rel_gene_r.split("..")[0].replace(">","").replace("<","")):
                                    rel_strand = 1
                                    rel_seg_coord = (int(min(rel_gene_l.split("_p_")[0].split(".."))),int(max(rel_gene_r.split("_p_")[0].split(".."))))
                                    
                                elif int(rel_gene_l.split("..")[0].replace(">","").replace("<","")) > int(rel_gene_r.split("..")[0].replace(">","").replace("<","")):
                                    rel_strand = -1
                                    rel_seg_coord = (int(min(rel_gene_r.split("_p_")[0].split("..")).replace(">","").replace("<","")),int(max(rel_gene_l.split("_p_")[0].split("..")).replace(">","").replace("<","")))
                                
                                else:
                                    pass
                                                    
                                if check_synteny(gene_dict[rel_spe],rel_gene_l,rel_gene_r,rel_strand) == True:
                                    #finding functional gene between the flanking genes
                                    gene_align_coord = []
                                    for rel_gene in gene_dict[rel_spe]:
                                        
                                        rel_gene_start = min(int(rel_gene.split("..")[0].replace(">","").replace("<","")),int(rel_gene.split("_p_")[0].split("..")[1].replace(">","").replace("<","")))
                                        rel_gene_end = max(int(rel_gene.split("..")[0].replace(">","").replace("<","")),int(rel_gene.split("_p_")[0].split("..")[1].replace(">","").replace("<","")))
                                        
                                        if is_pseudo1(rel_gene) == False:
                    
                                            if (min(rel_seg_coord) < rel_gene_start) and (rel_gene_end < max(rel_seg_coord)) and (rel_gene_end - rel_gene_start > 150):
                                                gene_align_coord.append((rel_gene.split("..")[0],rel_gene.split("_p_")[0].split("..")[1]))
                                                
                                            else:
                                                pass
                                        else:
                                            pass
                                        
                                    if len(gene_align_coord) == 1:
                                        species_block = (rel_spe,rel_seg_coord,rel_strand,gene_align_coord[0])
                                        spe_align_blocks.append(species_block)
                                    else:
                                        pass
                                else:
                                    pass
                                
                            else:
                                pass
                            
                        else:
                            pass
                        
                    if len(spe_align_blocks) == 2:
                        ref_coord = (species,min(gene_l.split("_p_")[0].split("..")),max(gene_r.split("_p_")[0].split("..")))
                        spe_align_blocks.insert(0,ref_coord)
                                            
                        align_blocks.append(spe_align_blocks)
                        #print align_blocks
                            
                    else:
                        pass
                    
                else:
                    pass

    oup.append("query species\tseg start\tseg end\tsubject species\tseg start\tseg end\tstrand\tgene start\tgene end\tsubject species\tseg start\tseg end\tstrand\tgene start\tgene end\n")

    for block in align_blocks:
        print block
        oup.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%r\t%r\t%r\t%s\t%s\n"\
    %(block[0][0],block[0][1],block[0][2],block[1][0],block[1][1][0],block[1][1][1],block[1][2],block[1][3][0],block[1][3][1],block[2][0],block[2][1][0],block[2][1][1],block[2][2],block[2][3][0],block[2][3][1]))
        #oup.write("%r\n"%block)
        #oup.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%r\t%r\t%r\t%s\t%s\n"\
    #%(block[0][0],block[0][1],block[0][2],block[1][0],block[1][1][0],block[1][1][1],block[1][2],block[1][3][0],block[1][3][1],block[2][0],block[2][1][0],block[2][1][1],block[2][2],block[2][3][0],block[2][3][1]))

    return oup
    
    
def parse_synteny(line):
    info = []
    line = line.rstrip("\n")
    line= line.split("\t")
    
    if line[0] != "query species":
            
        pseudo_spe = line[0]
        pseudo_seg = (int(line[1]),int(line[2]))
        sub_spe1 = line[3]
        sub_seg1 = (int(line[4]),int(line[5]))
        sub_seg1_strand = int(line[6])
        sub_gene1 = (int(line[7]),int(line[8]))
        sub_spe2 = line[9]
        sub_seg2 = (int(line[10]),int(line[11]))
        sub_seg2_strand = int(line[12])
        sub_gene2 = (int(line[13]),int(line[14]))
        
        info = [pseudo_spe,pseudo_seg,sub_spe1,sub_seg1,sub_seg1_strand,sub_gene1,\
            sub_spe2,sub_seg2,sub_seg2_strand,sub_gene2]
        
        
    else:
        info = []
    
    return info

def is_pseudo(gene_name):
    try:
        pseudo = int(gene_name.split("_p_")[1])
        #pseudo = "_p_" in gene_name
    except IndexError:
        pseudo = 0
    return pseudo

def extract_coordinates(gene):
    strand = 0
    gene= gene.split("_p_")[0]
    start = int(gene.split("..")[0])
    end = int(gene.split("..")[1])
    
    
    if start < end:
        strand = 1
    elif start > end:
        strand = -1
    else:
        pass
    gene_coord = (min(start,end),max(start,end),strand)
    
    return gene_coord
    
def run_muscle(ref_seq,seqlist):
    #seqlist is a list of sequences. 
    #refseq is the sequence used as reference
    ref_seq.id
    seqlist.append(ref_seq)
    
    SeqIO.write(seqlist,"temp.fa","fasta")
    #os.system("muscle -in temp.fa -out temp.out.clw -clwstrict -quiet")
    os.system("einsi --clustalout temp.fa > temp.out.clw")
    align = AlignIO.read("temp.out.clw","clustal")
    
    return align

    
    
            

    

folder = os.listdir(genome_folder)
print folder
os.system("mkdir temp")
print os.listdir(".")

if "prot_db.fa" not in os.listdir("."):
    os.system("mkdir genomes")
    os.system("mkdir fasta_files")

    for inl in folder:
        print "Reading genome file %s"%inl
        shortname = inl.split(".gbk")[0][-15:]
        print shortname
        os.system("union -feature -osf genbank -outseq temp/%s %s/%s"%(shortname,genome_folder,inl))
        
    foldercat = os.listdir("temp")

    for genome in foldercat:
        #removing pseudogenes from the annotation
        #genome files have to be concatenated

        record = SeqIO.read("temp/%s"%genome,"genbank")
        seq = record.seq
        featlist = []
        seqrec = SeqRecord(seq,record.id,description=record.description)
        
        for feat in record.features:
            try:
                pseudo = feat.qualifiers["pseudo"]
            except KeyError:
                featlist.append(feat)
                
        seqrec.features = featlist
        SeqIO.write(seqrec,"genomes/%s"%genome,"genbank")
        
    #extract protein sequences for orthomcl run

    os.system("rm -fR temp")

    genome_folder = os.listdir("genomes")

    os.system("mkdir temp")

    for inl in genome_folder:
        name = inl.split(".gbk")[0]+"_prot.fa"
        prots = extractprots("genomes/%s"%inl)
        SeqIO.write(prots,"temp/%s"%name,"fasta")
        
    os.system("cat temp/*.fa > prot_db.fa")
    os.system("rm -fR temp")
    
else:
    pass


if "pseudo_genomes" not in os.listdir("."):
    
    print "running blastx"
    #running blastx on intergenic regions to find potential pseudogenes
    os.system("mkdir pseudo_genomes")

    for x in genome_folder:
        pseudo_genome = ortho_intergenic("genomes/%s"%x,"prot_db.fa")
        SeqIO.write(pseudo_genome,"pseudo_genomes/%s_prot.gbk"%x.split(".gbk")[0],"genbank")
        
    #preparin fasta files for the orthomcl run
    os.system("rm -fR genomes")

else:
    pass

fasta_folder = os.listdir("fasta_files")
print fasta_folder

if fasta_folder == []:

    pseudo_folder =  os.listdir("pseudo_genomes")
    print "translating CDS"

    for i in pseudo_folder:
        
        pseudo_prot_list = translate_CDS("pseudo_genomes/%s"%i)
        SeqIO.write(pseudo_prot_list,"fasta_files/%s.fa"%(i.split(".gbk")[0]),"fasta")
        
else:
    pass

if "one_to_one.txt" not in os.listdir("."):
    print "preparing files for orthomcl"
    fasta_folder = os.listdir("fasta_files")
    print fasta_folder
    os.system("ln -s ~/ORTHOMCL/orthomcl.pl")
    os.system("ln -s ~/ORTHOMCL/orthomcl_module.pm")
    
    for fasta_in in fasta_folder:
        print fasta_in
        path = os.path.abspath("fasta_files/%s"%fasta_in)
        print path
        #creating symbolic links into the orthomcl sample folder
        os.system("ln -sf %s -t ~/ORTHOMCL/sample_data"%path)

        
    #running orthomcl
    print "running orthomcl"
    process = subprocess.Popen("perl orthomcl.pl --mode 1 --fa_files %s"%(",".join(fasta_folder)),\
        shell= True,stderr=subprocess.PIPE)
                            
    #stdout = process.stderr.readlines()[1].replace(" ","")
    stdout = process.stderr.readlines()
    #print stdout

    path2 = stdout[1].split("\n")[0].replace(" ","") #working directory of orthomcl
    print "processing Orthomcl output"

    os.system("ln -s %sall_orthomcl.out"%path2)
    os.system("grep \"(3 genes,3 taxa)\" all_orthomcl.out > one_to_one.txt")
    
else:
    pass

print "Finding syntenic blocks"
synt_blocks = synteny_blocks("one_to_one.txt","fasta_files")

synt_out = open("syntenic_blocks.txt","w")

for line in synt_blocks:
    synt_out.write(line)

synt_out.close()

#it's necessary to reset the mafft binaries before running the program

try:
    del os.environ["MAFFT_BINARIES"]
except KeyError:
    pass


os.system("mkdir alignment_whole")
os.system("mkdir alignments")

#folder = os.listdir(genome_folder)
    
inl = open("syntenic_blocks.txt","rU")

line = inl.readline()

blocks = []
seq_dict = {}

while line != "":
    info_genes = parse_synteny(line)
    
    if info_genes != []:
        blocks.append(info_genes)
    else:
        pass
    
    line = inl.readline()
    
pseudo_folder =  os.listdir("pseudo_genomes")
    
for inl in pseudo_folder:
    seq_dict[inl.split(".gbk")[0]] = SeqIO.read("pseudo_genomes/%s"%(inl),"genbank").seq

print "extracting sequences from syntenic blocks"
for item in blocks:

    pseudo_spe = item[0]
    pseudo_seg = item[1] # is a tuple (start,end)
    sub_spe1 = item[2]
    sub_seg1 = item[3] # is a tuple (start,end)
    sub_seg1_strand = item[4]
    sub_gene1 = item[5] # is a tuple (start,end)
    sub_spe2 = item[6] 
    sub_seg2 = item[7] # is a tuple (start,end)
    sub_seg2_strand = item[8]
    sub_gene2 = item[9] # is a tuple (start,end)
    species_list = [pseudo_spe,sub_spe1,sub_spe2]
    
    sub_gene_coord1 = ()
    sub_gene_coord2 = ()
    sub_sequences = []
    coord_dict = {}
    #print "new item"
    
    pseudo_seq = seq_dict[pseudo_spe][pseudo_seg[0]:pseudo_seg[1]]
    name = "%s_%d..%d"%(pseudo_spe,pseudo_seg[0],pseudo_seg[1])
    
    pseudo_slice = SeqRecord(pseudo_seq,id = pseudo_spe, name = pseudo_spe)
    
    sub_seq1 = seq_dict[sub_spe1][sub_seg1[0]:sub_seg1[1]]
    
    if sub_seg1_strand == -1:
        sub_seq1 = sub_seq1.reverse_complement()
        sub_gene_coord1 = (sub_seg1[1] - max(sub_gene1), sub_seg1[1] - min(sub_gene1))
        
    elif sub_seg1_strand == 1:
        sub_gene_coord1 = (min(sub_gene1) - sub_seg1[0], max(sub_gene1) - sub_seg1[0])
        
    else:
        print "error strand sub seq1"
        pass
    
    sub_slice1 = SeqRecord(sub_seq1,id= sub_spe1,name=sub_spe1)
    sub_sequences.append(sub_slice1)
    coord_dict[sub_spe1] = sub_gene_coord1
    
    sub_seq2 = seq_dict[sub_spe2][sub_seg2[0]:sub_seg2[1]]
    
    if sub_seg2_strand == -1:
        sub_seq2 = sub_seq2.reverse_complement()
        sub_gene_coord2 = (sub_seg2[1] - max(sub_gene2), sub_seg2[1] - min(sub_gene2))
    
    elif sub_seg2_strand == 1:
        sub_gene_coord2 = (min(sub_gene2) - sub_seg2[0], max(sub_gene2) - sub_seg2[0])
        
    else:
        print "error strand sub seq2"
        pass
    sub_slice2 = SeqRecord(sub_seq2,id= sub_spe2,name=sub_spe2)
    sub_sequences.append(sub_slice2)
    coord_dict[sub_spe2] = sub_gene_coord2
    
    print "gene coord",sub_gene_coord1,sub_gene_coord2,sub_gene1,sub_gene2
    

    print "running alignment"
    align = run_muscle(pseudo_slice,sub_sequences)
    AlignIO.write(align,"alignment_whole/%s.aln"%name,"clustal")
    
    start_coord = []
    end_coord = []
    align_slice = []
    
    poor_alignment = False
    
    for x in align:
        # print len(x.seq)
        if x.seq.count("-") > 0.3*len(x.seq):
            poor_alignment = True
        else:
            pass
            
    if poor_alignment == False:
        for record in align:
                
            seq_rec = record.seq
            seq_id = record.id
            
            # Below is a necessary step since sequence names are often truncated in clustal
            # format
            for spe in species_list:
                if seq_id[:10] == spe[:10]:
                    seq_id = spe
                else:
                    pass
            
            if seq_id != pseudo_spe:
            
                slice_start = slice_end = 0
            
                # for coord in sub_slice_coords:
                    #print coord[0],seq_id
                coord = coord_dict[seq_id]
            
                slice_start = coord[0]
                slice_end = coord[1]
            
                #Calculating where the pseudogene alignment starts. This works by counting the characters in the alignment 
                #without the gaps.
                for i,j in enumerate(seq_rec):
                    #print i, seq_rec[:i].count("-"), i - seq_rec[:i].count("-")
                
                
                    if i - seq_rec[:i].count("-") == slice_start:
                        start_coord.append(i)
                    
                    elif i - seq_rec[:i].count("-") == slice_end:
                        #print i
                        end_coord.append(i)
                    else:
                        pass
                else:
                    pass
        
        if start_coord != end_coord:
            align_slice = [max(start_coord),min(end_coord)]
                
            align = align[:,align_slice[0]:align_slice[1]]

            for record in align:
                seq_id = record.id        
                seq = str(record.seq)
                
            AlignIO.write(align,"alignments/%s.aln"%name,"clustal")
        else:
            pass
    else:
        pass




    


    
