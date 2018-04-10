#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Search for deletion motifs in a folder containing the alignment files in CLUSTAL format
"""

def calc_distance(seq1, seq2):
    """
    Calculates hamming distance between two sequences
    """
    distance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            distance += 1
    return(distance)

def collapse_gaps(gapslist):
    """
    Converts  a list of sequence positions of gap ('-') characters in a string,
    to list containing the start and end positions of gaps
    """
    if gapslist:
        gaps = []
        gap = (0,0)
        for element in sorted(gapslist):
            if gap == (0,0):
                gap = [element, element]
            else:
                if element-gap[1] == 1:
                    gap[1] = element
                else:
                    gap[1] += 1
                    gaps.append(gap)
                    gap = [element,element]
        gap[1] += 1
        gaps.append(gap)
        return(gaps)
    else:
        return([])

def test_forward(ref, pseudo, gap):
    '''
    Search for 'forward' motifs in a certain gap
    Forward motif: motif where the start of the deleted sequence matches with
    the sequence right after the gap
    '''
    maxlen = gap[1]-gap[0]
    mismatch, i = 0,0
    subref = ''
    subpse = ''
    while mismatch < 2:
        if gap[1]+i < len(pseudo):
            refnuc = ref[gap[0]+i]
            psenuc = pseudo[gap[1]+i]
            if refnuc == psenuc:
                subref += refnuc
                subpse += psenuc
            elif mismatch == 0:
                subref += refnuc
                subpse += psenuc
                mismatch += 1
            elif mismatch == 1:
                mismatch += 1
            elif refnuc == '' or psenuc == '':
                break
            i += 1
        else:
            break
    return(subref[:maxlen], subpse[:maxlen])

def test_reverse(ref, pseudo, gap):
    '''
    Search for 'reverse' motifs in a certain gap
    Reverse motif: motif where the end of the deleted sequence matches with
    the sequence right before the gap
    '''
    maxlen = gap[1]-gap[0]
    mismatch, i = 0,0
    subref = ''
    subpse = ''
    while mismatch < 2:
        if gap[1]-i >= 0:
            refnuc = ref[gap[1]-1-i]
            psenuc = pseudo[gap[0]-1-i]
            if refnuc == psenuc:
                subref = refnuc+subref
                subpse = psenuc+subpse
            elif mismatch == 0:
                subref = refnuc+subref
                subpse = psenuc+subpse
                mismatch += 1
            elif mismatch == 1:
                mismatch += 1
            elif refnuc == '' or psenuc == '':
                break
            i += 1
        else:
            break
    return(subref[len(subref)-maxlen:], subpse[len(subpse)-maxlen:])
           

import os
total_gaps = 0
small_gaps = 0
found_motifs_f = 0
found_motifs_r = 0
fw_motifs = []
rv_motifs = []
nr_motifs = 0
all_motifs = []
motif_deletions = []
all_seqs = []
del_length = 0

for file in sorted(os.listdir()):
    if file.endswith('.aln'):
        print('===== ' + file + ' =====\n')
        seqdict = {}
        genes = []
        for line in open(file):
            line = line.upper()
            if 'CLUSTAL' not in line and line.strip() and '*' not in line:
                gene, sequence = line.split()
                seqdict[gene] = seqdict.get(gene, '')
                seqdict[gene] += sequence
                if gene not in genes:
                    genes.append(gene)
        pseudo = gene
        gaps = []
        for i in range(len(seqdict[pseudo])):
            if seqdict[pseudo][i] == '-':
                test = ''
                for key in seqdict.keys():
                    test += seqdict[key][i]
                if test.count('-') >= 1:
                    gaps.append(i)
        deletions = collapse_gaps(gaps)
        for deletion in deletions:
            if deletion and deletion[0] != 0 and deletion[1] != len(seqdict[genes[-1]]):
                del_status = 0
                for gene in genes:
                    all_seqs.append(seqdict[gene].replace('-', ''))
                    if '-' in seqdict[gene][deletion[0]:deletion[1]]:
                        del_status += 1
                if del_status == 1:
                    total_gaps += 1
                    del_length += deletion[1] - deletion[0]
                    if deletion[1] - deletion[0] <= 15:
                        small_gaps += 1
                    for gene in genes:
                        print(gene[:8], seqdict[gene][max(0, deletion[0]-10):deletion[1]+10])
                    print()
                    print('FORWARD')
                    f_motifs = []
                    for i in [0, 1]:
                        a, b = test_forward(seqdict[genes[i]], seqdict[pseudo], deletion)
                        if len(a) >= 4:
                            print(genes[i][:8], a)
                            print(pseudo[:8], b)
                            print()
                            if not f_motifs:
                                found_motifs_f += 1
                            f_motifs.append(b)
                        elif len(a) < 4 and len(a) == deletion[1]-deletion[0] and calc_distance(a,b) == 0:
                            print(genes[i][:8], a)
                            print(pseudo[:8], b)
                            print()
                            if not f_motifs:
                                found_motifs_f += 1
                            f_motifs.append(b)
                    if f_motifs:
                        fw_motifs.append(max(f_motifs, key=len))
                        all_motifs.append(max(f_motifs, key=len))
                    print('REVERSE')
                    r_motifs = []
                    for i in [0, 1]:
                        a, b = test_reverse(seqdict[genes[i]], seqdict[pseudo], deletion)
                        if len(a) >= 4:
                            print(genes[i][:8], a)
                            print(pseudo[:8], b)
                            print()
                            if not r_motifs:
                                found_motifs_r += 1
                            r_motifs.append(b)
                        elif len(a) < 4 and len(a) == deletion[1]-deletion[0] and calc_distance(a,b) == 0:
                            print(genes[i][:8], a)
                            print(pseudo[:8], b)
                            print()
                            if not r_motifs:
                                found_motifs_r += 1
                            r_motifs.append(b)
                    if r_motifs:
                        rv_motifs.append(max(r_motifs, key=len))
                        all_motifs.append(max(r_motifs, key=len))
                    if r_motifs or f_motifs:
                        nr_motifs += 1
                        motif_deletions.append(deletion[1]-deletion[0])

                    print()
                    print('-'*50+'\n')

print('Deletions: ', total_gaps)
print('Deletions <= 15:', small_gaps)
print('Forward motifs found: ', found_motifs_f)
print('Reverse motifs found: ', found_motifs_r)
print('Deletions with motifs found: ', nr_motifs)