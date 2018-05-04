#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script that calculates intervals and their size distribution between the occurence of kmers in a genome

Usage: python3 kmer_separation.py genome kmer-length
Output: Text files containing counts of interval lengths and png-files with the plotted distribution,
        either for each k-mer separately or pooled over all k-mers with the specified length
"""

from sys import argv
import matplotlib.pyplot as plt

script, genome, kmer = argv
kmer = int(kmer)

genomedict = {}
sequence = ''
for line in open(genome):
    if line.startswith('>'):
        contig = line.strip().lstrip('>')
        if sequence:
            genomedict[contig] = sequence
            sequence = ''
    else:
        sequence += line.strip().upper()
if sequence:
    genomedict[contig] = sequence

kmerdict = {}
for contig in genomedict.keys():
    sequence = genomedict[contig]
    for i in range(len(sequence)-kmer):
        subseq = sequence[i:i+kmer]
        kmerdict[subseq] = kmerdict.get(subseq, [])
        kmerdict[subseq].append(i)

all_counts = [0 for i in range(0, 100)]
for subseq in sorted(kmerdict.keys()):
    if 'N' not in subseq:
        print(subseq)
        starts = kmerdict[subseq]
        intervals = [starts[i]-starts[i-1] for i in range(1, len(starts))]
        all_counts = [x+y for x,y in zip(all_counts, [intervals.count(i) for i in range(0, 100)])]
        lengths = range(0, 100)
        t_lengths = [i for i in lengths if i%3 == 0]
        nt_lengths = [i for i in lengths if i%3 != 0]
        bins = [intervals.count(i) for i in range(0, 100)]
        three_bins = [intervals.count(i) for i in range(0,100) if i%3 == 0]
        non_three_bins = [intervals.count(i) for i in range(0, 100) if i%3 != 0]
#        plt.figure(figsize=(10,8))
#        plt.scatter(t_lengths, three_bins, color='red')
#        plt.scatter(nt_lengths, non_three_bins)        
#        plt.title(subseq)
#        plt.savefig(genome.split('.')[0]+'_'+subseq+'.png', format='png')
#        plt.close()
#        OUT = open(genome.split('.')[0]+'_'+subseq+'.txt', 'w+')
#        print('length,count,3n', file=OUT)
#        for i in range(kmer, 100):
#            print(i, intervals.count(i), i%3==0, sep=',', file=OUT)
#        OUT.close()
OUT = open(genome.split('/')[-1].split('.')[0]+'_pooled.txt', 'w+')
print('length,count,3n', file=OUT)
for i, j in zip(range(0, 100), all_counts):
    print(i, j, i%3==0, sep=',', file=OUT)
    
t_counts = [j for i,j in enumerate(all_counts) if i%3 == 0]
nt_counts = [j for i,j in enumerate(all_counts) if i%3 != 0]   
    
plt.figure(figsize=(10,8))
plt.scatter(t_lengths, t_counts, color='red')
plt.scatter(nt_lengths, nt_counts)        
plt.title('Pooled')
plt.savefig(genome.split('/')[-1].split('.')[0]+'_pooled_{}.png'.format(str(kmer)), format='png')
plt.close()