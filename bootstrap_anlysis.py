#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reads count tables of insertion or deletion
Performs bootstrapping analysis and creates plots
Creates datafiles used in regression analysis in R
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import sys
 
def bootstrap_resample(pop):
    '''
    Randomly resample a given population (with replacement)
    '''
    new_pop = []
    for i in range(len(pop)):
        new_pop.append(np.random.choice(pop))
    return(new_pop)
    
def pop_to_counts(pop, c):
    '''
    Counts occurence of of indels of certain length
    '''
    return([pop.count(i+1) for i in range(c)])

def pwl(x, c, m, c0):
    '''
    Power-law function
    '''
    return c0 + np.power(x,m) * c

def subdata(data):
    '''
    Removes 3n indel data for fitting power-law curve
    '''
    new = []
    for i, x in enumerate(data):
        if (i+1)%3 != 0:
            new.append(x)
    return(new)

script, table, indeltype = sys.argv
if indeltype == 'Deletion':
    short = 'del'
elif indeltype == 'Insertion':
    short = 'ins'
else:
    sys.exit('Invalid indel type')

data = pd.read_table(table, header=0, index_col=0)
print(data)

#Bootstrap parameters
p = .05    #p-value for confidence intervals
n = 10000    #nr. of bootstrap replicates
l = int((p/2)*n)    #confidence interval boundaries
c = 15  #highest indel length to consider in the analysis

species = list(data.index)
data = data.iloc[:,:c]
for s in range(len(species)):
    print('\n'+species[s]+'\n')
    plt.figure(figsize=(10,8))
    test = data.iloc[s,]
    population = []
    for i in range(c):
        for j in range(test[i]):
            population.append(i+1)
    #fit power-law
    xdata = subdata(list(range(1,c+1)))
    ydata = subdata(list(test))   
    curve = optimize.curve_fit(pwl, xdata, ydata, maxfev=5000, method='trf')
    #create figure
    plt.plot(list(np.arange(1,c+1,0.1)), pwl(list(np.arange(1,c+1,0.1)), *curve[0]))
    if 'B.' in species[s]:
        plt.title('Ca. '+species[s], style='italic')
    else:
        plt.title(species[s])
    plt.xlabel('{} size (nt)'.format(indeltype))
    plt.ylabel('Count')
    #perform bootstrap analysis
    diff = {}
    count = {}
    for i in range(n):
        new_pop = bootstrap_resample(population)
        for j in range(1,c+1):
            diff[j] = diff.get(j, [])
            diff[j].append(new_pop.count(j) - max(pwl(j, *curve[0]), 0))
            count[j] = count.get(j, [])
            count[j].append(new_pop.count(j))
    plt.errorbar(diff.keys(), [np.mean(count[x]) for x in sorted(diff.keys())], yerr=[[np.mean(count[x])-sorted(count[x])[l] for x in sorted(diff.keys())], [sorted(count[x])[-l]-np.mean(count[x]) for x in sorted(diff.keys())]], fmt='o', capsize=3)
    plt.savefig('{}_{}.png'.format(species[s], short), dpi=600, format='png')
    OUT = open('bootstrap.out', 'w+')
	for key in sorted(diff.keys()):
        print('Count: {}\tAverage Difference: {:.3f}\t95% CI: [{:.2f}, {:.2f}]]'.format(key, np.mean(diff[key]), sorted(diff[key])[l], sorted(diff[key])[-l]), file=OUT)
	OUT.close()
    