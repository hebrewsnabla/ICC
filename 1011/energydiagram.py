# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 00:39:29 2019

@author: wsr
"""

import matplotlib.pyplot as plt
import sys
import numpy as np

rxn = 'cl.txt' 
#rxn = sys.argv[1]

kcal = True
s2s = 2.0
width = 0.5


with open(rxn,'r') as f:
    rxndata = f.read()
    
energies_tot = []
for item in rxndata.split('--'):
    if '#' not in item: continue
    item = item.split('\n')
    energies = item[2].split()
    energies = np.array([float(i) for i in energies])
    energies -= energies[0]
    if kcal: energies *= 627.509
    energies = [round(i,1) for i in energies]
    energies_tot.append(energies)
energies_tot = np.hstack(energies_tot)
erange = energies_tot.max() - energies_tot.min()

legends = []
#energies_tot = []
for item in rxndata.split('--'):
    if '#' not in item: continue
    item = item.split('\n')
    legend = item[0].split('#')[1]
    legends.append(legend)
    color = item[0].split('#')[2].strip()
    #print(color)
    states = item[1].split()
    energies = item[2].split()
    if len(states)!=len(energies):
        print("numbers of states and energies do not match in %s" % legend)
        exit(-1)
    energies = np.array([float(i) for i in energies])
    energies -= energies[0]
    if kcal: energies *= 627.509
    energies = [round(i,1) for i in energies]
    #energies_tot.append(energies)
            
    fig = plt.figure()
    x = np.array(range(1,len(states)+1)) * s2s
    stages = [[xi-width/2,xi+width/2] for xi in x]
    slopes = [[x[i]+width/2, x[i+1]-width/2] for i in range(len(x)-1)]
    print(stages, slopes)
    for j in range(len(stages)): 
        plt.plot(stages[j], [energies[j], energies[j]], color, linewidth=3.5)
        if 'TS' in states[j]:
            flt = -4
        else:
            flt = -4
        plt.text(x[j] - len(states[j])*0.14*width, energies[j] +(-1.5 + flt)*erange/46, states[j], fontdict={'size':16,'color':color})
        plt.text(x[j] - len(str(energies[j]))*0.12*width, energies[j] +(-1.5 - flt)*erange/46, energies[j], fontdict={'size':16,'color':color})
    for k in range(len(slopes)):
        plt.plot(slopes[k], [energies[k], energies[k+1]], color, dashes=[2,2])
    
    ax=plt.gca() #gca=get current axis
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')


plt.ylim((energies_tot.min() - erange*0.3, energies_tot.max() + erange*0.3))
plt.xticks([])
plt.ylabel(r"$\Delta G_{298K}$(kcal/mol)")
plt.legend(labels=legends,loc='best')