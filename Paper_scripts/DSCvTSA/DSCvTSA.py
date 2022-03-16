import os, sys
import numpy as np
from matplotlib import pyplot as plt

labels = ['WT','S152R','S153R','A166G','A166V','T182M','L183P','L183R','S205T']
TSA_tm = np.array([(50.5,0.8),(51.1,0.1),(51.8,0.1),(46.4,0.5),(56.3,0.9),(46.4,0.1),(40.0,0.6),(33.6,1.1),(52.7,1.7)])
DSC_tm = np.array([(57.8,0.3),(57.0,1),(58.6,1),(55.1,1),(61.9,1),(51.5,1),(48.2,1),(41.7,1),(58.3,1)])

DSC = {}
TSA = {} 
for i,res in enumerate(labels):
    DSC[res] = DSC_tm[i]
    TSA[res] = TSA_tm[i]

markers = ['o','.','v','^','1','s','P','*','x']
fig,ax = plt.subplots(dpi=300)
for i,label in enumerate(labels):
   ax.scatter(DSC[label][0],TSA[label][0],marker=markers[i],label=label)
   ax.errorbar(DSC[label][0],TSA[label][0],xerr=DSC[label][1],yerr=TSA[label][1],linewidth=0.8,alpha=0.5)
ax.set_xlabel(r'$Tm_{DSC}(C^\circ)$')
ax.set_ylabel(r'$Tm_{TSA}(C^\circ)$')
fig.legend(loc='lower right',ncol=2, borderaxespad=5)
fig.tight_layout
fig.savefig('DSCvTSA.png')


