import numpy as np
import pandas as pd
import pyemma
import mdtraj as md
import matplotlib.pyplot as plt
from itertools import combinations

mutations = [152,153,166,182,183,205]
# Plotting Function
def plt_heat(df):
    plt.rcParams.update({'font.size': 20})
    fig, axes = plt.subplots(3,1,figsize=(14,21),sharex=True, dpi=300)
    labels = df.columns[2::]
    for i,ax in enumerate(axes.flat):
        x = df[df.columns[0]]
        y = df[df.columns[1]]
        z = df[df.columns[i+2]]
        im = ax.hist2d(x,y,bins=50,vmin=0, vmax=1, weights=z, cmap='bone_r')
        cbar1 = fig.colorbar(im[3], ax=ax)
        ax.text(220,170,labels[i],fontsize=30)
        ax.set_ylabel('Residue index',fontsize=20)
        ax.set_xticks(np.arange(150,251,25))
        ax.set_yticks(np.arange(150,251,25))
        ax.set_xlim(150,251)
        ax.set_ylim(150,251)
        ax.grid(alpha=0.3)
        ax.set_aspect(1)
        cbar1.set_label('Contact Frequency')
        for mut in mutations:
            ax.axvline(mut,c='black',linestyle=':')
    axes.flat[2].set_xlabel('Residue Index',fontsize=20)
    fig.tight_layout()
    fig.savefig('Nat_con_freq.png')

probabilities = pd.read_pickle('probabilites.pkl')

plt_heat(probabilities)
