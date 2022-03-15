import pyemma
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# set the font name for a font family
plt.rcParams.update({'font.sans-serif':'Arial'})

import pandas as pd
"""
tic_list = np.loadtxt('columns_tic.txt')[:,0]
traj = md.load('tic1.pdb')
dssp = md.compute_dssp(traj)
df = pd.DataFrame(dssp).T 
df.columns=tic_list
"""
mutations = [152,153,166,182,183,205]

df = pd.read_pickle('tic1_dssp.pkl')
states = np.load('states.npy', allow_pickle=True)

state_1 = states[0]
state_2 = states[1]
state_3 = states[2]

refdf = pd.read_pickle('ref_dssp.pkl')
helix = np.array([i for i in refdf.loc[refdf[0] == 'H'].index])
sheet = np.array([i for i in refdf.loc[refdf[0] == 'E'].index])
coils = np.array([i for i in refdf.loc[refdf[0] == 'C'].index])

#df = pd.read_pickle('/mnt/d/plots/Sec_struc_state_plot/Tic1_dssp.pkl')
#tic_1 = np.loadtxt('/mnt/d/owlsnest/500_cluster_4/tic1_cc.txt',dtype=int)
#df.columns=tic_1

ss_1 = df[state_1]

ss_2 = df[state_2]

ss_3 = df[state_3]

sts_list = [ss_3,ss_1,ss_2]


def ss_counter(df,helix,sheet,coil):
    cnt_array = np.zeros(100)
    for i,struc in enumerate(df.T):
        if i in helix:
            row = df.iloc[i]
            cnt_array[i] = row[row =='H'].count()
        elif i in sheet:
            row = df.iloc[i]
            cnt_array[i] = row[row =='E'].count()
        elif i in coil:
            row = df.iloc[i]
            cnt_array[i] = row[row =='C'].count()
    return cnt_array

### Define Gaussian Smoothing
def gaussian_smooth(x,y,sigma):
    pi=np.pi
    smoothed_vals = np.zeros(y.shape)
    for i,x_0 in enumerate(x):
        kernel = np.sqrt(2*pi*sigma**2)**(-1)*np.exp(-(x-x_0)**(2)/2*sigma**2) 
        kernel = kernel/np.sum(kernel)
        smoothed_vals[i] = np.sum(y*kernel)
    return smoothed_vals

def plt_distro(list):
    fig,ax = plt.subplots(figsize=(6,2.5), dpi=300, sharey=True)
    states = [state_3,state_1,state_2]
    colors = ['red','gray','blue']
    labels = ['U','I','N']
    helix1 = np.arange(15,26)
    helix2 = np.arange(33,44)
    turn = np.arange(45,48)
    helix3 = np.arange(53,70)
    beta4 = np.arange(73,77)
    beta5 = np.arange(86,90)
    structure_list = [helix1,helix2,turn,helix3,beta4,beta5]
    structure_colors = ['r', 'pink', 'grey', 'cyan', 'y', 'y']
    structure_str = [r'$\alpha1$',r'$\alpha2$',r'$\alpha T$',r'$\alpha3$',r'$\beta4$',r'$\beta5$']
    ax.set_ylabel('Probability')
    ax.set_xlabel('Residue')
    for i, struc in enumerate(list):
        structure = ss_counter(struc,helix,sheet,coils)
        x = np.arange(len(structure)) + 150
        y = gaussian_smooth(x,structure/len(states[i]),0.5) 
        ax.plot(x,y, color=colors[i], label=labels[i], linewidth=2)
    ax.legend()
    for i , struc in enumerate(structure_list):
        ax.axvspan(np.min(struc)+150, np.max(struc)+150, alpha=0.2, color=structure_colors[i])
        ax.text( np.median([np.min(struc)+150,np.max(struc)+150])-1,1.08,structure_str[i])
    for mut in mutations:
        ax.axvline(mut,c='black',linestyle=':',alpha=0.5)
    fig.tight_layout()
    fig.savefig('sec_struc_test.pdf')

plt_distro(sts_list)
