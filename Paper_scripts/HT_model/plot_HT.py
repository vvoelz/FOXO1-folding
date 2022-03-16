import numpy as np
import pandas as pd
import pyemma
import mdtraj as md
import matplotlib.pyplot as plt

# SASA Data import and Sifting

# Import SASA
sasa = np.load('all_sasa.npy', allow_pickle=True)

# Get Metastable States
states = np.load('states.npy',allow_pickle=True)
state_1 = states[0]
state_2 = states[1]
state_3 = states[2]

# Assign CC Samples to States
list_1 = []
list_2 = []
list_3 = []
for i,array in enumerate(sasa):
    if i in state_1:
        list_1.append(array) 
    elif i in state_2:
        list_2.append(array)
    elif i in state_3:
        list_3.append(array)

# Get Average per Residue SASA per State
array_1 = np.zeros((len(list_1),100))
array_2 = np.zeros((len(list_2),100))
array_3 = np.zeros((len(list_3),100))

for i,array in enumerate(list_1):
    df = pd.DataFrame(array)
    mean = df.mean(axis=0)
    for j, res in enumerate(mean):
        array_1[i,j] = res
    df_1 = pd.DataFrame(array_1).mean(axis=0)
    sd_1 = pd.DataFrame(array_1).std(axis=0)
    
for i,array in enumerate(list_2):
    df = pd.DataFrame(array)
    mean = df.mean(axis=0)
    for j, res in enumerate(mean):
        array_2[i,j] = res
    df_2 = pd.DataFrame(array_2).mean(axis=0)
    sd_2 = pd.DataFrame(array_2).std(axis=0)
    
for i,array in enumerate(list_3):
    df = pd.DataFrame(array)
    mean = df.mean(axis=0)
    for j, res in enumerate(mean):
        array_3[i,j] = res
    df_3 = pd.DataFrame(array_3).mean(axis=0)
    sd_3 = pd.DataFrame(array_3).std(axis=0)

# Lists etc.
U_N = pd.DataFrame(df_2.values-df_3.values)
U_Nsd = pd.DataFrame(sd_2.values-sd_3.values)
mut_list = [('S','R'),('S','R'),('A','G'),('A','V'),('T','M'),('L','R'),('L','P'),('S','T')]
mut_loci = [152-150,153-150,166-150,166-150,182-150,183-150,183-150,205-150]

# Define simple model function -------------------------------------------------------------------------------
def ddG_H(wild_dsasa, mutant_loci, mut_list):
    
    # wild_sasa = list of SASA values to consider for mutation (float)
    # mut_list = list of wild type residue and mutant strings [(N,M)]
    
    #Empirical SASA Normalization values for Amino Acide from Tien et al. PLos ONE 2013 {A^2)  *0.01 nm^2/A^2
    A0_G = 97*0.01
    A0_A = 121*0.01
    A0_S = 143*0.01
    A0_R = 265*0.01
    A0_P = 154*0.01
    A0_T = 163*0.01
    A0_M = 203*0.01
    A0_V = 165*0.01
    A0_L = 191*0.01
    
    A0 = {'S': A0_S,'G': A0_G,'A': A0_A,'R': A0_R,'P': A0_P,'T': A0_T,'M': A0_M,'V': A0_V,'L': A0_L}
    
    # Relative Hydrophobicities from Eisenberg 1984 Table 1 in kcal/mol
    H_G = 0.48
    H_A = 0.62
    H_S = -0.18
    H_R = -2.53
    H_P = 0.12
    H_T =  -0.05
    H_M = 0.64
    H_V = 1.08
    H_L = 1.06
    
    H = {'S': H_S, 'G': H_G,'A': H_A,'R': H_R,'P': H_P,'T': H_T,'M': H_M,'V': H_V,'L': H_L}
    
    ddG = np.zeros(len(mut_list))
    
    for i,(j,k) in enumerate(mut_list):
        N = j
        M = k
        A0_N = A0[N]
        H_N = H[N]
        H_M = H[M]
        loci = mut_loci[i]
        print(loci,j,k)
        dA = wild_dsasa[loci]
        ddG[i] = (H_M-H_N)*(dA/A0_N)
        
    return ddG

# Plot ----------------------------------------------------------------------------------

import matplotlib.patches as mpatches
ddG_UN = ddG_H(np.concatenate(U_N.values),mut_loci,mut_list)
np.save('ddG_UN',ddG_UN)
mutations = ['S152R', 'S153R' , 'A166G', 'A166V', 'T182M', 'L183R','L183P', 'S205T']
fig,ax = plt.subplots(dpi=300)
cc =[]
for i,val in enumerate(ddG_UN):
    if val < 0:
        cc.append('blue')
    elif val >= 0:
        cc.append('red')
ax.bar(mutations,ddG_UN,color=cc)
ax.set_ylabel(r' $\Delta\Delta G_H^{U \rightarrow N}$ (kcal/mol)')
red_patch = mpatches.Patch(color='red', label='Destabilizing')
blue_patch = mpatches.Patch(color='blue', label='Stabilizing')
plt.legend(handles=[red_patch,blue_patch])
fig.savefig('HT_model.png')

