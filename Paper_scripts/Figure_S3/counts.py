# Imports
import os
import pickle
import matplotlib.pyplot as plt
import numpy as np
import pyemma
import pandas as pd


names = np.load('traj_list.npy')
names = [os.path.splitext(os.path.basename(i))[0] for i in names]

if os.path.exists('analyzed.txt'):
    with open('analyzed.txt', "r") as f:
        analyzed = f.read()
else:
    analyzed = []
# Get Metastable States
states = np.load('states.npy',allow_pickle=True)
I = states[0]
N = states[1]
U = states[2]

dtrajs = np.load('distances_cluster_1000_4_dtrajs.npy',allow_pickle=True)
transition_labels = ['U-I','U-N','I-N','I-U','N-I','N-U']

if os.path.exists("transitions_dict.pkl"):
    with open("transitions_dict.pkl", "rb") as f:
        transitions = pickle.load(f)
else:
    transitions = {}
    for label in transition_labels:
        transitions[f'{label}_trajs'] = []
        transitions[f'{label}_counts'] = 0

for i,traj in enumerate(dtrajs):
    if names[i] in analyzed:
        print(f'{names[i]} Already Analyzed...skipping.')
        pass
    else:
        print(f'{names[i]}')
        for j,cc in enumerate(traj[1:]):
            if cc in U:
                if traj[1:][j-1] in U:
                    pass
                elif traj[1:][j-1] in I:
                    transitions['U-I_counts']+= 1
                    if names[i] not in transitions['U-I_trajs']: 
                        transitions['U-I_trajs'].append(names[i])
                elif traj[1:][j-1] in N:
                    transitions['U-N_counts'] += 1
                    if names[i] not in transitions['U-N_trajs']: 
                        transitions['U-N_trajs'].append(names[i])
            elif cc in I:
                if traj[1:][j-1] in I:
                    pass
                elif traj[1:][j-1] in N:
                    transitions['I-N_counts'] += 1
                    if names[i] not in transitions['I-N_trajs']: 
                        transitions['I-N_trajs'].append(names[i])
                elif traj[1:][j-1] in U:
                    transitions['I-U_counts'] += 1
                    if names[i] not in transitions['I-U_trajs']: 
                        transitions['I-U_trajs'].append(names[i])
            elif cc in N:
                if traj[1:][j-1] in N:
                    pass
                elif traj[1:][j-1] in I:
                    transitions['N-I_counts'] += 1
                    if names[i] not in transitions['N-I_trajs']: 
                        transitions['N-I_trajs'].append(names[i])
                elif traj[1:][j-1] in U:
                    transitions['N-U_counts'] += 1
                    if names[i] not in transitions['N-U_trajs']: 
                        transitions['N-U_trajs'].append(names[i])

        with open('analyzed.txt', "a+") as f:
            f.write(f'{names[i]}\n')

        with open("transitions_dict.pkl","wb") as f:
            pickle.dump(transitions,f)

counts = [transitions[f'{label}_counts']for label in transition_labels]
fig,ax = plt.subplots(dpi=300)
ax.bar(transition_labels,counts)
ax.set_yscale('log')
ax.set_ylabel('Transition Counts')
fig.savefig('counts.png')
plt.close()
