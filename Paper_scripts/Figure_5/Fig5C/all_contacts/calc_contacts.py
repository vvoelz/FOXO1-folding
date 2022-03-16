import numpy as np
import pandas as pd
import pyemma
import mdtraj as md
import matplotlib.pyplot as plt
from itertools import combinations

mutations = [152,153,166,182,183,205]

tic_list = np.loadtxt('columns_tic.txt')[:,0]
traj = md.load('tic1.pdb')# Cluster Centers aligned along tIC1


states = np.load('states.npy',allow_pickle=True) # Get metastable macrostates

state_1 = states[0] # I
state_2 = states[1] # N
state_3 = states[2] # U

counts= [len(state_3),len(state_1),len(state_2)] # Count of microstates assigned to each macro
state_list = [state_3,state_1,state_2]

contacts = md.compute_contacts(traj,scheme='closest-heavy') # Calculate Cluster Center FNC

res1 = contacts[1][:,0]
res2 = contacts[1][:,1]
dist = contacts[0]*10
df = pd.DataFrame([res1,res2]+[i for i in dist],index=['Res1','Res2']+[f'CC_{i.astype(int)}' for i in tic_list]).T
df[df.columns[0]] = df[df.columns[0]].astype(int)
df[df.columns[1]] = df[df.columns[1]].astype(int)


con_1 =  df[[f'CC_{i.astype(int)}' for i in state_1]] # FNC I
con_2 = df[[f'CC_{i.astype(int)}' for i in state_2]] # FNC N
con_3 = df[[f'CC_{i.astype(int)}' for i in state_3]] # FNC U

sts_list = [con_3,con_1,con_2]

"""
For each microstate assigned to each macrostate, for each inter-residue contact, values  < 4.5 Angstrom are summed and divided by the total number of micostates assigned to each macro
"""
z1 = con_1[con_1 < 4.5].count(axis=1)
z2 = con_2[con_2 < 4.5].count(axis=1)
z3 = con_3[con_3 < 4.5].count(axis=1)

contact_sum = [z3,z1,z2]
contact_prob = [contact_sum[i]/counts[i] for i in range(3)]
x = res1+150
y = res2+150
probabilities = pd.DataFrame([x.astype(int),y.astype(int)]+contact_prob,index=['Res1','Res2','U','I','N']).T
probabilities.Res1 = probabilities.Res1.astype(int)
probabilities.Res2 = probabilities.Res2.astype(int)
probabilities.to_pickle('probabilites.pkl')

