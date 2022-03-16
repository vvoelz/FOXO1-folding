# Imports
import matplotlib.pyplot as plt
import numpy as np
import pyemma
import pandas as pd

# Import SASA
sasa = np.load('all_sasa.npy', allow_pickle=True)
tripep = pd.read_pickle('Tripep_SASA.pkl')

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

for i,array in enumerate(list_2):
    df = pd.DataFrame(array)
    mean = df.mean(axis=0)
    for j, res in enumerate(mean):
        array_2[i,j] = res
    df_2 = pd.DataFrame(array_2).mean(axis=0)
    
for i,array in enumerate(list_3):
    df = pd.DataFrame(array)
    mean = df.mean(axis=0)
    for j, res in enumerate(mean):
        array_3[i,j] = res
    df_3 = pd.DataFrame(array_3).mean(axis=0)

# Plot SASA Bar Plot for Mutations
mut_list = [152-150,153-150,166-150,182-150,183-150,205-150]
labels = ['S152','S153','A166','T182','L183','S205']
states = ['U_tp','U','I','N']
rel_sasa = []
for i in range(len(mut_list)):
    mut = [tripep[labels[i]].values[0],df_3[mut_list[i]],df_1[mut_list[i]],df_2[mut_list[i]]]
    rel_sasa_i = mut[0]/mut[1]
    rel_sasa.append(rel_sasa_i)
    print(f'{labels[i]} relative SASA (/nm^2)\n {states[0]}/{states[1]} = {mut[0]}/{mut[1]} = {rel_sasa_i}')
print(f'Mean = {np.mean(rel_sasa)}\nSD = {np.std(rel_sasa)}')
