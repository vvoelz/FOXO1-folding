# Imports
import matplotlib.pyplot as plt
import tqdm 
import numpy as np
import pyemma
import os, sys

names = np.load('traj_list.npy')
names = [os.path.splitext(os.path.basename(i))[0] for i in names]
# print('names', names)

# Import TICA data
distances_tica_output=np.load(f'distances_tica2.5_output.npy',allow_pickle=True)

distances_tica_concatenated = np.concatenate(distances_tica_output)
print('distances_tica_concatenated.shape', distances_tica_concatenated)


# make a contour plot for the tICA landscape
bins = np.arange(-40, 50, 1)
zbins = np.arange(-40, 51, 1)
X, Y = np.meshgrid(bins, bins)
Z, edges1, edges2 = np.histogram2d(distances_tica_concatenated[:,0], distances_tica_concatenated[:,1],
                                   bins=zbins, normed=True)
Z = Z.transpose()



### Find the indices for EXAMPLE trajectories
U_to_U = 'RUN1_CLONE110'
U_to_I = 'RUN12_CLONE451'  # runner up: 'RUN12_CLONE54' 
U_to_N = None
I_to_U = 'RUN0_CLONE67'
I_to_I = 'RUN0_CLONE223'
I_to_N = 'RUN0_CLONE21'
N_to_U = 'RUN18_CLONE48'
N_to_I = 'RUN18_CLONE80'  # runners up: 'RUN18_CLONE429', RUN18_CLONE167
N_to_N = 'RUN3_CLONE8'

panel_labels  = ['$\\bf{U}$ $\\rightarrow$ $\\bf{U}$', '$\\bf{U}$ $\\rightarrow$ $\\bf{I}$', None]
panel_labels += ['$\\bf{I}$ $\\rightarrow$ $\\bf{U}$', '$\\bf{I}$ $\\rightarrow$ $\\bf{I}$', '$\\bf{I}$ $\\rightarrow$ $\\bf{N}$']
panel_labels += ['$\\bf{N}$ $\\rightarrow$ $\\bf{U}$', '$\\bf{N}$ $\\rightarrow$ $\\bf{I}$', '$\\bf{N}$ $\\rightarrow$ $\\bf{N}$']

# panel_labels = ['X->Y']*9
traj_indices = [None]*9

for i in range(len(names)):

    if names[i] == U_to_U:
        traj_indices[0] = i
    elif names[i] == U_to_I:
        traj_indices[1] = i

    elif names[i] == I_to_U:
        traj_indices[3] = i
    elif names[i] == I_to_I:
        traj_indices[4] = i
    elif names[i] == I_to_N:
        traj_indices[5] = i

    elif names[i] == N_to_U:
        traj_indices[6] = i
    elif names[i] == N_to_I:
        traj_indices[7] = i
    elif names[i] == N_to_N:
        traj_indices[8] = i


print(traj_indices)

### Make the plot

plt.figure( figsize=(10, 10) )

cmap = plt.cm.rainbow

for j in range(len(traj_indices)):

    i = traj_indices[j] 
    panel = j+1

    if i != None:

        print(f'Plotting {names[i]} in panel {panel}...')

        plt.subplot(3,3,panel)

        # plot the contours of the tICA landscapes 
        cs = plt.contour(X, Y, Z, levels=[1e-5, 1e-4, 0.001, 0.005, 0.01, 0.05, 0.1],
                 colors=['#808080', '#A0A0A0', '#C0C0C0'], extend='both')

        # plot a rainbow traj on top of the contour background
        traj = distances_tica_output[i]  # shape (nsnapshots, nTICs=4)
        nframes, ntics = traj.shape
        for t in range(nframes-1):
            plt.plot([traj[t][0], traj[t+1][0]], [traj[t][1], traj[t+1][1]], color=cmap(t/nframes) )
        plt.xlabel('tIC1')
        plt.ylabel('tIC2')
        plt.title(names[i] + f' ({0.5*nframes} ns)')
        plt.text(-30,30, panel_labels[j], fontsize=18)

#plt.legend(loc='best', ncol=2, fontsize='xx-small', mode='marker')

plt.tight_layout()
outfile = 'example_trajectories.pdf'
plt.savefig(outfile)
print('Wrote:', outfile)
plt.close()


