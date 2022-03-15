# Imports
import matplotlib.pyplot as plt
import numpy as np
import pyemma

# Import TICA data
distances_tica_output=np.load(f'distances_tica2.5_output.npy',allow_pickle=True)
distances_tica_concatenated = np.concatenate(distances_tica_output)

# Import TICA CLuster Centers
cluster_centers = np.load(f'distances_cluster_1000_4_centers.npy',allow_pickle=True)

#Import MSM
states = np.load('states.npy',allow_pickle=True)
highest_membership = np.load('highest_membership.npy',allow_pickle=True)
coarse_state_centers = np.load('coarse_state_centers.npy',allow_pickle=True)
coarse_state_centers[0] = [np.mean(cluster_centers[states[0]][:,j]) for j in np.arange(4)]

inverse_mfpt = np.load('inverse_mfpt.npy',allow_pickle=True)
mfpt = np.load('mfpt.npy',allow_pickle=True)
nstates = 3

#Plot
markersize=[]
p = np.load('p.npy',allow_pickle=True)
for i, p in enumerate(p):
    markersize.append(p * 10000)
fig, ax = plt.subplots(figsize=(10,5),dpi=300)
pyemma.plots.plot_network(inverse_mfpt,pos=coarse_state_centers[:, :2],
    figpadding=0,
    state_labels=['I','N','U'],                      
    arrow_label_format='%.1f $\\mu$s',
    arrow_labels=None,
    arrow_curvature=4,
    show_frame=True, ax=ax, xticks=True, yticks=True)
pyemma.plots.plot_free_energy(
    *distances_tica_concatenated[:, :2].T,
    weights=np.load('weights.npy',allow_pickle=True),
    ax=ax,
    legacy=False, alpha=0.5, zorder=0)
ax.set_xlim(-45, 55)
ax.set_ylim(-40, 55)
#ax.scatter(coarse_state_centers[:,0],coarse_state_centers[:,1], s=markersize, c='yellow',edgecolors='black', alpha=0.5)
ax.scatter(*cluster_centers[:, :2].T, s=1, c='black',alpha=0.4)
ax.set_xlabel('IC 1')
ax.set_ylabel('IC 2')
ax.text(11,9,f'{np.round(mfpt[0][2],decimals=1)} $\\mu$s')
ax.text(-10,10,f'{np.round(mfpt[1][2],decimals=1)} $\\mu$s')
ax.text(10,-30,f'{np.round(mfpt[2][0],decimals=1)} $\\mu$s')
ax.text(-14,-18,f'{np.round(mfpt[2][1],decimals=1)} $\\mu$s')
ax.text(-18,2,f'{np.round(mfpt[1][0],decimals=1)} $\\mu$s')
ax.text(-36,-10,f'{np.round(mfpt[0][1],decimals=1)} $\\mu$s')
fig.tight_layout()
fig.savefig('FOX01_FES.png')
