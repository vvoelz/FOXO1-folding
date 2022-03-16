# Imports
import matplotlib.pyplot as plt
import numpy as np
import pyemma
import os, sys

FES_dir = 'fes_plot'
its_dir = 'its'

print('Importing TICA data...')
distances_tica_output = np.load(f'{FES_dir}/distances_tica2.5_output.npy',allow_pickle=True)
distances_tica_concatenated = np.concatenate(distances_tica_output)

print('Importing TICA CLuster Centers...')
cluster_centers = np.load(f'{FES_dir}/distances_cluster_1000_4_centers.npy',allow_pickle=True)
cluster_dtrajs = list(np.load(f'{FES_dir}/distances_cluster_1000_4_dtrajs.npy',allow_pickle=True))
distances_dtrajs_concatenated = np.concatenate(np.load(f'{FES_dir}/distances_cluster_1000_4_dtrajs.npy',allow_pickle=True))

print('Importing the MSM object...')
states = np.load(f'{FES_dir}/states.npy',allow_pickle=True)
highest_membership = np.load(f'{FES_dir}/highest_membership.npy',allow_pickle=True)
coarse_state_centers = np.load(f'{FES_dir}/coarse_state_centers.npy',allow_pickle=True)
coarse_state_centers[0] = [np.mean(cluster_centers[states[0]][:,j]) for j in np.arange(4)]
print('Importing MFPTs...')
inverse_mfpt = np.load(f'{FES_dir}/inverse_mfpt.npy',allow_pickle=True)
mfpt = np.load(f'{FES_dir}/mfpt.npy',allow_pickle=True)
nstates = 3


# make the Plot

from matplotlib.gridspec import GridSpec

# create objects
fig = plt.figure(figsize=(7,6), constrained_layout=False, dpi=300)
gs = GridSpec(3, 3, figure=fig)
 
# create sub plots as grid
ax1 = fig.add_subplot(gs[0:2, :])   # TICA landscape
ax2 = fig.add_subplot(gs[2,0])      # its
ax3 = fig.add_subplot(gs[2,1])      # eigenvec 1
ax4 = fig.add_subplot(gs[2,2])      # eigenvec 2
 
### TICA landscape
if (True):
    markersize=[] 
    plist = np.load(f'{FES_dir}/p.npy',allow_pickle=True)
    for i, p in enumerate(plist):
        markersize.append(p * 10000)
    #fig, ax = plt.subplots(figsize=(10,5),dpi=300)

    pyemma.plots.plot_network(inverse_mfpt,pos=coarse_state_centers[:, :2],
        figpadding=0,
        state_labels=['I','N','U'],                      
        arrow_label_format='%.1f $\\mu$s',
        arrow_labels=None,
        arrow_curvature=4,
        show_frame=True, ax=ax1, xticks=True, yticks=True)
    
    pyemma.plots.plot_free_energy(
        *distances_tica_concatenated[:, :2].T,
        weights=np.load(f'{FES_dir}/weights.npy',allow_pickle=True),
        ax=ax1,
        legacy=False, alpha=0.5, zorder=0)

    ax1.set_xlim(-45, 55)
    ax1.set_ylim(-40, 55)
    #ax.scatter(coarse_state_centers[:,0],coarse_state_centers[:,1], s=markersize, c='yellow',edgecolors='black', alpha=0.5)
    ax1.scatter(*cluster_centers[:, :2].T, s=1, c='black',alpha=0.4)
    ax1.set_xlabel('IC1')
    ax1.set_ylabel('IC2')
    ax1.text(11,9,f'{np.round(mfpt[0][2],decimals=1)} $\\mu$s')
    ax1.text(-10,10,f'{np.round(mfpt[1][2],decimals=1)} $\\mu$s')
    ax1.text(10,-30,f'{np.round(mfpt[2][0],decimals=1)} $\\mu$s')
    ax1.text(-14,-18,f'{np.round(mfpt[2][1],decimals=1)} $\\mu$s')
    ax1.text(-18,2,f'{np.round(mfpt[1][0],decimals=1)} $\\mu$s')
    ax1.text(-36,-10,f'{np.round(mfpt[0][1],decimals=1)} $\\mu$s')

### implied timescales
if (True):
    lags = [5,10,50,100]
    # plt.rcParams.update({'font.size': 15})
    dtrajs = list(np.load(f'{its_dir}/distances_cluster_1000_4_dtrajs.npy',allow_pickle=True)) 
    lag = 50
    dt = 0.5
    its_file = f'{its_dir}/its_timescales_{lag}.npy'
    timescales = np.load(f'{its_dir}/its_timescales_{lag}.npy',allow_pickle=True)
    lags = np.load(f'{its_dir}/its_lagtimes_{lag}.npy',allow_pickle=True) 
    if not os.path.exists(its_file):
        its = pyemma.msm.its(dtrajs, lags=lag, nits=10)
        its.save(its_file, model_name=f'ITS_{lag}')
    else:
        its = np.load(its_file)
    colors = ['blue', 'red', 'green', 'cyan', 'purple', 'orange', 'violet']
    xmax = np.max(lags)
    srt = np.argsort(lags)
    dt = [dt] * 2
    units = ['ns']*2
    for i in range(len(colors)):
        ax2.plot(lags[srt] * dt[0], timescales[..., i][srt] * dt[1],color=colors[i % len(colors)])
#    pyemma.plots.plot_implied_timescales(its, units='ns', dt=0.5, ax=ax2) # Plot and save ITS plot
    ax2.plot(lags[srt] * dt[0], lags[srt] * dt[1], linewidth=2, color='black')
    ax2.set_xlim([1.0 * dt[0], xmax * dt[0]])
    ax2.fill_between(
        lags[srt] * dt[0], ax2.get_ylim()[0]*np.ones(len(lags))*dt[1], lags[srt] * dt[1],
        alpha=0.5, color='grey')
    # formatting
    ax2.set_xlabel('lag time / %s' % units[0])
    ax2.set_ylabel('timescale / %s' % units[1])
    ax2.set_yscale('log')
    ax2.set_box_aspect(1)
    ax2.set_yticks([10,100,1000,10000])

### eigenvecs 1 and 2
if (True):

    import matplotlib as mpl
    from mpl_toolkits.axes_grid1 import make_axes_locatable


    eigvec = np.load('eigvec_right.npy',allow_pickle=True)
    cbars = [False,False,False]
    # evec 1

    i = 0
    mappable = pyemma.plots.plot_contour(
        *distances_tica_concatenated[:, :2].T,
        eigvec[distances_dtrajs_concatenated, i + 1],
        ax=ax3,
        cbar=cbars[i],
        cax=ax3,
        cmap='PiYG',
        mask=True,
        norm=mpl.colors.TwoSlopeNorm(vmin=-1 , vcenter=0 ,vmax=1)
        )
    cbarticks = np.array([np.min(eigvec[distances_dtrajs_concatenated, i + 1]), 0, np.max(eigvec[distances_dtrajs_concatenated, i + 1])])
    cbarlabels = (['-', '0', '+'])
    ax3.set_xlabel('IC1')
    ax3.scatter(*cluster_centers[:, :2].T, s=20, c='white', edgecolors='black')
    ax3.set_ylabel('IC2')

    # evec 2

    i = 1
    mappable = pyemma.plots.plot_contour(
        *distances_tica_concatenated[:, :2].T,
        eigvec[distances_dtrajs_concatenated, i + 1],
        ax=ax4,
        cbar=cbars[i],
        cax=None,
        cmap='PiYG',
        mask=True,
        norm=mpl.colors.TwoSlopeNorm(vmin=-1 , vcenter=0 ,vmax=1)
        )
    #cbarticks = np.array([np.min(eigvec[distances_dtrajs_concatenated, i + 1]), 0, np.max(eigvec[distances_dtrajs_concatenated, i + 1])])
    #cbarlabels = (['-', '0', '+'])
    ax4.set_xlabel('IC1')
    ax4.scatter(*cluster_centers[:, :2].T, s=20, c='white', edgecolors='black')
    # ax4.set_ylabel('IC2')

    if cbars[i]:
        divider = make_axes_locatable(ax)
        cax = fig.add_axes([0.99, 0.2, 0.01, 0.65])#divider.append_axes("right", size="5%", pad=0.1)
        cax.yaxis.set_ticks(cbarticks)
        cax.yaxis.set_ticklabels(['-', '0', '+'])
    else:
        cax=None    

fig.tight_layout()
fig.savefig('Figure4_nolabels_test.png')


