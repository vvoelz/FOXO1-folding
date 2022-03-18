import numpy as np
import pyemma
import matplotlib.pyplot as plt

cluster_scores = np.load('FOXO1_cluster_scores.npy')
Vamp2_scores = np.load('FOXO1_VAMP2_scores.npy')
Vamp2_errors = np.load('FOXO1_VAMP2_errors.npy')
timestep=0.5
stride = 1
# Vamp2 Scores by Lag
lags = [5, 10, 20, 50]
dims = [i + 1 for i in range(10)]
fig, (ax1,ax2) = plt.subplots(2,1,figsize=(10,7),dpi=300)
for i, lag in enumerate(lags):
    scores = Vamp2_scores[i]
    errors = Vamp2_errors[i]
    color = 'C{}'.format(i)
    ax1.fill_between(dims, scores - errors, scores + errors, alpha=0.3, facecolor=color)
    ax1.plot(dims, scores, '--o', color=color, label=f'lag={lag*timestep*stride}ns')
ax1.legend()
ax1.set_aspect(1.0/ax1.get_data_ratio(), adjustable='box')
ax1.set_xlabel('number of dimensions')
ax1.set_ylabel('VAMP2 score')
#fig.tight_layout()
#fig.savefig(f'{name}_lag_pick.png')

# Vamp2 Scores for Clustering
n_clustercenters = [100,500,1000,2000,2500]
lower, upper = pyemma.util.statistics.confidence_interval(cluster_scores.T.tolist(), conf=0.9)
ax2.fill_between(n_clustercenters, lower, upper, alpha=0.3)
ax2.plot(n_clustercenters, np.mean(cluster_scores, axis=1), '-o')
ax2.semilogx()
ax2.set_aspect(1.0/ax2.get_data_ratio(), adjustable='box')
ax2.set_xlabel('number of cluster centers')
ax2.set_ylabel('VAMP-2 score')
ax1.text(0.91,0.05,'A',transform=ax1.transAxes,fontsize=20)
ax2.text(0.91,0.05,'B',transform=ax2.transAxes,fontsize=20)
fig.tight_layout()
fig.savefig(f'Figure_SI8.png')
"""
