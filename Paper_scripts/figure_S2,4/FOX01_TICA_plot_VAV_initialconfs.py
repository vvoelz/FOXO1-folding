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



### Find the indices for the first CLONE of each RUN
traj_indices = 20*[None]    # index for each RUN
for i in range(len(names)):
    fields = names[i].split('_')
    run = int(fields[0].replace('RUN',''))
    if traj_indices[run] == None:
        traj_indices[run] = i


### Make the plot

plt.figure( figsize=(4, 4) )

# plot the contours of the tICA landscapes 
cs = plt.contour(X, Y, Z, levels=[1e-5, 1e-4, 0.001, 0.005, 0.01, 0.05, 0.1],
                 colors=['#808080', '#A0A0A0', '#C0C0C0'], extend='both')

# Plot the first frame of each of these trajectories in different colors and markers
markers = ["o",     # m02     circle
           "v",     # m03	triangle_down
           "^",     # m04	triangle_up
           "<",     # m05	triangle_left
           ">",     # m06	triangle_right
#           "1",     # m07	tri_down
#           "2",     # m08	tri_up
#           "3",     # m09	tri_left
#           "4",     # m10	tri_right
           "8",     # m11	octagon
           "s",     # m12	square
           "p",     # m13	pentagon
           "P",     # m23	plus (filled)
           "*",     # m14	star
           "h",     # m15	hexagon1
           "H",     # m16	hexagon2
           "+",     # m17	plus
           "x",     # m18	x
           "X",     # m24	x (filled)
           "D" ]    # m19	diamond
markers = markers*2    # recycle the markers

for run in range(20):

    cmap = plt.cm.rainbow
    try:
        i =  traj_indices[run]
        # plot a rainbow traj on top of the contour background
        traj = distances_tica_output[i]  # shape (nsnapshots, nTICs=4)
        #nframes, ntics = traj.shape
        plt.plot(traj[0][0], traj[0][1], lw=0, marker=markers[run], color=cmap(run/20.0), label='RUN '+str(run) )
    except:
        print('...Failed.')

plt.xlabel('tIC1')
plt.ylabel('tIC2')
plt.legend(loc='best', ncol=2, fontsize='xx-small', mode='marker')

plt.tight_layout()
outfile = 'initial_confs.pdf'
plt.savefig(outfile)
print('Wrote:', outfile)
plt.close()


sys.exit(1)

