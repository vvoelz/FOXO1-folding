import numpy as np
from matplotlib import pyplot as plt
import sklearn.metrics
import sys

labels = ['S152R', 'S153R' , 'A166G', 'A166V', 'T182M', 'L183R','L183P', 'S205T']
labels2 = labels
dTm_DSC = [-0.8,0.8,-2.7,4.2,-6.3,-16.1,-9.6,0.5]

Use_OLD_Values = False
if (Use_OLD_Values):
    FEP = [-3.7,-5.8,2.6,-3.6,-1.3,25.7,28.6,7.3]
    FEP_err = [6.5,1.5,0.4,0.5,0.4,0.9,4.2,7.3]
    FEP = np.array([i/4.184 for i in FEP])
    FEP_err = np.array([i/4.184 for i in FEP_err])

else:
    FEP = np.array([0.19, -1.27, 0.79, -0.67, -0.33, 7.04, 7.07, 1.79])
    FEP_err = np.array([0.1528, 0.1815, 0.0239, 0.0263, 0.0621, 0.129, 0.1075, 0.0287])


sim = np.load('ddG_UN.npy',allow_pickle=True)
sim_err = np.load('ddG_err.npy',allow_pickle=True)
exp = np.load('simpleddg.npy',allow_pickle=True) 
exp_err = np.load('simpleddg_err.npy',allow_pickle=True)
 
FOLDX = np.load('foldx_ddg.npy',allow_pickle=True)

m_0, b_0 = np.polyfit(exp, sim, 1)
m_1, b_1 = np.polyfit(exp, FEP, 1)
m_2, b_2 = np.polyfit(exp, FOLDX, 1)


plt.figure(figsize=(11, 3.75))


markers = ['.','v','^','1','s','P','*','x']

RMSE_sim = np.sqrt(sklearn.metrics.mean_squared_error(sim, exp))
RMSE_sim = np.round(RMSE_sim,2)
RMSE_FEP = np.sqrt(sklearn.metrics.mean_squared_error(FEP, exp))
RMSE_FEP = np.round(RMSE_FEP,2)
RMSE_FOLDX = np.sqrt(sklearn.metrics.mean_squared_error(FOLDX, exp))
RMSE_FOLDX = np.round(RMSE_FOLDX,2)

"""
ax.text(0,-1.5,r'$\Delta \Delta G_{Exp.} = '+f'{np.round(m_0,2)}'+ r'\Delta\Delta G_{Sim.} + '+ f'{np.round(b_0,2)} : RMSE={RMSE_sim}$\n'+ r'$\Delta \Delta G_{Exp.} = '+ 
f'{np.round(m_1,2)}'+ r'\Delta\Delta G_{FEP}$ + '+f'{np.round(b_1,2)} : $RMSE={RMSE_FEP}$\n'
r'$\Delta \Delta G_{Exp.} = '+f'{np.round(m_2,2)}'+ r'\Delta\Delta G_{FOLDX} + '+ f'{np.round(b_2,2)} : RMSE={RMSE_FOLDX}$',fontsize='xx-small')
"""

##### Simulation vs. exp

plt.subplot(1,3,1)

for i in range(len(sim)):
    plt.scatter(exp[i], sim[i], color='purple',marker=markers[i],label=labels[i])
    plt.errorbar(exp[i], sim[i], color='purple', xerr=exp_err[i], yerr=sim_err[i], alpha=0.5, capsize=3)  
    # plt.legend(loc="lower center", fontsize='x-small',bbox_to_anchor=(0.5, -0.4), mode='marker',ncol = 4)
plt.legend(loc='best', fontsize='x-small', mode='marker')
plt.text(1.5, -0.5, f'RMSE={RMSE_sim}\nkcal/mol')

    
line = np.arange(-2,8)
plt.plot(line, m_0*line + b_0,c='purple',linestyle=':',label='Sim.',alpha=0.3)
plt.plot(line, line, 'k-')
plt.xlim(-1,3)
plt.ylim(-1,3)
#plt.yticks(np.arange(-2, 8, 2))
plt.xlabel(r'$\Delta\Delta G_{exp}$ (kcal/mol)')
plt.ylabel(r'$\Delta\Delta G_{pred}$ (kcal/mol) ')
plt.title('Simulation + HT model')
plt.tight_layout()

############## FOLDX vs. exp
plt.subplot(1,3,2)

for i in range(len(FOLDX)):
    plt.scatter(exp[i], FOLDX[i], color='blue',marker=markers[i],label=labels[i])
    plt.errorbar(exp[i], FOLDX[i], color='green', xerr=exp_err[i], alpha=0.5, capsize=3)
#plt.legend(loc="lower center", fontsize='x-small',bbox_to_anchor=(0.5, -0.4), mode='marker',ncol = 4)
plt.legend(loc='best', fontsize='xx-small', mode='marker')
plt.text(1.5, -0.5, f'RMSE={RMSE_FOLDX}\nkcal/mol')

line = np.arange(-2,8)
plt.plot(line, m_2*line + b_2,c='blue',linestyle=':',label='FOLDX',alpha=0.3)
plt.plot(line, line, 'k-')
plt.xlim(-1,3)
plt.ylim(-1,3)
plt.xlabel(r'$\Delta\Delta G_{exp}$ (kcal/mol)')
plt.ylabel(r'$\Delta\Delta G_{pred}$ (kcal/mol) ')
plt.title('FOLDX')

plt.tight_layout()


##### FEP vs. exp

plt.subplot(1,3,3)

for i in range(len(FEP)):
    plt.scatter(exp[i], FEP[i], color='green', marker=markers[i], label=labels[i])
    plt.errorbar(exp[i], FEP[i], color='green', xerr=exp_err[i], yerr=FEP_err[i], alpha=0.5, capsize=3)
#plt.legend(loc="lower center", fontsize='x-small',bbox_to_anchor=(0.5, -0.4), mode='marker',ncol = 4)
plt.legend(loc='best', fontsize='x-small', mode='marker')
plt.text(1.5, -1.0, f'RMSE={RMSE_FEP}\nkcal/mol')

line = np.arange(-2,8)
plt.plot(line, m_1*line + b_1,c='green',linestyle=':',label='FEP',alpha=0.3)
plt.plot(line, line, 'k-')
plt.xlim(-1,3)
plt.ylim(-2,8)
#plt.yticks(np.arange(-2, 8, 2))
plt.xlabel(r'$\Delta\Delta G_{exp}$ (kcal/mol)')
plt.ylabel(r'$\Delta\Delta G_{pred}$ (kcal/mol) ')
plt.title('FEP')
plt.tight_layout()



# ax.set_box_aspect(1)

# fig.subplots_adjust(bottom=0.25)
plt.tight_layout()

outfile = 'three_panel_dG_comp.pdf'
plt.savefig('three_panel_dG_comp.pdf')

print('Wrote:', outfile)

