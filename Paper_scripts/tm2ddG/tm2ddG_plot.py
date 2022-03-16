import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches

Tm_DSC = np.array((57.8,57.0,58.6,55.1,61.9,51.5,41.7,48.2,58.3))
labels = ['S152R','S153R','A166G','A166V','T182M','L183R','L183P','S205T']

def ddG_fit(T,Tm,N):
    ddG_T = np.zeros(len(Tm)-1)
    Tm0 = Tm[0]+273
    for i, Tm in enumerate(Tm[1::]):
        Tm = Tm+273
        dCp = -0.0139*N
        dHm = -0.698*N #0.698 for 60 degrees 1.26 for 100
        print(Tm,dCp,dHm)
        ddG_T[i] = T*(dCp*np.log(Tm/Tm0)+dHm*(Tm0**(-1)-Tm**(-1)))-dCp*(Tm-Tm0)
    return ddG_T
test = ddG_fit(298,Tm_DSC,100)
np.save('simpleddg',test)
fig,ax = plt.subplots(dpi=300)
cc =[]
for i,val in enumerate(test):
    if val < 0:
        cc.append('blue')
    elif val >= 0:
        cc.append('red')
ax.bar(labels,test,color=cc)
ax.set_ylabel(r' $\Delta\Delta G_{N_{res}}^{U \rightarrow N}$ (kcal/mol)')
red_patch = mpatches.Patch(color='red', label='Destabilizing')
blue_patch = mpatches.Patch(color='blue', label='Stabilizing')
plt.legend(handles=[red_patch,blue_patch])
fig.savefig('simpleg')

#----------------------------------------------------------------------------------------------------
dTm_DSC = [-0.8,0.8,-2.7,4.2,-6.3,-16.1,-9.6,0.5]

m_0, b_0 = np.polyfit(test, dTm_DSC, 1)


fig,ax = plt.subplots(dpi=300)
ax.set_ylabel('$\Delta T_m (^oC)$')
ax.set_xlabel(r'$\Delta\Delta G_{N_{res}}^{U \rightarrow N}$ (kcal/mol) ')
ax.set_title(r'$\Delta T_m^{DSC}$ v. $\Delta\Delta G_{N_{res}}^{U \rightarrow N}$')

ax.plot(test, m_0*test + b_0,c='black', alpha=0.3,linestyle=':')

markers = ['.','v','^','1','s','P','*','x']

correlation_matrix = np.corrcoef(test,dTm_DSC)
correlation_xy = correlation_matrix[0,1]
r2_DSC = np.round(correlation_xy**2,2)

ax.text(0.2,3,f'$\Delta T_m={np.round(m_0,2)}$'+r'$\Delta\Delta G_{N_{res}}^{U \rightarrow N}$'+f'+{np.round(b_0,2)}\n $r^2={r2_DSC}$',fontsize=12)
c_list = []
for i in test:
    if i > 0:
        c_list.append('r')
    if i < 0:
        c_list.append('b')

for i in range(len(test)):
    ax.scatter(test[i],dTm_DSC[i],color=c_list[i],marker=markers[i],label=labels[i])
    ax.legend()
fig.tight_layout()
fig.savefig('simple_G_comp')
