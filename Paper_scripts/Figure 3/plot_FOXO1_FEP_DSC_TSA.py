import numpy as np
import matplotlib
from matplotlib import pyplot as plt

import pandas as pd

### Set Arial font fomily

# first method
# matplotlib.rcParams['font.family'] = ['Source Han Sans TW', 'sans-serif']

# second method
print(matplotlib.rcParams.keys())

matplotlib.rcParams['font.family'] = ['sans-serif']
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica'] # 'Source Han Sans TW', ...]

# read in the data
df = pd.read_csv('FOXO1_FEP_DSC_TSA.csv')
print(df)

# ddG = df['ddG(kJ/mol)']/4.184  # convert to kcal mol^{-1}
ddG = df['ddG_FEP(kcal/mol)'] # from Table 4 of paper 1/08/2022!
ddG_sigma = df['sigma']

dTm_TSA = df['TSA(C)']   # in C 
dTm_DSC = df['DSC(C)']   # in C 

plt.figure(figsize=(3.25, 2.75))
plt.plot(ddG, dTm_DSC, 'rs', label='DSC')
plt.plot(ddG, dTm_TSA, 'ko', label='TSA')
for i in range(len(ddG)):
    plt.plot([ddG[i]-ddG_sigma[i], ddG[i]+ddG_sigma[i]], [dTm_TSA[i]]*2, 'k-')
plt.xlim(-2, 8)
plt.ylim(-20,10)
plt.xlabel('$\Delta \Delta G_{FEP}$ (kcal mol$^{-1}$)')
plt.ylabel('$\Delta T_m$ ($^{\circ}$C)')
plt.legend(loc='best')
plt.tight_layout()

plt.savefig('Figure3B.pdf')
print('...wrote Figure3B.pdf.')








