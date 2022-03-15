import numpy as np

# All Tm_WT and Tm_mut are from DSC measurements
Tm_WT_C  = 57.8     # C
Tm_WT    = Tm_WT_C + 273.15   # convert to K
dTm_WT   = 0.3      # uncertainty of Tm_WT

Tm_mut_labels = ['S152R', 'S153R', 'A166G', 'A166V', 'T182M', 'L183R', 'L183P', 'S205T'] 
Tm_mut_C  = np.array((57.0, 58.6, 55.1, 61.9, 51.5, 41.7, 48.2, 58.3))
Tm_mut    = Tm_mut_C + 273.15  # convert to K
dTm_mut   =  1.0    # uncertainty of Tm_mut


def ddG_fold_uncertainty(Tm_WT, dTm_WT, Tm_mut, dTm_mut, N=100, T=298.15, debug=True):
    """Return an estimate of the uncertainty of ddG_fold in kcal/mol.

    INPUTS
    Tm_WT      WT melting temp in K
    dTm_WT     experimental uncertainty in Tm_WT
    Tm_mut     np.array of mutant melting temps in K
    dTm_mut    experimental uncertainty in Tm_mut

    PARAMETERS
    N          number of residues in the protein
    T          temperature for ddG_fold (Default: 298.15 K std temp)
    """

    # Compute estimates of dC_p and dH_fold based on
    # Robertson and Murphy (2013) https://doi.org/10.1016/j.abb.2012.09.008
    dC_p    = -0.0139*N    # kcal mol^{-1} K^{-1}
    dH_fold = -0.698*N     # kcal mol^{-1}


    term1  =  ( -dH_fold/T*(T/Tm_WT)**2 + dC_p*(1.0 - T/Tm_WT) )**2  * dTm_WT**2
    if debug:
        print('term1', term1)   # should be a scalar
 
    term2  =  ( dH_fold/T*(T/Tm_mut)**2 - dC_p*(1.0 - T/Tm_mut) )**2  * dTm_mut**2
    if debug:
        print('term2', term2)   # should be an array

    result = np.sqrt( term1 + term2 )
    return result

                     
#########

result = ddG_fold_uncertainty(Tm_WT, dTm_WT, Tm_mut, dTm_mut)
print('result', result)
print()

ddG_fold_values = [0.04, -0.04, 0.16, -0.17, 0.43, 1.51, 0.73, -0.03]  # in kcal/mol
Tm_mut_labels = ['S152R', 'S153R', 'A166G', 'A166V', 'T182M', 'L183R', 'L183P', 'S205T']

print('mut\tddG_fold(kcal/mol)l\t+/-uncertainty')
for i in range(len(Tm_mut_labels)):
    print(f'{Tm_mut_labels[i]}\t{ddG_fold_values[i]}\t{result[i]}')

outfile = 'simpleddg_err.npy'
np.save(outfile, result)
print('Wrote:', outfile)



