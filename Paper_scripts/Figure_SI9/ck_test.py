import numpy as np
import pyemma
import matplotlib.pyplot as plt

cktest = pyemma.load('cktest.h5')
# cktest
pyemma.plots.plot_cktest(cktest, dt=0.5, units='ns')
plt.savefig('FOX01_cktest.png')
"""
