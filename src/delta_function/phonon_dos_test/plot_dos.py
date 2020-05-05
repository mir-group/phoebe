import numpy as np
import matplotlib.pyplot as plt
#import scipy.integrate as scint

# file paths?
plotfilepath = "./"
outfilepath = "./"

#phoebe tetra test
w_n = np.loadtxt(plotfilepath + 'phonon_dos',usecols=(0,))
dos_n = np.loadtxt(plotfilepath + 'phonon_dos',usecols=(1,))

#shengbte
w_s = np.loadtxt(plotfilepath + 'sheng.dos_wBAs22c',usecols=(0,))/2/np.pi #change to THz
dos_s = np.loadtxt(plotfilepath + 'sheng.dos_wBAs22c',usecols=(1,))*2*np.pi #change to 1/THz

plt.plot(w_s,dos_s,'b.',label='shengbte: adaptive gaussian')
plt.plot(w_n,dos_n,'r-',label='phoebe: tetrahedron')

plt.legend(loc='upper left', frameon=False)    
plt.xlabel('frequency (THz)',fontsize=14)
plt.ylabel('density of states (THz$^{-1}$)',fontsize=14)

plt.savefig(outfilepath + 'dos.pdf',dpi=72)
#plt.clf()
