import json
import matplotlib.pyplot as plt
import numpy as np

# load in the json output
#jfile = open('phonon_dos.json') # works for phonons too
jfile = open('electron_dos.json')
data = json.load(jfile)

# unpack the json file
energies = np.array(data['energies'])
dos = np.array(data['dos'])
mu = data['chemicalPotential']

# plot chemical potential
plt.figure(figsize=(5,3))
if(data['particleType'] == "electron"):
        plt.axvline(0, color='grey', ls='--')

# plot the dos
plt.plot(energies - mu, dos, lw=2, color='royalblue')

# plot aesthetics
plt.xlabel(r'Energy [' + data['energyUnit'] +']')
plt.ylabel(r'DoS [' + data['dosUnit'] + ']')
plt.xlim(np.min(energies),np.max(energies))
plt.ylim(None, np.max(dos)+np.max(dos)*0.1)

plt.show()
