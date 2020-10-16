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

# if fermiLevel as set in the input file,
# we can also read it in and plot it
plt.figure(figsize=(6,2))
if('fermiLevel' in data):
        mu = data['fermiLevel']
        plt.axvline(0, color='grey', ls='--')
        energyLabel = r'E-E$_F$ [' + data['energyUnit'] +']'
# if it wasn't set, we won't subtract
# anything from the band energies
else:
        mu = 0
        energyLabel = r'Energy [' + data['energyUnit'] +']'

# plot the bands
plt.plot(energies - mu, dos, lw=2, color='royalblue')

# plot aesthetics
plt.xlabel(energyLabel,fontsize=12)
plt.ylabel(r'DoS [' + data['dosUnit'] + ']',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(np.min(energies),np.max(energies))
plt.ylim(None, np.max(dos)+np.max(dos)*0.1)

plt.tight_layout()
plt.show()
