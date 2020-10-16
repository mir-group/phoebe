import json
import matplotlib.pyplot as plt
import numpy as np

# load in the json output
# change to phonon_bands.json, the same script works for phonons
jfile = open('electron_bands.json')
data = json.load(jfile)

# unpack the json file
points = np.array(data['wavevectorIndices'])
numBands = data['numBands']
energies = np.array(data['energies'])
pathLabels = data['highSymLabels']
pathTicks = data['highSymIndices']

# plot some vertical lines at high sym points
plt.figure(figsize=(5.5,5))
for i in pathTicks:
        plt.axvline(i, color='grey')

# if fermiLevel as set in the input file,
# we can also read it in and plot it
if('fermiLevel' in data):
        mu = data['fermiLevel']
        plt.axhline(0, color='grey', ls='--')
        energyLabel = r'E-E$_F$ [' + data['energyUnit'] +']'
# if it wasn't set, we won't subtract
# anything from the band energies
else:
        mu = 0
        energyLabel = r'Energy [' + data['energyUnit'] +']'

# plot the bands
for i in range(numBands):
        plt.plot(points, energies[:,i] - mu, lw=2, color='royalblue')

plt.xticks(pathTicks,pathLabels,fontsize=12)
plt.yticks(fontsize=12)
plt.ylabel(energyLabel,fontsize=14)
plt.ylim(None, None)
plt.xlim(points[0],points[-1])

plt.show()
