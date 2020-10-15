import json
import matplotlib.pyplot as plt
import numpy as np

# load in the json output
# change to phonon_bands.json, the same script works for phonons
jfile = open('electron_bands.json')
data = json.load(jfile)

# unpack the json file
points = np.array(data['pathIndices'])
numBands = data['numBands']
energies = np.array(data['energies'])
mu = data['chemicalPotential']
pathLabels = data['highSymLabels']
pathTicks = data['highSymIndices']

# plot some vertical lines at high sym points
for i in pathTicks:
        plt.axvline(i, color='grey')

# plot chemical potential
if(data['particleType'] == "electron"):
        plt.axhline(0, color='grey', ls='--')

# plot the bands
for i in range(numBands):
        plt.plot(points, energies[:,i] - mu, lw=2, color='royalblue')

plt.xticks(pathTicks,pathLabels)
plt.ylabel(r'Energy [' + data['energyUnit'] +']')
plt.ylim(None, None)
plt.xlim(points[0],points[-1])

plt.show()
