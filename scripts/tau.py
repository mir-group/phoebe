import json
import matplotlib.pyplot as plt
import numpy as np

# load in the json output
jfile = open('rta_relaxation_times.json')
data = json.load(jfile)

# unpack the json file
energies = np.array(data['energies'])      # dimensions (nCalculation, dim, ik, ib)
tau = np.array(data['relaxationTimes'])    # dimensions (nCalculation, dim, ik, ib)
mu = np.array(data['chemicalPotentials'])
T = np.array(data['temperatures'])

# for now, let's select one calculation
# the index used to select the calculation index
# also corresponds to the index for the temperature
# and chemical potential of this calculation as stored
# in those arrays.
calcIndex = 0
energies = energies[calcIndex]
tau = tau[calcIndex]
mu = mu[calcIndex]
nbands = energies.shape[2]

# plot chemical potential
plt.figure(figsize=(5,5))
if(data['particleType'] == "electron"):
        plt.axvline(mu, color='grey', ls='--')

# plot the lifetimes
colors = plt.get_cmap('winter')(np.linspace(0,1,nbands))
for ib in range(nbands):
        plt.scatter(energies[:,:,ib] - mu, tau[:,:,ib], marker='x', s=18, alpha=0.25, color=colors[ib])

# plot aesthetics
plt.xlabel(r'Energy [' + data['energyUnit'] +']',fontsize=12)
plt.ylabel(r'$\tau_{' + data['particleType'] + '}$ [' + data['relaxationTimeUnit'] + ']',fontsize=12)
plt.xlim(np.min(energies),np.max(energies))
plt.ylim(None, np.max(tau)+np.max(tau)*0.1)

plt.show()
