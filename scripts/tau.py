import json
import matplotlib.pyplot as plt
import numpy as np

# load in the json output
jfile = open('rta_relaxation_times.json')
data = json.load(jfile)

# unpack the json file
energies = np.array(data['energies'])      # dimensions (iCalc, ik, ib)
nbands = energies.shape[2]
tau = np.array(data['relaxationTimes'])    # dimensions (iCalc, ik, ib, idim)
mu = np.array(data['chemicalPotentials'])
T = np.array(data['temperatures'])

# Could also load in group velocities, or wavevectors
#vels = np.array(data['velocities'])         # dimensions: (iCalc, ik, ib, dim)
#wavevectors = np.array(data['wavevectorCoordinates']) # dimensions: (ik, dim)

# for now, let's select one calculation
# the index used to select the calculation
# also corresponds to the index for the temperature
# and chemical potential of that calculation as stored
# in those arrays.
calcIndex = 0
energies = energies[calcIndex]
tau = tau[calcIndex]
mu = mu[calcIndex]
print("Calculation Temperature: ", T[calcIndex])

# plot the lifetimes, colored by band, for all dimensions
plt.figure(figsize=(5,5))
colors = plt.get_cmap('winter')(np.linspace(0,1,nbands))
for ib in range(nbands):
        for idim in range(tau.shape[2]):
                plt.scatter(energies[:,ib] - mu, tau[:,ib,idim], marker='x', s=18, alpha=0.25, color=colors[ib])

# plot aesthetics
plt.yscale('log')
plt.xlabel(r'Energy [' + data['energyUnit'] +']',fontsize=12)
plt.ylabel(r'$\tau_{' + data['particleType'] + '}$ [' + data['relaxationTimeUnit'] + ']',fontsize=12)
plt.xlim(0,None)
plt.ylim(1e3, np.max(tau)+np.max(tau)*0.1)

plt.tight_layout()
plt.show()
