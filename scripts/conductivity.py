import json
import matplotlib.pyplot as plt
import numpy as np

# load in the json output
jfile = open('rta_phonon_thermal_cond.json')
data = json.load(jfile)

# unpack the json file
kappa = np.array(data['thermalConductivity']) # dimensions: (temperature, dim1, dim2)
T = np.array(data['temperatures'])

# plot the thermal conductivity (here we just plot xx component)
plt.plot(T, kappa[:,0,0], lw=2, mew=1.5, ms=8, marker='x', color='royalblue')

plt.xlabel('Temperature [' + data['temperatureUnit'] + ']',fontsize=12)
plt.ylabel(r'$\kappa_\mathrm{' + data['particleType'] + '}$ [' + data['thermalConductivityUnit'] +']',fontsize=12)
plt.ylim(None, np.max(kappa) + np.max(kappa)*0.1)
plt.xlim(None, np.max(T) + np.max(T)*0.1)

plt.tight_layout()
plt.show()
