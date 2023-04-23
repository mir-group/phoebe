#!/usr/bin/env python3
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

calcIndex = 0     #

colors = ["green", "red"]

for idx, i in enumerate(['N','U']):

    # load in the json output
    jfileName = "rta_ph_relaxation_times.json"
    with open(jfileName) as jfile:
        data = json.load(jfile)

    particleType = data['particleType']

    # unpack the json file
    if i == 'U':
        tau = np.array(data['umklappRelaxationTimes'])    # dimensions (iCalc, ik, ib)
    elif i == 'N':
        tau = np.array(data['normalRelaxationTimes'])    # dimensions (iCalc, ik, ib)

    tau[np.where(tau==None)] = 0   # remove None values (from gamma pt acoustic ph)
    energies = np.concatenate(data['energies'][calcIndex]).ravel()
    T = np.array(data['temperatures'])

    tau = np.concatenate(tau[calcIndex]).ravel() #np.asarray(tau[calcIndex].flatten())

    # plot the lifetimes
    plt.scatter(energies, tau, marker='o', s=10, alpha=0.25, label=i,color=colors[idx])

# plot aesthetics
plt.yscale('log')
plt.xlabel(r'Energy [' + data['energyUnit'] +']',fontsize=12)
units = ' [' + data['relaxationTimeUnit'] + ']'
plt.ylabel(r'$\tau_{ph}$ ' + units, fontsize=12)

# set limits
plt.tight_layout()
plt.ylim(1000, 1000000000)
#plt.xlim(0,60)
plt.legend(frameon=False)
plt.savefig("UN.png")


