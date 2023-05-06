#!/usr/bin/env python3
import json
import matplotlib.pyplot as plt
import numpy as np
import argparse

if __name__ == "__main__":

    # load in the json output for the phonon thermal conductivity
    jfileName = "rta_onsager_coefficients.json"
    with open(jfileName) as jfile:
        data = json.load(jfile)

    try:
        kappa_e = np.array(data['electronicThermalConductivity'])
        # size: (temperature, dim1, dim2)
    except KeyError:
        raise KeyError("Phel Phonon thermal conductivity not found.")


    # load in the json output for the phonon thermal conductivity
    jfileName = "rta_phonon_thermal_cond.json"
    with open(jfileName) as jfile:
        data = json.load(jfile)

    try:
        kappa_ph = np.array(data['thermalConductivity'])
        # size: (temperature, dim1, dim2)
    except KeyError:
        raise KeyError("Phonon thermal conductivity not found.")

    T = np.array(data['temperatures'])

    # plot the thermal conductivity (here we just plot xx component)
    plt.figure(figsize=(5,4))
    plt.plot(T, kappa_ph[:,0,0], lw=2, mew=1.5, ms=8, marker='x', color='grey',label="$\kappa_{ph}$")
    plt.plot(T, kappa_e[:,0,0], lw=2, mew=1.5, ms=8, marker='x', ls='--', color='grey',label="$\kappa_{e}$")
    plt.plot(T, kappa_ph[:,0,0] + kappa_e[:,0,0], lw=2, mew=1.5, ms=8, marker='x', ls='--', color='mediumpurple',label="$\kappa_{total}$")

    plt.legend(frameon=False)
    plt.xlabel('Temperature [' + data['temperatureUnit'] + ']',fontsize=12)
    plt.ylabel(r'$\kappa$ [' + data['thermalConductivityUnit'] +']',fontsize=12)
    #plt.ylim(None, np.max(kappa) + np.max(kappa)*0.1)
    plt.xlim(None, np.max(T) + np.max(T)*0.1)

    plt.tight_layout()
    plotFileName = "kappa_tot.png"
    plt.savefig(plotFileName,dpi=200)
