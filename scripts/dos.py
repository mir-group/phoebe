#!/usr/bin/env python3
import json
import matplotlib.pyplot as plt
import numpy as np
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot Density of States from"
                                     " a Json file generated by Phoebe")
    parser.add_argument("INPUT",
                        help="Name of the JSON file with the DoS")
    args = parser.parse_args()
    jfileName = args.INPUT

    # load in the json output
    with open(jfileName) as jfile:
        data = json.load(jfile)

    # unpack the json file
    try:
        dos = np.array(data['dos'])
    except KeyError:
        raise KeyError("DoS not found."
                       "Are you using the correct input json file?")

    energies = np.array(data['energies'])

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
    plt.xlim(np.min(energies-mu),np.max(energies-mu))
    plt.ylim(None, np.max(dos)+np.max(dos)*0.1)

    plt.tight_layout()
    plotFileName = "./" + jfileName.rstrip(".json")+".pdf"
    plt.savefig(plotFileName)
    plt.show(block=False)
