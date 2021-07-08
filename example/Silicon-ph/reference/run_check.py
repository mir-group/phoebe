#!/usr/bin/env python3

import glob
import os
import json
import sys
import numpy

if __name__ == "__main__":

    listOfJsons = glob.glob("*.json")
    
    for filename in listOfJsons:

        filename2 = os.path.join("reference", filename)
    
        with open(filename) as f1:
            data1 = json.load(f1)
        with open(filename2) as f2:
            data2 = json.load(f2)

        print(filename)
        print(filename2)
        print(" ")
            
        if "thermal_cond" in filename:
            k1 = numpy.array(data1['thermalConductivity'])
            k2 = numpy.array(data2['thermalConductivity'])
            diff = ((k1-k2)**2).sum()
            if diff > 0.25:
                print(diff)
                print(filename)
                sys.exit(1)

        if "specific_heat" in filename:
            k1 = numpy.array(data1['specificHeat'])
            k2 = numpy.array(data2['specificHeat'])
            diff = (k1-k2)**2
            if diff > 0.1:
                print(diff)
                print(filename)
                sys.exit(1)

        if "specific_heat" in filename:
            k1 = numpy.array(data1['specificHeat'])
            k2 = numpy.array(data2['specificHeat'])
            diff = (k1-k2)**2
            if diff > 0.1:
                print(diff)
                print(filename)
                sys.exit(1)

        if "_bands." in filename:
            k1 = numpy.array(data1['energies'])
            k2 = numpy.array(data2['energies'])
            diff = ((k1-k2)**2).sum()
            if diff > 0.1:
                print(diff)
                print(filename)
                sys.exit(1)
            
        if "_dos." in filename:
            k1 = numpy.array(data1['dos'])
            k2 = numpy.array(data2['dos'])
            diff = ((k1-k2)**2).sum()
            if diff > 0.1:
                print(diff)
                print(filename)
                sys.exit(1)
            
        if "path_" in filename and "_relaxation_times" in filename:
            k1 = numpy.array(data1['linewidths'])
            k2 = numpy.array(data2['linewidths'])
            diff = ((k1-k2)**2).sum()
            if diff > 0.1: 
                print(diff)
                print(filename)
                sys.exit(1)
            
    print("Reference checks Done")
    sys.exit(0)
