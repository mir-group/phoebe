#!/usr/bin/python
import numpy as np
import sys
import json

loc = "./"

#Get k and mode index from command line
if(len(sys.argv) != 6):
        print("""To run: python plotEigendisplacements.py <qx> <qy> <qz> <branch> <outputFile.xsf>
(Generates a file read by VESTA, which displays an eigendisplacement for the specified branch (indexed from zero) at wavevector q.""")
        exit(1)

# parse arguments
q = [float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])]
band = int(sys.argv[4])
outputFile = sys.argv[5]

# load in the json output
jfileName = 'phonon_bands.json' #args.INPUT
with open(jfileName) as jfile:
        data = json.load(jfile)

# read in and look up this qpoints
qlist = data["wavevectorCoordinates"]
qidx = qlist.index(q)
print("Found qpoint: ",qlist[qidx])

# read in eigendisplacements
# these should be in the shape (qpoints, bands, atoms, 3, 2) where 3 are x,y,z displacements
# and 2 is the real and complex parts
if "phononEigendisplacements" not in data:
  raise ValueError("JSON file does not contain phononEigendisplacements." +
        " You must first run phononBands app with outputEigendisplacements = true");

eigendisp = np.array(data["phononEigendisplacements"])[qidx,band]
# read in lattice
lattice = np.array(data["latticeVectors"])
# read in atom positions -- already output in cartesian
atomPos = np.array(data["atomPositions"])
numAtoms = atomPos.shape[0]

# atom species
atomSpecies = np.array(data["atomSpecies"])

Zdict = {
        'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16,
        'Cl':17, 'Ar':18, 'K':19, 'Ca':20, 'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32,
        'As':33, 'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y':39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48,
        'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I':53, 'Xe':54, 'Cs':55, 'Ba':56, 'La':57, 'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64,
        'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71, 'Hf':72, 'Ta':73, 'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80,
        'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 'Fr':87, 'Ra':88, 'Ac':89, 'Th':90, 'Pa':91, 'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96
}

# write to file
fp = open(outputFile, 'w')
fp.write('CRYSTAL\nPRIMVEC\n')
for i in range(3):
        fp.write('{} {} {}\n'.format(lattice[0,i], lattice[1,i], lattice[2,i]))
fp.write('PRIMCOORD\n')
fp.write(str(numAtoms) + ' 1\n')
# print the real part of the displacement to file as arrows
# Last index of eigendisp is the real/imag part
for iat in range(numAtoms):
        fp.write('{}  {} {} {}  {} {} {}\n'.format(
                Zdict[atomSpecies[iat]], atomPos[iat,0], atomPos[iat,1], atomPos[iat,2],
                eigendisp[iat,0,0], eigendisp[iat,1,0], eigendisp[iat,2,0]))

