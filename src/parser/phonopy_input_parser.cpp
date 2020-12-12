#include <algorithm> // to use .remove_if
#include <fstream>
#include <iomanip> // to declare istringstream
#include <iostream>
#include <math.h>   // round()
#include <stdlib.h> // abs()
#include <string>
#include <vector>

#include "constants.h"
#include "eigen.h"
#include "exceptions.h"
#include "particle.h"
#include "periodic_table.h"
//#include "qe_input_parser.h"
#include "utilities.h"
#include "full_points.h"

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

std::tuple<Crystal, PhononH0> P3pyParser::parsePhHarmonic(Context &context) {
  //  Here we read the dynamical matrix of interatomic force constants
  //    in real space.

  // First, we read in crystal information from phono3py.yaml
  // TODO have to call this file something else
  std::string fileName = context.getPhD2FileName();
  if (fileName == "") {
    Error e("Must provide a phono3py.yaml file.", 1);
  }
  // open input file
  if (not infile.is_open()) {
    Error e("phono3py.yaml file not found", 1);
  }

  // first line will always be natoms in supercell
  std::getline(infile, line);
  int numAtoms = std::stoi(line.substr(line.find(" ") ,line.back()));

  int numElements = std::stoi(lineSplit[0]);

  // read the rest of the file to find atomic positions, 
  // lattice vectors, and species
  Eigen::Matrix3d directUnitCell;
  Eigen::Vector3i qCoarseGrid;

  Eigen::MatrixXd atomicPositions(numAtoms, 3);
  Eigen::VectorXi atomicSpecies(numAtoms);

  std::vector<std::string> speciesNames;
  Eigen::VectorXd speciesMasses(numElements);

  int ilatt = 3;
  int ipos = 0;
  // TODO watchout, this is reading the primitive cell from phono3py. 
  // we might want to read in the unit cell, which could be different
  // because of some conversions they do internally. 
  // Unit cell is also written to this fine in the same way as read
  // below.
  while(infile) {
    getline(infile, line);

    if(line.find("dim: ") != std::string::npos) {
      std::string temp = line.substr(line.find("\""), line.find("\"\n"));
      std::istringstream iss(line);
      iss >> qCoarseGrid[0] >> qCoarseGrid[1] >> qCoarseGrid[2];
    }
    // if this line has a species, save it 
    if(line.find("symbol: ") != std::string::npos) {
      speciesNames.push_back(line.substr(13,line.find("#")-1));
    }
    // if this line has a mass, save it 
    if(line.find("mass: ") != std::string::npos) {
      speciesMasses(ipos) = std::stod(line.substr(10)); // TODO convert to ry?
    }
    // if this is a cell position, save it
    if(line.find("coordinates: ") != std::string::npos) {
      std::string temp = line.substr(19,59); // just the positions
      int idx1 = temp.find(",");
      atomicPositions(ipos,0) = std::stod(temp.substr(0,idx1));
      int idx2 = temp.find(",", idx1+1);
      atomicPositions(ipos,1) = std::stod(temp.substr(idx1+1,idx2));
      atomicPositions(ipos,2) = std::stod(temp.substr(idx2+1));
      ipos++;
    }
    // parse lattice vectors
    if(ilatt < 3) { // count down lattice lines 
      std::string temp = line.substr(5,62); // just the elements
      int idx1 = temp.find(",");
      directUnitCell(ilatt,0) = std::stod(temp.substr(0,idx1));
      int idx2 = temp.find(",", idx1+1);
      directUnitCell(ilatt,1) = std::stod(temp.substr(idx1+1,idx2));
      directUnitCell(ilatt,2) = std::stod(temp.substr(idx2+1));
      ilatt++;
    }
    if(line.find("lattice:") != std::string::npos) {
      ilatt = 0;
    }
    // only read lines before this point
    if(line.find("unit_cell:") != std::string::npos) {
      break;
    }
  }
  infile.close();

  // in the file, they appear in crystal coordinates
  // we convert from crystal to cartesian coordinates
  // TODO can we do eigen mat*vector better than this?
  for(int i = 0; i<ipos; i++) {
    Eigen::Vector3d temp(
        atomicPositions(i,0),atomicPositions(i,1),atomicPositions(i,2));
    Eigen::Vector3d temp2 = lattice * temp;
    temp2 = temp2 / distanceBohrToAng; // lattice vectors are in angstrom
    atomicPositions(i,0) = temp2(0);
    atomicPositions(i,1) = temp2(1);
    atomicPositions(i,2) = temp2(2);
  }

  // build the atomicSpecies list 
  // this is a list of integers specifying which species
  // number each element is 
  std::vector<std::string> species;
  for(int i = 0; i<ipos; i++) {
    int speciesIdx = std::find(species.begin(), species.end(), speciesNames[ipos]);
    atomicSpecies(i) = speciesIdx;
    // species was not in the list
    if(speciesIdx == species.end()) {
      species.push_back(speciesNames(ipos));
    }
  }

  //  Read if hasDielectric
  hasDielectric = false; // TODO for now, we just say no dielectric
  Eigen::Matrix3d dielectricMatrix;
  dielectricMatrix.setZero();
  Eigen::Tensor<double, 3> bornCharges(numAtoms, 3, 3);
  bornCharges.setZero();

  // Parse the fc2.hdf5 file and read in the dynamical matrix 
  #ifndef HDF5_AVAIL
    Error e("Phono3py HDF5 output cannot be read if Phoebe is not built with HDF5.");
    //return void;
  #else

  // now we parse the coarse q grid
  fileName = context.getPhD2FileName();
  if (fileName == "") {
    Error e("Must provide a D2 file name, like fc2.hdf5", 1);
  }

  // Open the hdf5 file
  HighFive::File file(fileName, HighFive::File::ReadOnly);

  // Set up hdf5 datasets
  HighFive::DataSet difc2 = file.getDataSet("/fc2");

  // set up buffer to read entire matrix
  // have to use this because the phono3py data is shaped as a
  // 4 dimensional array, and eigen tensor is not supported by highFive
  std::vector<std::vector<std::vector<std::vector<double>>>> ifc2;

  // read in the ifc3 data
  difc2.read(ifc2);

  Eigen::Tensor<double, 7> forceConstants(
      3, 3, qCoarseGrid[0], qCoarseGrid[1], qCoarseGrid[2], numAtoms, numAtoms);

  // for the second atom, we must loop over all possible
  // cells in the supercell containing copies of these 
  // unit cell atoms
  for (int r3 = 0; r3 < qCoarseGrid[2]; r3++) {
    for (int r2 = 0; r2 < qCoarseGrid[1]; r2++) {
      for (int r1 = 0; r1 < qCoarseGrid[0]; r1++) {

        // NOTE we do this because phonopy has an 
        // "old" and "new" format for supercell files, 
        // and there's not an eay way for us to tell which 
        // one a user might have loaded in. Therefore, 
        // we can't intelligently guess the ordering of 
        // atoms in the supercell, and instead use R
        // to find their index. 

        // build the R vector associated with this 
        // cell in the supercell
        Eigen::Vector3d R;
        R = ir1 * directUnitCell(0) +
            ir2 * directUnitCell(1) +
            ir3 * directUnitCell(2);

        // use the find cell function to determine
        // the index of this cell in the list of 
        // R vectors (named cellPositions from above)
        //TODO copy this function in from ifc3
        int ir =

        // loop over the first atoms. Because we consider R1=0, 
        // these are only primitive unit cell atoms.
        for (int iat = 0; iat < numAtoms; iat++) {
          for (int jat = 0; jat < numAtoms; jat++) {
   
            // Need to convert jat to supercell index 
            // Atoms in supercell are ordered so that
            // there is a unit cell atom followed by 
            // numUnitcell-in-supercell-#-of-atoms
            // TODO lets think of an atom in the 8th cell. 
            // the 8th cell. the 8th cell atoms are every idx 
            // 7, 15, 23, 31.... that's ir + iat*numAtom.
            int jsat = ir + numAtoms * jat; 

            // loop over cartesian directions
            for (int ic : {0,1,2}) {
              for (int jc : {0,1,2}) {

                // here, cellMap tells us the position of this
                // unit cell atom in the supercell of phonopy
                forceConstants(ic, jc, r1, r2, r3, iat, jat) = 
                      ifc2[cellMap{iat}][jsat][ic][jc];

              }
            }
          }
        }
      }
    }
  }
  #endif

  // Now we do postprocessing
  long dimensionality = context.getDimensionality();
  Crystal crystal(directUnitCell, atomicPositions, atomicSpecies, speciesNames,
                  speciesMasses, dimensionality);

  if (qCoarseGrid(0) <= 0 || qCoarseGrid(1) <= 0 || qCoarseGrid(2) <= 0) {
    Error e("qCoarseGrid smaller than zero", 1);
  }

  PhononH0 dynamicalMatrix(crystal, dielectricMatrix, bornCharges,
                           forceConstants, context.getSumRuleD2());

  return {crystal, dynamicalMatrix};
};
