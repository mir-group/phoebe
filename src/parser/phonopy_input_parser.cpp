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
#include "phonopy_input_parser.h"
#include "utilities.h"
#include "points.h"

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

int findRIndex(Eigen::MatrixXd &cellPositions2, Eigen::Vector3d &position2) {
  int ir2 = -1;
  for (int i = 0; i < cellPositions2.cols(); i++) {
    if ((position2 - cellPositions2.col(i)).norm() < 1.e-12) {
      ir2 = i;
      return ir2;
    }
  }
  if (ir2 == -1) {
    Error e("index not found");
  }
  return ir2;
}

std::tuple<Crystal, PhononH0> PhonopyParser::parsePhHarmonic(Context &context) {

  // Here we read the real space dynmat of interatomic force constants

  // This should be the directory containing all phonopy D2 files.
  // It should include fc2.hdf5, disp_fc2/3.yaml,and phono3py_disp.yaml
  std::string directory = context.getPhD2FileName();
  std::string line;
  if (directory == "") {
    Error e("Must provide a path to a directory containing fc2.hdf5, \n"
        "phono3py_disp.yaml, and either disp_fc3.yaml (if p3py force constants \n"
        "were generated with same dim for fc2 and fc3 supercells) \n"
        "or disp_fc2.yaml (if supercell dims were different).", 1);
  }

  // Read disp_fc2.yaml or disp_fc3.yaml file
  // ====================================================
  // which has the fc2 supercell regardless of whether or not forces
  // were generated with different fc2 vs fc3 supercells, or if both
  // supercells were the same, from disp_fc3.yaml.
  // If both used the same supercell, disp_fc2.yaml will not have been
  // created.

  // open input file
  std::ifstream infile(directory+"/disp_fc2.yaml");
  // if there's no disp_fc2 file, use disp_fc3 instead
  if(!infile.is_open()) {
    infile.clear();
    infile.open(directory+"/disp_fc3.yaml");
    if(!infile.is_open()) {
      Error e("disp_fc2.yaml/disp_fc3.yaml file not found at path " + directory, 1);
    }
  }

  // first line of disp.yaml will be natoms in supercell
  std::getline(infile, line);
  int numSupAtoms = std::stoi(line.substr(line.find(" ") ,line.back()));

  // read the rest of the file to get supercell positions
  Eigen::MatrixXd supPositions(numSupAtoms, 3);
  Eigen::MatrixXd supLattice(3, 3);
  int ilatt = 3;
  int ipos = 0;
  while(infile) {

    getline(infile, line);

    // if this is a cell position, save it
    if(line.find("position: ") != std::string::npos) {
      std::string temp = line.substr(14,57); // just the positions
      int idx1 = temp.find(",");
      supPositions(ipos,0) = std::stod(temp.substr(0,idx1));
      int idx2 = temp.find(",", idx1+1);
      supPositions(ipos,1) = std::stod(temp.substr(idx1+1,idx2));
      supPositions(ipos,2) = std::stod(temp.substr(idx2+1));
      ipos++;
    }
    if(ilatt < 3) { // count down lattice lines
      // convert from angstrom to bohr
      std::string temp = line.substr(5,62); // just the elements
      int idx1 = temp.find(",");
      supLattice(ilatt,0) = std::stod(temp.substr(0,idx1))/distanceBohrToAng;
      int idx2 = temp.find(",", idx1+1);
      supLattice(ilatt,1) = std::stod(temp.substr(idx1+1,idx2))/distanceBohrToAng;
      supLattice(ilatt,2) = std::stod(temp.substr(idx2+1))/distanceBohrToAng;
      ilatt++;
    }
    if(line.find("lattice:") != std::string::npos) {
      ilatt = 0;
    }
  }
  infile.close();

  // Read phono3py_disp.yaml file
  // ===================================================
  // Read in unit cell crystal information from phono3py_disp.yaml
  // This file is created at the start of a phono3py run, and contains
  // information specifying if fc2 has the same supercell dims as fc3

  // open input file
  infile.clear();
  infile.open(directory+"/phono3py_disp.yaml");

  // open input file
  if (not infile.is_open()) {
    Error e("phono3py_disp.yaml file not found in " + directory, 1);
  }

  // read in the dimension information.
  // we have to do this first, because we need to use this info
  // to allocate the below data storage.
  Eigen::Vector3i qCoarseGrid;
  while(infile) {
    getline(infile, line);

    // In the case where force constants where generated with different
    // supercells for fc2 and fc3, the label we need is dim_fc2.
    // Otherwise, if both are the same, it's just dim.
    if(line.find("dim_fc2: ") != std::string::npos) {
      std::string temp = line.substr(line.find("\"")+1, line.find("\"\n")-3);
      std::istringstream iss(temp);
      iss >> qCoarseGrid[0] >> qCoarseGrid[1] >> qCoarseGrid[2];
      break;
    }
    // this comes up first in the file, so if dim_fc2 is present, we'll overwrite it
    if(line.find("dim: ") != std::string::npos) {
      std::string temp = line.substr(10,5);
      std::istringstream iss(temp);
      iss >> qCoarseGrid(0) >> qCoarseGrid(1) >> qCoarseGrid(2);
    }
    // pause reading here
    if(line.find("physical_unit:") != std::string::npos) break;
  }

  if (qCoarseGrid(0) <= 0 || qCoarseGrid(1) <= 0 || qCoarseGrid(2) <= 0) {
    Error e("Phonon supercell dims read as 0 or less.\n"
        "Something is wrong with your input.", 1);
  }

  // set number of unit cell atoms
  int numAtoms = numSupAtoms/(qCoarseGrid(0)*qCoarseGrid(1)*qCoarseGrid(2));

  // read the rest of the file to find atomic positions,
  // lattice vectors, and species
  Eigen::Matrix3d directUnitCell;

  Eigen::MatrixXd atomicPositions(numAtoms, 3);
  Eigen::VectorXi atomicSpecies(numAtoms);

  std::vector<std::string> speciesNames;
  PeriodicTable pt;

  ilatt = 3;
  ipos = 0;
  // TODO watchout, this is reading the primitive cell from phono3py.
  // we might want to read in the unit cell, which could be different
  // because of some conversions they do internally.
  // Unit cell is also written to this fine in the same way as read
  // below.
  while(infile) {
    getline(infile, line);

    // if this line has a species, save it
    if(line.find("symbol: ") != std::string::npos) {
      std::string temp = line.substr(line.find("symbol: ")+8,line.find("#")-13);
      // remove any trailing whitespaces
      temp.erase(std::remove_if(temp.begin(), temp.end(), ::isspace), temp.end());
      if(std::find(speciesNames.begin(), speciesNames.end(), temp) == speciesNames.end()) {
        speciesNames.push_back(temp);
      }
      // save the atom number of this species
      atomicSpecies(ipos) = std::find(speciesNames.begin(), speciesNames.end(), temp) - speciesNames.begin();
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
      std::string temp = line.substr(9,67); // just the elements
      int idx1 = temp.find(",");
      directUnitCell(ilatt,0) = std::stod(temp.substr(0,idx1))/distanceBohrToAng;
      int idx2 = temp.find(",", idx1+1);
      directUnitCell(ilatt,1) = std::stod(temp.substr(idx1+1,idx2))/distanceBohrToAng;
      directUnitCell(ilatt,2) = std::stod(temp.substr(idx2+1,temp.find("]")))/distanceBohrToAng;
      ilatt++;
    }
    if(line.find("lattice:") != std::string::npos) {
      ilatt = 0;
    }
    // this signals we are done reading primitive cell info
    if(line.find("reciprocal_lattice:") != std::string::npos) {
      break;
    }
  }
  infile.close();

  // Process the information that has been read in
  // =================================================
  // calculate species mass values (in Ry) for later use
  Eigen::VectorXd speciesMasses(speciesNames.size());
  int count = 0;
  for ( auto i : speciesNames ) {
    speciesMasses[count] = pt.getMass(i) * massAmuToRy;
    count++;
  }

  // convert supercell positions to cartesian, in bohr
  for(int i = 0; i<numSupAtoms; i++) {
    Eigen::Vector3d temp(supPositions(i,0),supPositions(i,1),supPositions(i,2));
    Eigen::Vector3d temp2 = supLattice.transpose() * temp; // supLattice already converted to Bohr
    supPositions(i,0) = temp2(0);
    supPositions(i,1) = temp2(1);
    supPositions(i,2) = temp2(2);
  }

  // convert unit cell positions to cartesian, in bohr
  for(int i = 0; i<numAtoms; i++) {
    Eigen::Vector3d temp(atomicPositions(i,0),atomicPositions(i,1),atomicPositions(i,2));
    Eigen::Vector3d temp2 =  directUnitCell.transpose() * temp; // lattice already in Bohr
    atomicPositions(i,0) = temp2(0);
    atomicPositions(i,1) = temp2(1);
    atomicPositions(i,2) = temp2(2);
  }

  // Determine the list of possible R2, R3 vectors
  // the distances from the unit cell to a supercell
  // nCells here is the number of unit cell copies in the supercell
  int nCells = qCoarseGrid(0)*qCoarseGrid(1)*qCoarseGrid(2);
  Eigen::MatrixXd cellPositions2(3, nCells);
  cellPositions2.setZero();
  for (int icell = 0; icell < nCells; icell++) {
      // find the non-WS cell R2 vectors which are
      // position of atomPosSupercell - atomPosUnitCell = R
      cellPositions2(0,icell) = supPositions(icell,0) - supPositions(0,0);
      cellPositions2(1,icell) = supPositions(icell,1) - supPositions(0,1);
      cellPositions2(2,icell) = supPositions(icell,2) - supPositions(0,2);
  }

  // If hasDielectric, look for BORN file
  // ============================================================
  // the BORN file contains the dielectric matrix on the first line, 
  // and the BECs of unique atoms on the following lines
  bool hasDielectric = false;
  Eigen::Matrix3d dielectricMatrix;
  dielectricMatrix.setZero();
  Eigen::Tensor<double, 3> bornCharges(numAtoms, 3, 3);
  bornCharges.setZero();

  infile.clear();
  infile.open(directory+"/BORN");
  if(infile.is_open()) {
    hasDielectric = true;
    if(mpi->mpiHead()) std::cout << "Using BORN file found in D2 directory." << std::endl;
    getline(infile,line); 

    // NOTE need to extract the list of which atoms are listed in this file
    // from the comment on the first line of the file. Unfortuantely, there is
    // not a better way to do this.
    std::string temp = line.substr(line.find("atoms")+5);
    std::vector<int> becList;
    std::istringstream iss(temp);
    int i;
    while(iss>>i) { becList.push_back(i); }

    int atomType = -1;
    while(infile) {
      getline(infile, line);
      // make sure it's not a comment
      if(line.find("#") != std::string::npos) continue;
      // make sure it's not blank
      if(line.find_first_not_of(' ') == std::string::npos) continue;
  
      // the first non-comment line in the file is the dielectric matrix
      if(atomType == -1) {
        std::istringstream iss(line);
        iss >> dielectricMatrix(0,0) >> dielectricMatrix(0,1) >> dielectricMatrix(0,2) >>  
          dielectricMatrix(1,0) >> dielectricMatrix(1,1) >> dielectricMatrix(1,2) >> 
          dielectricMatrix(2,0) >> dielectricMatrix(2,1) >> dielectricMatrix(2,2);
        atomType++;
      }
      else {
        int numDupes;
        if(atomType+1 >= becList.size()){ numDupes = numAtoms - (becList[atomType] - 1); }
        else{ numDupes = becList[atomType+1] - becList[atomType]; }
        for(int i = 0; i < numDupes; i++) {
          std::istringstream iss(line);
     	  int iat = i + (becList[atomType] - 1);
	  iss >> bornCharges(iat,0,0) >> bornCharges(iat,0,1) >> bornCharges(iat,0,2) >>  
             bornCharges(iat,1,0) >> bornCharges(iat,1,1) >> bornCharges(iat,1,2) >> 
             bornCharges(iat,2,0) >> bornCharges(iat,2,1) >> bornCharges(iat,2,2);
        }
        atomType++;
      }
    }
  }

// print BECs
/*  for(int iat = 0; iat< numAtoms; iat++) {
	  std::cout << bornCharges(iat,0,0) << " " << bornCharges(iat,0,1) << " " << bornCharges(iat,0,2) << " " <<  
             bornCharges(iat,1,0) << " " << bornCharges(iat,1,1) << " " << bornCharges(iat,1,2) << " " <<
             bornCharges(iat,2,0) << " " << bornCharges(iat,2,1) << " " << bornCharges(iat,2,2) << std::endl;
  }
*/

  // Parse the fc2.hdf5 file and read in the dynamical matrix
  // ==========================================================
  #ifndef HDF5_AVAIL
  Error e("Phono3py HDF5 output cannot be read if Phoebe is not built with HDF5.");
  #else

  // set up buffer to read entire matrix
  // have to use this because the phono3py data is shaped as a
  // 4 dimensional array, and eigen tensor is not supported by highFive
  std::vector<std::vector<std::vector<std::vector<double>>>> ifc2;
  std::vector<int> cellMap;

  try {
    // Open the hdf5 file
    HighFive::File file(directory+"/fc2.hdf5", HighFive::File::ReadOnly);

    // Set up hdf5 datasets
    HighFive::DataSet difc2 = file.getDataSet("/force_constants");
    HighFive::DataSet dcellMap = file.getDataSet("/p2s_map");

    // read in the ifc3 data
    difc2.read(ifc2);
    dcellMap.read(cellMap);
  }
  catch(std::exception& error) {
    Error e("Issue reading fc2.hdf5 file. Make sure it exists in " + directory +
        "\n and is not open by some other persisting processes.");
  }

  Eigen::Tensor<double, 7> forceConstants(
      3, 3, qCoarseGrid[0], qCoarseGrid[1], qCoarseGrid[2], numAtoms, numAtoms);

  // phonopy force constants are in ev/ang^2, convert to atomic
  double conversion = pow(distanceBohrToAng, 2) / energyRyToEv;

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
        R = r1 * directUnitCell.row(0) +
            r2 * directUnitCell.row(1) +
            r3 * directUnitCell.row(2);

        // use the find cell function to determine
        // the index of this cell in the list of
        // R vectors (named cellPositions from above)
	int ir = findRIndex(cellPositions2, R);

        // loop over the first atoms. Because we consider R1=0,
        // these are only primitive unit cell atoms.
        for (int iat = 0; iat < numAtoms; iat++) {
          for (int jat = 0; jat < numAtoms; jat++) {

            // Need to convert jat to supercell index
            // Atoms in supercell are ordered so that there is a
            // unit cell atom followed by numUnitcell-in-supercell-#-of-atoms,
            // which are representations of the unit cell atom in other cells
            // we find this by using cellMap to tell us where the
            // first atom of this type is, then adding ir, which tells us
            // which cell it's in.
            int jsat = cellMap[jat] + ir;

            // loop over cartesian directions
            for (int ic : {0,1,2}) {
              for (int jc : {0,1,2}) {

                // here, cellMap tells us the position of this
                // unit cell atom in the supercell of phonopy
                forceConstants(ic, jc, r1, r2, r3, iat, jat) =
                      ifc2[cellMap[iat]][jsat][ic][jc] * conversion;

              }
            }
          }
        }
      }
    }
  }
  #endif

  if(mpi->mpiHead()) std::cout << "Successfully parsed harmonic phonopy files." << std::endl;

  // Now we do postprocessing
  int dimensionality = context.getDimensionality();
  Crystal crystal(context, directUnitCell, atomicPositions, atomicSpecies,
                  speciesNames, speciesMasses, dimensionality);

  PhononH0 dynamicalMatrix(crystal, dielectricMatrix, bornCharges,
                           forceConstants, context.getSumRuleD2());

  return {crystal, dynamicalMatrix};
};
