#include <iostream>
#include <fstream>
#include "ifc3_parser.h"
#include "eigen.h"
#include "constants.h"

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

Interaction3Ph IFC3Parser::parse(Context &context, Crystal &crystal) {
    auto fileName = context.getPhD3FileName();

    // Open IFC3 file
    //std::ifstream infile(fileName);

    //if (not infile.is_open()) {
    //    Error e("D3 file not found", 1);
    //}

    // ShengBTE has a single integer on the first line
    // QE has more (9). We use this to differentiate the formats
    // TODO should just add an input variable for this
    //std::string line;
    //std::getline(infile, line);
    //std::istringstream iss(line);
    //std::string item;
    //int counter = 0;
    //while (iss >> item) {
    //    counter++;
   // }

    return parseFromPhono3py(context,crystal);
    //if (counter == 1) {
    //    return parseFromShengBTE(context, crystal);
    //} else {
    //    return parseFromQE(context, crystal);
    //}
}

Interaction3Ph IFC3Parser::parseFromPhono3py(Context &context, Crystal &crystal) {

#ifndef HDF5_AVAIL

  Error e("Phono3py HDF5 output cannot be read if Phoebe is not built with HDF5.");
  //return void;

#else

  // Notes about p3py's fc3.hdf5 file:
  // 3rd order fcs are listed as (num_atom, num_atom, num_atom, 3, 3, 3)
  // stored in eV/Angstrom3
  // Look here for additional detals
  // https://phonopy.github.io/phono3py/output-files.html#fc3-hdf5

  // the information we need outside the ifcs3 is in a file called
  // disp_fc3.yaml, which contains the atomic positions of the
  // supercell and the displacements of the atoms
  // The mapping between the supercell and unit cell atoms is 
  // set up so that the first nAtoms*dimSup of the supercell are unit cell 
  // atom #1, and the second nAtoms*dimSup are atom #2, and so on.

  // First, read in the information form disp_fc3.yaml
  int numAtoms = crystal.getNumAtoms();
  int numSupAtoms;
  int numBands = numAtoms * 3;

  // Open disp_fc3 file, read supercell positions, nSupAtoms
  // TODO we need to supply a path rather than a filename,
  // since in this case there's two files...
  std::ifstream infile("disp_fc3.yaml");
  std::string line;
  if (not infile.is_open()) {
      Error e("Phono3py disp_fc3.yaml file not found");
  }

  // first line will always be natoms in supercell
  std::getline(infile, line);
  numSupAtoms = std::stoi(line.substr(line.find(" ") ,line.back()));

  // read the rest of the file to look for supercell positions
  Eigen::MatrixXd supPositions(numSupAtoms, 3);
  Eigen::MatrixXd lattice(3, 3);
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
      std::string temp = line.substr(5,62); // just the elements
      int idx1 = temp.find(",");
      lattice(ilatt,0) = std::stod(temp.substr(0,idx1));
      int idx2 = temp.find(",", idx1+1);
      lattice(ilatt,1) = std::stod(temp.substr(idx1+1,idx2));
      lattice(ilatt,2) = std::stod(temp.substr(idx2+1));      
      ilatt++;
    }
    if(line.find("lattice:") != std::string::npos) {
      ilatt = 0;
    }
  }
  infile.close();

  // convert positions to cartesian, remember to make in bohr
  for(int i = 0; i<ipos; i++) {
    Eigen::Vector3d temp(supPositions(i,0),supPositions(i,1),supPositions(i,2));
    Eigen::Vector3d temp2 = lattice * temp;
    temp2 = temp2 * ang2Bohr;
    supPositions(i,0) = temp2(0);
    supPositions(i,1) = temp2(1);
    supPositions(i,2) = temp2(2);
  }


  // Open the hdf5 file containing the IFC3s
  auto fileName = context.getPhD3FileName();
  HighFive::File file(fileName, HighFive::File::ReadOnly);

  // Set up hdf5 datasets
  HighFive::DataSet difc3 = file.getDataSet("/fc3");
  HighFive::DataSet dcellMap = file.getDataSet("/p2s_map");

  // set up buffer to read entire matrix
  // have to use this monstrosity because the phono3py data is shaped as a
  // 6 dimensional array, and eigen tensor is not supported by highFive
  std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> ifc3Tensor;
  std::vector<int> cellMap;

  // read in the ifc3 data
  difc3.read(ifc3Tensor);
  dcellMap.read(cellMap); 

  // Determine the list of unique triplets
  std::vector<std::tuple<int,int>> triplets;
  std::vector<std::tuple<int,int>> supTriplets;
  int numTriplets = 0;
  for (int isa2 = 0; isa2 < numSupAtoms; isa2++) { 
    for (int isa3 = 0; isa3 < numSupAtoms; isa3++) {

      // Convert supercell atoms 2 and 3 to unit cell 2 and 3
      // these indices are the index of each atom in the
      // unit cell -- not the p3py supercell. 
      // To get that information, use cellMap[ia2]
      int ia2 = floor(isa2 / numAtoms);
      int ia3 = floor(isa3 / numAtoms);
      std::tuple<int,int> triplet = {ia2,ia3}; 

      // check to see if this one is already in the list
      if(std::find(triplets.begin(), triplets.end(), triplet) != triplets.end()) {
        triplets.push_back(triplet);
        supTriplets.push_back({isa2,isa3});
        numTriplets++;
      }
    }
  }
  // num triplets goes over unit cell atoms, plus all the r2, r3 atom
  // combinations found above
  numTriplets *= numAtoms;

  // Allocate final storage of read in quantities
  // TODO this needs to be last two dims changed to #2nd atom moves, #3rd atom moves.
  Eigen::Tensor<double, 5> D3(numBands, numBands, numBands, numTriplets, numTriplets);
  D3.setZero();
  Eigen::MatrixXd cellPositions2(3, numTriplets);
  Eigen::MatrixXd cellPositions3(3, numTriplets);
  cellPositions2.setZero();
  cellPositions3.setZero();
  Eigen::MatrixXd unitCellPos = crystal.getAtomicPositions(); // TODO are these direct?

  // reshape to D3 format

  // For the first atom index, we only go over atoms in the unit cell, 
  // as we consider R1 = 0. 
  double conversion = pow(distanceBohrToAng, 3) / energyRyToEv;
  for (int ia1 = 0; ia1 < numAtoms; ia1++) {
    // loop over unique R2 and R3 atom triplets
    for ( int idxTriplet = 0; idxTriplet<numTriplets; idxTriplet++) {

      int ia2 = std::get<0>(triplets[idxTriplet]);
      int ia3 = std::get<1>(triplets[idxTriplet]); 
      int isa2 = std::get<0>(supTriplets[idxTriplet]);
      int isa3 = std::get<1>(supTriplets[idxTriplet]);

      // find the vectors R2, R3 for this triplet
      // position of atomPosSupercell - atomPosUnitCell = R
      cellPositions2.col(idxTriplet) = supPositions.row(isa2) - unitCellPos.row(ia2);
      cellPositions3.col(idxTriplet) = supPositions.row(isa3) - unitCellPos.row(ia3);
        
      for (int ic1 : { 0, 1, 2 }) {
        for (int ic2 : { 0, 1, 2 }) {
          for (int ic3 : { 0, 1, 2 }) {

            // mux the cartesian indices and unit cell atom indices 
            // for the first three indices of D3
            auto ind1 = compress2Indeces(ia1, ic1, numAtoms, 3);
            auto ind2 = compress2Indeces(ia2, ic2, numAtoms, 3);
            auto ind3 = compress2Indeces(ia3, ic3, numAtoms, 3);              

            // indices here are atom+cart index, then indices for displaced atoms 
            // in this triplet -- TODO, if there are issues, 
            // check the last two indices. 
            D3(ind1, ind2, ind3, idxTriplet, idxTriplet) // convert from ev/ang to atomic 
              = ifc3Tensor[cellMap[ia1]][isa2][isa3][ic1][ic2][ic3] * conversion;
          }
        }
      }
    }
  }

  // Create interaction3Ph object
  Interaction3Ph interaction3Ph(crystal, D3, cellPositions2, cellPositions3);

  return interaction3Ph;

#endif
}

long findIndexRow(Eigen::MatrixXd &cellPositions2, Eigen::Vector3d &position2) {
  long ir2 = -1;
  for (int i = 0; i < cellPositions2.cols(); i++) {
    if ((position2 - cellPositions2.col(i)).norm() == 0.) {
      ir2 = i;
      return ir2;
    }
  }
  if (ir2 == -1) {
    Error e("index not found");
  }
  return ir2;
}

Interaction3Ph IFC3Parser::parseFromShengBTE(Context &context,
        Crystal &crystal) {

  auto fileName = context.getPhD3FileName();

  // Open IFC3 file
  std::ifstream infile(fileName);
  std::string line;

  if (not infile.is_open()) {
      Error e("D3 file not found");
  }

  // Number of triplets
  std::getline(infile, line);
  long numTriplets = std::stoi(line);

  // Allocate readables
  Eigen::Tensor<double, 4> ifc3Tensor(3, 3, 3, numTriplets);
  ifc3Tensor.setZero();
  Eigen::Tensor<double, 3> cellPositions(numTriplets, 2, 3);
  cellPositions.setZero();
  Eigen::Tensor<long, 2> displacedAtoms(numTriplets, 3);
  displacedAtoms.setZero();

  for (long i = 0; i < numTriplets; i++) {	// loop over all triplets

    // empty line
    std::getline(infile, line);
    // line with a counter
    std::getline(infile, line);

    // Read position of 2nd cell
    std::getline(infile, line);
    std::istringstream iss(line);
    std::string item;
    int j = 0;
    while (iss >> item) {
      cellPositions(i, 0, j) = std::stod(item) / distanceBohrToAng;
      j++;
    }

    // Read position of 3rd cell
    std::getline(infile, line);
    std::istringstream iss2(line);
    j = 0;
    while (iss2 >> item) {
      cellPositions(i, 1, j) = std::stod(item) / distanceBohrToAng;
      j++;
    }

    // Read triplet atom indices
    std::getline(infile, line);
    std::istringstream iss3(line);
    j = 0;
    long i0;
    while (iss3 >> i0) {
      displacedAtoms(i, j) = i0 - 1;
      j++;
    }

    // Read the 3x3x3 force constants tensor
    long i1, i2, i3;
    double d4;
    double conversion = pow(distanceBohrToAng, 3) / energyRyToEv;
    for (long a : { 0, 1, 2 }) {
      for (long b : { 0, 1, 2 }) {
        for (long c : { 0, 1, 2 }) {
          std::getline(infile, line);
          std::istringstream iss4(line);
          while (iss4 >> i1 >> i2 >> i3 >> d4) {
            ifc3Tensor(c, b, a, i) = d4 * conversion;
          }
        }
      }
    }
  }

  // Close IFC3 file
  infile.close();

  //TODO Round cellPositions to the nearest lattice vectors

  // start processing ifc3s into D3 matrix
  int numAtoms = crystal.getNumAtoms();
  int numBands = numAtoms * 3;
  int nr2 = 0; 
  int nr3 = 0;
  std::vector<Eigen::Vector3d> tmpCellPositions2, tmpCellPositions3;

  // TODO can't we simplify this? why are we setting
  // found2 and found3?
  for (long it = 0; it < numTriplets; it++) {

    // load the position of the 2 atom in the current triplet
    Eigen::Vector3d position2, position3;
    for (int ic : {0, 1, 2}) {
      position2(ic) = cellPositions(it, 0, ic);
      position3(ic) = cellPositions(it, 1, ic);
    }

    // now check if this element is in the list.
    bool found2 = false;
    if (std::find(tmpCellPositions2.begin(), tmpCellPositions2.end(),
                  position2) != tmpCellPositions2.end()) {
      found2 = true;
    }
    bool found3 = false;
    if (std::find(tmpCellPositions3.begin(), tmpCellPositions3.end(),
                  position3) != tmpCellPositions3.end()) {
      found3 = true;
    }

    if (!found2) {
      tmpCellPositions2.push_back(position2);
      nr2++;
    }
    if (!found3) {
      tmpCellPositions3.push_back(position3);
      nr3++;
    }
  }

  Eigen::MatrixXd cellPositions2(3, nr2);
  Eigen::MatrixXd cellPositions3(3, nr3);
  cellPositions2.setZero();
  cellPositions3.setZero();
  for (int i = 0; i < nr2; i++) {
    cellPositions2.col(i) = tmpCellPositions2[i];
  }
  for (int i = 0; i < nr3; i++) {
    cellPositions3.col(i) = tmpCellPositions3[i];
  }

  Eigen::Tensor<double, 5> D3(numBands, numBands, numBands, nr2, nr3);
  D3.setZero();

  for (long it = 0; it < numTriplets; it++) { // sum over all triplets
    long ia1 = displacedAtoms(it, 0);
    long ia2 = displacedAtoms(it, 1);
    long ia3 = displacedAtoms(it, 2);

    Eigen::Vector3d position2, position3;
    for (int ic : {0, 1, 2}) {
      position2(ic) = cellPositions(it, 0, ic);
      position3(ic) = cellPositions(it, 1, ic);
    }

    long ir2 = findIndexRow(cellPositions2, position2);
    long ir3 = findIndexRow(cellPositions3, position3);

    for (int ic1 : {0, 1, 2}) {
      for (int ic2 : {0, 1, 2}) {
        for (int ic3 : {0, 1, 2}) {

          auto ind1 = compress2Indeces(ia1, ic1, numAtoms, 3);
          auto ind2 = compress2Indeces(ia2, ic2, numAtoms, 3);
          auto ind3 = compress2Indeces(ia3, ic3, numAtoms, 3);

          D3(ind1, ind2, ind3, ir2, ir3) = ifc3Tensor(ic3, ic2, ic1, it);
        }
      }
    }
  }

  Interaction3Ph interaction3Ph(crystal, D3, cellPositions2, cellPositions3);

  return interaction3Ph;
}

Interaction3Ph IFC3Parser::parseFromQE(Context &context, Crystal &crystal) {

  auto fileName = context.getPhD3FileName();

  const double threshold = 1.0e-20;

  // Open IFC3 file
  std::ifstream infile(fileName);
  std::string line;

  if (not infile.is_open()) {
    Error e("D3 file not found");
  }

  long numAtoms = crystal.getNumAtoms();
  long numSpecies = crystal.getSpeciesMasses().size();

  // The first few lines contain info on the crystal, which we ignore
  std::getline(infile, line);
  for (int i = 0; i < numSpecies; i++) {
    std::getline(infile, line);
  }
  for (int i = 0; i < numAtoms; i++) {
    std::getline(infile, line);
  }
  // line with dielectric
  std::getline(infile, line);

  // read the mesh of triplets
  std::getline(infile, line);

  long numTriplets = 0;

  // in this first loop we count the number of triplets that have
  // non zero derivative. This allows us to decrease the size of the matrix
  for (long na1 = 0; na1 < numAtoms; na1++) {
    for (long na2 = 0; na2 < numAtoms; na2++) {
      for (long na3 = 0; na3 < numAtoms; na3++) {
        for (long j1 : { 0, 1, 2 }) {
          for (long j2 : { 0, 1, 2 }) {
            for (long j3 : { 0, 1, 2 }) {
              // suppress compiler warnings
              (void) j1;
              (void) j2;
              (void) j3;

              // read 6 integer which represent cartesian and
              // atomic basis indeces
              std::getline(infile, line);

              // read the # of elements in this block
              std::getline(infile, line);
              std::istringstream iss2(line);
              long nR;
              iss2 >> nR;

              for (long i = 0; i < nR; i++) {
                std::getline(infile, line);
                std::istringstream iss3(line);
                int i1, i2, i3, i4, i5, i6;
                double x1;
                iss3 >> i1 >> i2 >> i3 >> i4 >> i5 >> i6 >> x1;

                if (abs(x1) > threshold) {
                  numTriplets++;
                }
              }
            }
          }
        }
      }
    }
  }
  infile.close();
  infile.clear();

  // now that we know the number of triplets, we build the matrix

  // Allocate readables
  Eigen::Tensor<double, 4> ifc3Tensor(3, 3, 3, numTriplets);
  ifc3Tensor.setZero();
  Eigen::Tensor<double, 3> cellPositions(numTriplets, 2, 3);
  cellPositions.setZero();
  Eigen::Tensor<long, 2> displacedAtoms(numTriplets, 3);
  displacedAtoms.setZero();

  // Open IFC3 file
  infile.open(fileName);

  // The first four lines contain info on the crystal, which we ignore
  std::getline(infile, line);
  for (int i = 0; i < numSpecies; i++) {
    std::getline(infile, line);
  }
  for (int i = 0; i < numAtoms; i++) {
    std::getline(infile, line);
  }
  std::getline(infile, line);

  // read the mesh of triplets
  std::getline(infile, line);

  long it = 0;

  for (long na1 = 0; na1 < numAtoms; na1++) {
    for (long na2 = 0; na2 < numAtoms; na2++) {
      for (long na3 = 0; na3 < numAtoms; na3++) {
        for (long j1 : { 0, 1, 2 }) {
          for (long j2 : { 0, 1, 2 }) {
            for (long j3 : { 0, 1, 2 }) {

              // read 6 integer which represent cartesian and
              // atomic basis indeces
              std::getline(infile, line);

              // read the # of elements in this block
              std::getline(infile, line);
              std::istringstream iss2(line);
              long nR;
              iss2 >> nR;

              for (long i = 0; i < nR; i++) {

                std::getline(infile, line);
                std::istringstream iss3(line);
                int i1, i2, i3, i4, i5, i6;
                double x1;
                iss3 >> i1 >> i2 >> i3 >> i4 >> i5 >> i6 >> x1;

                if (abs(x1) > threshold) {
                  cellPositions(it, 0, 0) = double(i1);
                  cellPositions(it, 0, 1) = double(i2);
                  cellPositions(it, 0, 2) = double(i3);

                  cellPositions(it, 1, 0) = double(i4);
                  cellPositions(it, 1, 1) = double(i5);
                  cellPositions(it, 1, 2) = double(i6);

                  // this is an index over the atomic basis
                  // for the three atoms in the triplet
                  displacedAtoms(it, 0) = na1;
                  displacedAtoms(it, 1) = na2;
                  displacedAtoms(it, 2) = na3;

                  ifc3Tensor(j3, j2, j1, it) = x1;
                  it++;
                }
              }
            }
          }
        }
      }
    }
  }
  infile.close();

  // start processing ifc3s into D3 matrix
  int numBands = numAtoms * 3;
  int nr2 = 0;
  int nr3 = 0;
  std::vector<Eigen::Vector3d> tmpCellPositions2, tmpCellPositions3;

  // TODO can't we simplify this? why are we setting
  // found2 and found3?
  for (long it = 0; it < numTriplets; it++) {

    // load the position of the 2 atom in the current triplet
    Eigen::Vector3d position2, position3;
    for (int ic : {0, 1, 2}) {
      position2(ic) = cellPositions(it, 0, ic);
      position3(ic) = cellPositions(it, 1, ic);
    }

    // now check if this element is in the list.
    bool found2 = false;
    if (std::find(tmpCellPositions2.begin(), tmpCellPositions2.end(),
                  position2) != tmpCellPositions2.end()) {
      found2 = true;
    }
    bool found3 = false;
    if (std::find(tmpCellPositions3.begin(), tmpCellPositions3.end(),
                  position3) != tmpCellPositions3.end()) {
      found3 = true;
    }

    if (!found2) {
      tmpCellPositions2.push_back(position2);
      nr2++;
    }
    if (!found3) {
      tmpCellPositions3.push_back(position3);
      nr3++;
    }
  }

  Eigen::MatrixXd cellPositions2(3, nr2);
  Eigen::MatrixXd cellPositions3(3, nr3);
  for (int i = 0; i < nr2; i++) {
    cellPositions2.col(i) = tmpCellPositions2[i];
  }
  for (int i = 0; i < nr3; i++) {
    cellPositions3.col(i) = tmpCellPositions3[i];
  }

  Eigen::Tensor<double, 5> D3(numBands, numBands, numBands, nr2, nr3);
  D3.setZero();

  for (long it = 0; it < numTriplets; it++) { // sum over all triplets
    long ia1 = displacedAtoms(it, 0);
    long ia2 = displacedAtoms(it, 1);
    long ia3 = displacedAtoms(it, 2);

    Eigen::Vector3d position2, position3;
    for (int ic : {0, 1, 2}) {
      position2(ic) = cellPositions(it, 0, ic);
      position3(ic) = cellPositions(it, 1, ic);
    }

    long ir2 = findIndexRow(cellPositions2, position2);
    long ir3 = findIndexRow(cellPositions3, position3);

    for (int ic1 : {0, 1, 2}) {
      for (int ic2 : {0, 1, 2}) {
        for (int ic3 : {0, 1, 2}) {

          auto ind1 = compress2Indeces(ia1, ic1, numAtoms, 3);
          auto ind2 = compress2Indeces(ia2, ic2, numAtoms, 3);
          auto ind3 = compress2Indeces(ia3, ic3, numAtoms, 3);

          D3(ind1, ind2, ind3, ir2, ir3) = ifc3Tensor(ic3, ic2, ic1, it);
        }
      }
    }
  }

  Interaction3Ph interaction3Ph(crystal, D3, cellPositions2, cellPositions3);

  return interaction3Ph;
}
