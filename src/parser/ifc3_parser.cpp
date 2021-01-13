#include <iostream>
#include <fstream>
#include "ifc3_parser.h"
#include "eigen.h"
#include "constants.h"
#include "mpiHelper.h"

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

Interaction3Ph IFC3Parser::parse(Context &context, Crystal &crystal) {

    auto fileName = context.getPhD3FileName();

    // TODO would it be better to add one input variable called
    // phononFileType=[phonopy, shengbte, qe, etc] and then have
    // D3FileName and D2FileName changed to one phononInputDirectory?
    // TODO I left out input parsing for qe -- I think we discussed that
    // we weren't planning to support this?

    // check if this is a phono3py fc3.hdf5 file or a shengbte file
    if (fileName.find(".hdf5") != std::string::npos) {
      return parseFromPhono3py(context, crystal);
    }
    else {
      return parseFromShengBTE(context, crystal);
    }
}

Eigen::MatrixXd IFC3Parser::wsinit(Crystal &crystal, Eigen::Vector3i qCoarseGrid) {
  const int nx = 2;
  int index = 0;
  const int nrwsx = 200;
  Eigen::MatrixXd directUnitCell = crystal.getDirectUnitCell();
  //directUnitCell.transposeInPlace();

  Eigen::MatrixXd unitCell(3, 3);
  unitCell.col(0) = directUnitCell.col(0) * qCoarseGrid(0);
  unitCell.col(1) = directUnitCell.col(1) * qCoarseGrid(1);
  unitCell.col(2) = directUnitCell.col(2) * qCoarseGrid(2);

  Eigen::MatrixXd tmpResult(3, nrwsx);

  for (int ir = -nx; ir <= nx; ir++) {
    for (int jr = -nx; jr <= nx; jr++) {
      for (int kr = -nx; kr <= nx; kr++) {
        for (int i : {0,1,2}) {
          tmpResult(i, index) =
              unitCell(i, 0) * ir + unitCell(i, 1) * jr + unitCell(i, 2) * kr;
        }

        if (tmpResult.col(index).squaredNorm() > 1.0e-6) {
          index += 1;
        }
        if (index > nrwsx) {
          Error e("WSInit > nrwsx", 1);
        }
      }
    }
  }
  int nrws = index;

  Eigen::MatrixXd rws(3, nrws);
  for (int i = 0; i < nrws; i++) {
    rws.col(i) = tmpResult.col(i);
  }
  return rws;
}

double wsweight(const Eigen::VectorXd &r, const Eigen::MatrixXd &rws) {
  int nreq = 1;
  for (int ir = 0; ir < rws.cols(); ir++) {
    double rrt = r.dot(rws.col(ir));
    double ck = rrt - rws.col(ir).squaredNorm() / 2.;
    if (ck > 1.0e-6) {
      return 0.;
    }
    if (abs(ck) < 1.0e-6) {
      nreq += 1;
    }
  }
  double x = 1. / (double)nreq;
  return x;
}

int findIndexRow(Eigen::MatrixXd &cellPositions, Eigen::Vector3d &position) {
  int ir2 = -1;
  for (int i = 0; i < cellPositions.cols(); i++) {
    if ((position - cellPositions.col(i)).norm() < 1.e-12) {
      ir2 = i;
      return ir2;
    }
  }
  if (ir2 == -1) { Error("index not found"); }
  return ir2;
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::VectorXd> reorderDynamicalMatrix(Crystal &crystal, Eigen::Vector3i qCoarseGrid,
        Eigen::MatrixXd rws, Eigen::Tensor<double,5>& mat3R, Eigen::MatrixXd cellPositions,
        std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>> ifc3Tensor,
        std::vector<int> cellMap) {

  int numAtoms = crystal.getNumAtoms();
  Eigen::MatrixXd atomicPositions = crystal.getAtomicPositions();
  int nr1Big = qCoarseGrid(0) * 2;
  int nr2Big = qCoarseGrid(1) * 2;
  int nr3Big = qCoarseGrid(2) * 2;
  Eigen::MatrixXd directUnitCell = crystal.getDirectUnitCell();

  // Count the number of bravais vectors we need to loop over
  int numBravaisVectors = 0;
  // Record how many BV with non-zero weight
  // belong to each atom, so that we can index iR3
  // properly later in the code
  Eigen::VectorXi startBVForAtom(numAtoms);
  int numBVThisAtom = 0;

  for(int na = 0; na < numAtoms; na++) {
    startBVForAtom(na) = numBVThisAtom;
    for(int nb = 0; nb < numAtoms; nb++) {

      // loop over all possible indices for the first vector
      for (int nr1 = -nr1Big; nr1 < nr1Big; nr1++) {
        for (int nr2 = -nr2Big; nr2 < nr2Big; nr2++) {
          for (int nr3 = -nr3Big; nr3 < nr3Big; nr3++) {

            // calculate the R2 vector for atom nb
            Eigen::Vector3d r2;
            for (int i : {0, 1, 2}) {
              r2(i) = nr1 * directUnitCell(i, 0) + nr2 * directUnitCell(i, 1) +
                     nr3 * directUnitCell(i, 2);
              r2(i) = r2(i) - atomicPositions(na, i) + atomicPositions(nb, i);
            }
            double weightR2 = wsweight(r2, rws);
            // count this vector if it contributes
            if(weightR2 > 0) {
              numBravaisVectors += 1;
              numBVThisAtom += 1;
            }
          }
        }
      }
    }
  }

  // next, we reorder the dynamical matrix along the bravais lattice vectors
  // similarly to what was done in the harmonic class.
  Eigen::MatrixXd bravaisVectors2 = Eigen::MatrixXd::Zero(3, numBravaisVectors);
  Eigen::MatrixXd bravaisVectors3 = Eigen::MatrixXd::Zero(3, numBravaisVectors);
  Eigen::VectorXd weights2 = Eigen::VectorXd::Zero(numBravaisVectors);
  Eigen::VectorXd weights3 = Eigen::VectorXd::Zero(numBravaisVectors);

  mat3R.resize(3*numAtoms, 3*numAtoms, 3*numAtoms, numBravaisVectors, numBravaisVectors);
  mat3R.setZero();

  int iR2 = -1; int iR3 = -1;
  double conversion = pow(distanceBohrToAng, 3) / energyRyToEv;

  // TODO this is not efficient, redo it
    // TODO there is almost certainly a more favorable way to order these
    // loops or perform something like the wscache in harmonic/phonon_h0
    // TODO there's likely a way to easily avoid findIndexRow
    // TODO likely a way to avoid recalculating wsweight every time

  // loop over all possible indices for the first vector
  for (int nr3 = -nr3Big; nr3 < nr3Big; nr3++) {
    for (int nr2 = -nr2Big; nr2 < nr2Big; nr2++) {
      for (int nr1 = -nr1Big; nr1 < nr1Big; nr1++) {

        // loop over the atoms involved in the first vector
        for (int na = 0; na < numAtoms; na++) {
          for (int nb = 0; nb < numAtoms; nb++) {

            // calculate the R2 vector for atom nb
            Eigen::Vector3d r2;
            Eigen::Vector3d ra2;
            for (int i : {0, 1, 2}) {
              r2(i) = nr1 * directUnitCell(i, 0) + nr2 * directUnitCell(i, 1) +
                   nr3 * directUnitCell(i, 2);
              ra2(i) = r2(i) - atomicPositions(na, i) + atomicPositions(nb, i);
            }

            // skip vectors which do not contribute
            double weightR2 = wsweight(ra2, rws);
            if(weightR2 == 0) continue;
            iR2 +=1; // if selected increment the R2 index

            // the iR3 index should start based on how many R3 Bravais
      // vectors we covered in the earlier loops over na.
      // The vector startBVForAtom holds the starting
      // points in iR3 for each atom (minus 1 for zero based
      // indexing)
      iR3 = startBVForAtom(na) -1;

            // loop over atoms involved in R3 vector
            for (int nc = 0; nc < numAtoms; nc++) {

              // loop over all possible indices for the R3 vector
              for (int mr3 = -nr3Big; mr3 < nr3Big; mr3++) {
                for (int mr2 = -nr2Big; mr2 < nr2Big; mr2++) {
                  for (int mr1 = -nr1Big; mr1 < nr1Big; mr1++) {

                    // calculate the R3 vector for atom nc
                    Eigen::Vector3d r3, ra3;
                    for (int i : {0, 1, 2}) {
                      r3(i) = mr1 * directUnitCell(i, 0) + mr2 * directUnitCell(i, 1) +
                           mr3 * directUnitCell(i, 2);
                      ra3(i) = r3(i) - atomicPositions(na, i) + atomicPositions(nc, i);
                    }
                    // continue only if this vector matters
                    double weightR3 = wsweight(ra3, rws);
                    if(weightR3 == 0) continue;
                    iR3 +=1;

                    // save the weights and vectors for the next step
                    weights2(iR2) = weightR2;
                    weights3(iR3) = weightR3;
                    bravaisVectors2.col(iR2) = r2;
                    bravaisVectors3.col(iR3) = r3;
                    //std::cout << "chosen r3 weight " << iR3 << " " << weightR3 << " " << r3(0) << " " << r3(1) << " " << r3(2) << std::endl;
                    //std::cout << "chosen r2 weight "<< iR2 << " " << weightR2 << " " << r2(0) << " " << r2(1) << " " << r2(2) << std::endl;

                    // calculate the positive quadrant equivalents
                    int m1 = mod((nr1 + 1), qCoarseGrid(0));
                    if (m1 <= 0) { m1 += qCoarseGrid(0); };
                    int m2 = mod((nr2 + 1), qCoarseGrid(1));
                    if (m2 <= 0) { m2 += qCoarseGrid(1); };
                    int m3 = mod((nr3 + 1), qCoarseGrid(2));
                    if (m3 <= 0) { m3 += qCoarseGrid(2);};
                    m1 += -1; m2 += -1; m3 += -1;

                    int p1 = mod((mr1 + 1), qCoarseGrid(0));
                    if (p1 <= 0) { p1 += qCoarseGrid(0); };
                    int p2 = mod((mr2 + 1), qCoarseGrid(1));
                    if (p2 <= 0) { p2 += qCoarseGrid(1);};
                    int p3 = mod((mr3 + 1), qCoarseGrid(2));
                    if (p3 <= 0) { p3 += qCoarseGrid(2); };
                    p1 += -1; p2 += -1; p3 += -1;

                    // build the unit cell equivalent vectors to look up
                    Eigen::Vector3d rp2, rp3;
                    for (int i : {0, 1, 2}) {
                      rp2(i) = m1 * directUnitCell(i, 0) + m2 * directUnitCell(i, 1) +
                           m3 * directUnitCell(i, 2);
                      rp3(i) = p1 * directUnitCell(i, 0) + p2 * directUnitCell(i, 1) +
                           p3 * directUnitCell(i, 2);
                    }
                    // look up the index of the original atom
                    auto ir2 = findIndexRow(cellPositions, rp2);
                    auto ir3 = findIndexRow(cellPositions, rp3);

                    // index back to supercell positions
                    int sat2 = cellMap[nb] + ir2;
                    int sat3 = cellMap[nc] + ir3;

                    //std::cout << "pos r3 ir3 sat2 "<< iR3 << " "  <<  ir3 << " " << sat2 << " " << rp3(0) << " " << rp3(1) << " " << rp3(2) << std::endl;
                    //std::cout << "pos r2 ir2 sat3 "<< iR2 << " "  <<  ir2 << " " << sat3 << " " <<  rp2(0) << " " << rp2(1) << " " << rp2(2) << std::endl;

        // loop over the cartesian directions
                    for (int i : {0, 1, 2}) {
                      for (int j : {0, 1, 2}) {
                        for (int k : {0, 1, 2}) {

                          // compress the cartesian indices and unit cell atom indices
                          // for the first three indices of D3
                          auto ind1 = compress2Indices(na, i, numAtoms, 3);
                          auto ind2 = compress2Indices(nb, j, numAtoms, 3);
                          auto ind3 = compress2Indices(nc, k, numAtoms, 3);

                          //std::cout << "ind 1 2 3 iR2 iR3 " << ind1 << " " << ind2 << " " << ind3 << " " << iR2 << " " << iR3 << " ijk cellMap[na] na nb nc sat2 sat3 " << i << j << k << " " << cellMap[na] << " " << na << " " << nb << " " << nc << " " << sat2 << " " << sat3 << std::endl;

                          mat3R(ind1, ind2, ind3, iR2, iR3) +=
                                ifc3Tensor[cellMap[na]][sat2][sat3][i][j][k] * conversion;
                        } // close cartesian loops
                      }
                    }
                  } // r3 loops
                }
              }
            } // close nc loop
          } // close R2 loops
        }
      }
    } // close nb loop
  } // close na loop

  return {bravaisVectors2, weights2, bravaisVectors3, weights3};
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

  // Open disp_fc3 file, read supercell positions, nSupAtoms
  // TODO we need to supply a path rather than a filename,
  // since in this case there's two files...
  auto directory = context.getPhD2FileName();
  std::ifstream infile(directory + "/disp_fc3.yaml");
  std::string line;
  if (not infile.is_open()) {
      Error e("Phono3py disp_fc3.yaml file not "
  "found in directory "+directory+".");
  }

  // first line will always be natoms in supercell
  std::getline(infile, line);
  int numSupAtoms = std::stoi(line.substr(line.find(" ") ,line.back()));
  int nr = numSupAtoms/numAtoms;

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
      // convert from angstrom to bohr
      std::string temp = line.substr(5,62); // just the elements
      int idx1 = temp.find(",");
      lattice(ilatt,0) = std::stod(temp.substr(0,idx1))/distanceBohrToAng;
      int idx2 = temp.find(",", idx1+1);
      lattice(ilatt,1) = std::stod(temp.substr(idx1+1,idx2))/distanceBohrToAng;
      lattice(ilatt,2) = std::stod(temp.substr(idx2+1))/distanceBohrToAng;
      ilatt++;
    }
    if(line.find("lattice:") != std::string::npos) {
      ilatt = 0;
    }
  }
  infile.close();
  infile.clear();

  // convert positions to cartesian, in bohr
  for(int i = 0; i<ipos; i++) {
    Eigen::Vector3d temp(supPositions(i,0),supPositions(i,1),supPositions(i,2));
    Eigen::Vector3d temp2 = lattice.transpose() * temp;
    supPositions(i,0) = temp2(0);
    supPositions(i,1) = temp2(1);
    supPositions(i,2) = temp2(2);
  }

  // Open the hdf5 file containing the IFC3s
  auto fileName = context.getPhD3FileName();
  if(fileName == "") fileName = directory + "fc3.hdf5";
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

  // Read dimension of supercell from phono3py_disp file --------------------------------------
  infile.open(directory+"/phono3py_disp.yaml");

  if (not infile.is_open()) {
    Error e("phono3py_disp.yaml file not found in " + directory, 1);
  }

  // read in the dimension information.
  Eigen::Vector3i qCoarseGrid;
  while(infile) {
    getline(infile, line);
    if(line.find("dim: ") != std::string::npos) {
      std::string temp = line.substr(10,5);
      std::istringstream iss(temp);
      iss >> qCoarseGrid(0) >> qCoarseGrid(1) >> qCoarseGrid(2);
    }
  }

  // Determine the list of possible R2, R3 vectors
  // the distances from the origin unit cell to another
  // unit cell in the positive quadrant supercell
  Eigen::MatrixXd cellPositions(3, nr);
  cellPositions.setZero();

  // get all possible ws cell vectors
  Eigen::MatrixXd rws = wsinit(crystal,qCoarseGrid);

  for (int is = 0; is < nr; is++) {
      // find all possible vectors R2, R3, which are
      // position of atomPosSupercell - atomPosUnitCell = R
      cellPositions(0,is) = supPositions(is,0) - supPositions(0,0);
      cellPositions(1,is) = supPositions(is,1) - supPositions(0,1);
      cellPositions(2,is) = supPositions(is,2) - supPositions(0,2);
  }

  // reshape to D3 format
  Eigen::Tensor<double, 5> D3;

  auto tup = reorderDynamicalMatrix(crystal, qCoarseGrid, rws, D3, cellPositions, ifc3Tensor, cellMap);
  Eigen::MatrixXd bravaisVectors2 = std::get<0>(tup);
  Eigen::VectorXd weights2 = std::get<1>(tup);
  Eigen::MatrixXd bravaisVectors3 = std::get<2>(tup);
  Eigen::VectorXd weights3 = std::get<3>(tup);

  if(mpi->mpiHead()) std::cout << "Successfully parsed anharmonic phono3py files.\n" << std::endl;

  // Create interaction3Ph object
  Interaction3Ph interaction3Ph(crystal, D3, bravaisVectors2, bravaisVectors3, weights2, weights3);

  return interaction3Ph;

#endif
}

Interaction3Ph IFC3Parser::parseFromShengBTE(Context &context, Crystal &crystal) {

  auto fileName = context.getPhD3FileName();

  // Open IFC3 file
  std::ifstream infile(fileName);
  std::string line;

  if (not infile.is_open()) {
      Error e("D3 file not found");
  }

  // Number of triplets
  std::getline(infile, line);
  int numTriplets = std::stoi(line);

  // Allocate readables
  Eigen::Tensor<double, 4> ifc3Tensor(3, 3, 3, numTriplets);
  ifc3Tensor.setZero();
  Eigen::Tensor<double, 3> cellPositions(numTriplets, 2, 3);
  cellPositions.setZero();
  Eigen::Tensor<int, 2> displacedAtoms(numTriplets, 3);
  displacedAtoms.setZero();

  for (int i = 0; i < numTriplets; i++) {      // loop over all triplets

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
    int i0;
    while (iss3 >> i0) {
      displacedAtoms(i, j) = i0 - 1;
      j++;
    }

    // Read the 3x3x3 force constants tensor
    int i1, i2, i3;
    double d4;
    double conversion = pow(distanceBohrToAng, 3) / energyRyToEv;
    for (int a : { 0, 1, 2 }) {
      for (int b : { 0, 1, 2 }) {
        for (int c : { 0, 1, 2 }) {
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
  for (int it = 0; it < numTriplets; it++) {

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

  for (int it = 0; it < numTriplets; it++) { // sum over all triplets
    int ia1 = displacedAtoms(it, 0);
    int ia2 = displacedAtoms(it, 1);
    int ia3 = displacedAtoms(it, 2);

    Eigen::Vector3d position2, position3;
    for (int ic : {0, 1, 2}) {
      position2(ic) = cellPositions(it, 0, ic);
      position3(ic) = cellPositions(it, 1, ic);
    }

    int ir2 = findIndexRow(cellPositions2, position2);
    int ir3 = findIndexRow(cellPositions3, position3);

    for (int ic1 : {0, 1, 2}) {
      for (int ic2 : {0, 1, 2}) {
        for (int ic3 : {0, 1, 2}) {

          auto ind1 = compress2Indices(ia1, ic1, numAtoms, 3);
          auto ind2 = compress2Indices(ia2, ic2, numAtoms, 3);
          auto ind3 = compress2Indices(ia3, ic3, numAtoms, 3);

          D3(ind1, ind2, ind3, ir2, ir3) = ifc3Tensor(ic3, ic2, ic1, it);
        }
      }
    }
  }
  Eigen::VectorXd weights(cellPositions2.cols());
  weights.setZero();
  for (int i = 0; i< cellPositions2.cols(); i++) {
        weights(i) += 1;
  }

  Interaction3Ph interaction3Ph(crystal, D3, cellPositions2, cellPositions3, weights, weights);

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

  int numAtoms = crystal.getNumAtoms();
  int numSpecies = crystal.getSpeciesMasses().size();

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

  int numTriplets = 0;

  // in this first loop we count the number of triplets that have
  // non zero derivative. This allows us to decrease the size of the matrix
  for (int na1 = 0; na1 < numAtoms; na1++) {
    for (int na2 = 0; na2 < numAtoms; na2++) {
      for (int na3 = 0; na3 < numAtoms; na3++) {
        for (int j1 : { 0, 1, 2 }) {
          for (int j2 : { 0, 1, 2 }) {
            for (int j3 : { 0, 1, 2 }) {
              // suppress compiler warnings
              (void) j1;
              (void) j2;
              (void) j3;

              // read 6 integer which represent cartesian and
              // atomic basis indices
              std::getline(infile, line);

              // read the # of elements in this block
              std::getline(infile, line);
              std::istringstream iss2(line);
              int nR;
              iss2 >> nR;

              for (int i = 0; i < nR; i++) {
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
  Eigen::Tensor<int, 2> displacedAtoms(numTriplets, 3);
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

  int it = 0;

  for (int na1 = 0; na1 < numAtoms; na1++) {
    for (int na2 = 0; na2 < numAtoms; na2++) {
      for (int na3 = 0; na3 < numAtoms; na3++) {
        for (int j1 : { 0, 1, 2 }) {
          for (int j2 : { 0, 1, 2 }) {
            for (int j3 : { 0, 1, 2 }) {

              // read 6 integer which represent cartesian and
              // atomic basis indices
              std::getline(infile, line);

              // read the # of elements in this block
              std::getline(infile, line);
              std::istringstream iss2(line);
              int nR;
              iss2 >> nR;

              for (int i = 0; i < nR; i++) {

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
  for (int it = 0; it < numTriplets; it++) {

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

  for (int it = 0; it < numTriplets; it++) { // sum over all triplets
    int ia1 = displacedAtoms(it, 0);
    int ia2 = displacedAtoms(it, 1);
    int ia3 = displacedAtoms(it, 2);

    Eigen::Vector3d position2, position3;
    for (int ic : {0, 1, 2}) {
      position2(ic) = cellPositions(it, 0, ic);
      position3(ic) = cellPositions(it, 1, ic);
    }

    int ir2 = findIndexRow(cellPositions2, position2);
    int ir3 = findIndexRow(cellPositions3, position3);

    for (int ic1 : {0, 1, 2}) {
      for (int ic2 : {0, 1, 2}) {
        for (int ic3 : {0, 1, 2}) {

          auto ind1 = compress2Indices(ia1, ic1, numAtoms, 3);
          auto ind2 = compress2Indices(ia2, ic2, numAtoms, 3);
          auto ind3 = compress2Indices(ia3, ic3, numAtoms, 3);

          D3(ind1, ind2, ind3, ir2, ir3) = ifc3Tensor(ic3, ic2, ic1, it);
        }
      }
    }
  }
  Eigen::VectorXd weights(cellPositions2.cols());
  for (int i = 0; i< cellPositions2.cols(); i++) { weights(i) += 1; }
  Interaction3Ph interaction3Ph(crystal, D3, cellPositions2, cellPositions3, weights, weights);

  return interaction3Ph;
}
