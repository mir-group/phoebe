#include "crystal.h"
#include "constants.h"
#include "exceptions.h"
#include "spglib.h"
#include "utilities.h"
#include "eigen.h"
#include "mpi/mpiHelper.h"

double calcVolume(const Eigen::Matrix3d &directUnitCell) {
  Eigen::Vector3d a1 = directUnitCell.row(0);
  Eigen::Vector3d a2 = directUnitCell.row(1);
  Eigen::Vector3d a3 = directUnitCell.row(2);
  double volume;
  volume = abs(a1.dot((a2.cross(a3))));
  volume += abs(a2.dot((a3.cross(a1))));
  volume += abs(a3.dot((a1.cross(a2))));
  volume /= 3.;
  return volume;
}

Crystal::Crystal(Context &context, Eigen::Matrix3d &directUnitCell_,
                 Eigen::MatrixXd &atomicPositions_,
                 Eigen::VectorXi &atomicSpecies_,
                 std::vector<std::string> &speciesNames_,
                 Eigen::VectorXd &speciesMasses_, int &dimensionality_) {

  setDirectUnitCell(directUnitCell_); // sets both direct and reciprocal
  volumeUnitCell = calcVolume(directUnitCell);

  if (volumeUnitCell <= 0.) {
    Error e("Unexpected non positive volume");
  }

  dimensionality = dimensionality_;

  if (atomicSpecies_.size() != atomicPositions_.rows()) {
    Error e("atomic species and positions are not aligned");
  }
  if (atomicPositions_.cols() != 3) {
    Error e("atomic positions need three coordinates");
  }
  if ((int)speciesMasses_.size() != (int)speciesNames_.size()) {
    Error e("species masses and names are not aligned");
  }

  atomicSpecies = atomicSpecies_;
  atomicPositions = atomicPositions_;
  speciesMasses = speciesMasses_;
  speciesNames = speciesNames_;

  numAtoms = atomicPositions.rows();
  numSpecies = speciesNames.size();

  Eigen::VectorXd atomicMasses_(numAtoms);
  std::vector<std::string> atomicNames_(numAtoms);

  for (int i = 0; i < numAtoms; i++) {
    atomicMasses_(i) = speciesMasses(atomicSpecies(i));
    atomicNames_[i] = speciesNames[atomicSpecies(i)];
  }
  atomicMasses = atomicMasses_;
  atomicNames = atomicNames_;

  int maxSize = 50;
  int rotations[maxSize][3][3];
  double translations[maxSize][3];

  if (context.getUseSymmetries()) {

    // We now look for the symmetry operations of the crystal
    // in this implementation, we rely on spglib

    // Declare and allocate c-style arrays for spglib calls
    double latticeSPG[3][3];
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        // note: directUnitCell has lattice vectors along rows (a1 = cell.row(0))
        latticeSPG[i][j] = directUnitCell(i, j);
      }
    }

    // note: spglib wants fractional positions
    double(*positionSPG)[3];
    allocate(positionSPG, numAtoms);
    Eigen::Vector3d positionCrystal;
    Eigen::Vector3d positionCartesian;
    for (int i = 0; i < numAtoms; i++) {
      positionCartesian = atomicPositions.row(i);
      positionCrystal = directUnitCell.inverse() * positionCartesian;
      for (int j = 0; j < 3; j++) {
        positionSPG[i][j] = positionCrystal(j);
      }
    }

    // also wants integer types >= 1
    int *typesSPG;
    allocate(typesSPG, numAtoms);
    for (int i = 0; i < numAtoms; i++) {
      typesSPG[i] = atomicSpecies(i) + 1;
    }

    double symprec = 1e-5;
    numSymmetries =
        spg_get_symmetry(rotations, translations, maxSize, latticeSPG,
                         positionSPG, typesSPG, numAtoms, symprec);

    // need to explicitly deallocate allocated arrays.
    delete[] typesSPG;
    delete[] positionSPG;

    if (numSymmetries == 0) {
      Error e("SPGlib failed at recognizing symmetries");
    }

    if (mpi->mpiHead()) {
      std::cout << "Found " << numSymmetries << " symmetries\n";
    }

  } else { // if we disable symmetries, and just use the identity

    if (mpi->mpiHead()) {
      std::cout << "Disabling symmetries\n";
    }
    numSymmetries = 1;
    for ( int i : {0,1,2}) {
      for ( int j : {0,1,2}) {
        rotations[0][i][j] = 0.;
      }
    }
    for ( int i : {0,1,2}) {
      translations[0][i] = 0.;
      rotations[0][i][i] = 1.;
    }
  }

  // If the system has no symmetries, we disable them
  if (numSymmetries==1) {
    context.setUseSymmetries(false);
  }

  // store the symmetries inside the class
  // note: spglib returns rotation and translation in fractional coordinates

  for (int isym = 0; isym < numSymmetries; isym++) {
    Eigen::Vector3d thisTranslation;
    thisTranslation(0) = translations[isym][0];
    thisTranslation(1) = translations[isym][1];
    thisTranslation(2) = translations[isym][2];
    Eigen::Matrix3d thisMatrix;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        thisMatrix(j, i) = rotations[isym][i][j]; // note the transpose
      }
    }

    SymmetryOperation s = {thisMatrix, thisTranslation};
    symmetryOperations.push_back(s);
  }
}

// empty constructor
Crystal::Crystal() {
  volumeUnitCell = 0.;
  numAtoms = 0;
  numSpecies = 0;
  dimensionality = 0.;
  numSymmetries = 0;
}

// copy constructor
Crystal::Crystal(const Crystal &obj) {
  directUnitCell = obj.directUnitCell;
  reciprocalUnitCell = obj.reciprocalUnitCell;
  volumeUnitCell = obj.volumeUnitCell;
  numAtoms = obj.numAtoms;
  numSpecies = obj.numSpecies;
  dimensionality = obj.dimensionality;
  atomicPositions = obj.atomicPositions;
  atomicSpecies = obj.atomicSpecies;
  atomicNames = obj.atomicNames;
  atomicMasses = obj.atomicMasses;
  speciesNames = obj.speciesNames;
  speciesMasses = obj.speciesMasses;
  symmetryOperations = obj.symmetryOperations;
  numSymmetries = obj.numSymmetries;
}

// assignment operator
Crystal &Crystal::operator=(const Crystal &obj) {
  if (this != &obj) {
    directUnitCell = obj.directUnitCell;
    reciprocalUnitCell = obj.reciprocalUnitCell;
    volumeUnitCell = obj.volumeUnitCell;
    numAtoms = obj.numAtoms;
    numSpecies = obj.numSpecies;
    dimensionality = obj.dimensionality;
    atomicPositions = obj.atomicPositions;
    atomicSpecies = obj.atomicSpecies;
    atomicNames = obj.atomicNames;
    atomicMasses = obj.atomicMasses;
    speciesNames = obj.speciesNames;
    speciesMasses = obj.speciesMasses;
    symmetryOperations = obj.symmetryOperations;
    numSymmetries = obj.numSymmetries;
  }
  return *this;
}

Eigen::Matrix3d
Crystal::calcReciprocalCell(const Eigen::Matrix3d &directUnitCell) {
  Eigen::Matrix3d reciprocalCell = twoPi * directUnitCell.inverse().transpose();
  return reciprocalCell;
}

void Crystal::setDirectUnitCell(Eigen::Matrix3d directUnitCell_) {
  directUnitCell = directUnitCell_;
  reciprocalUnitCell = calcReciprocalCell(directUnitCell);
}

const Eigen::Matrix3d &Crystal::getDirectUnitCell() { return directUnitCell; }

const Eigen::Matrix3d &Crystal::getReciprocalUnitCell() {
  // note: reciprocalUnitCell is  in units of twoPi
  // i.e. must be multiplied by twoPi
  return reciprocalUnitCell;
}

const int &Crystal::getNumAtoms() { return numAtoms; }

double Crystal::getVolumeUnitCell(int dimensionality_) {
  double volume;
  if (dimensionality_ == 3) {
    volume = volumeUnitCell;
  } else if (dimensionality_ == 2) {
    volume = abs(directUnitCell(0, 0) * directUnitCell(1, 1) -
                 directUnitCell(0, 1) * directUnitCell(1, 0));
  } else {
    volume = directUnitCell(2, 2);
  }
  return volume;
}

const Eigen::MatrixXd &Crystal::getAtomicPositions() { return atomicPositions; }

const Eigen::VectorXi &Crystal::getAtomicSpecies() { return atomicSpecies; }

const std::vector<std::string> &Crystal::getAtomicNames() {
  return atomicNames;
}

const Eigen::VectorXd &Crystal::getAtomicMasses() { return atomicMasses; }

const std::vector<std::string> &Crystal::getSpeciesNames() {
  return speciesNames;
}

const Eigen::VectorXd &Crystal::getSpeciesMasses() { return speciesMasses; }

const std::vector<SymmetryOperation> &Crystal::getSymmetryOperations() {
  return symmetryOperations;
}

const int &Crystal::getNumSymmetries() { return numSymmetries; }

int Crystal::getDimensionality() { return dimensionality; }

int Crystal::getNumSpecies() { return numSpecies; }

std::tuple<Eigen::MatrixXd, Eigen::VectorXd>
Crystal::buildWignerSeitzVectors(const Eigen::Vector3i &grid,
                                 const int &supercellFactor) {

  int nx = supercellFactor;

  std::vector<Eigen::Vector3d> tmpVectors;
  std::vector<double> tmpDegeneracies;

  // what are we doing:
  // the "n" specifies a grid of lattice vectors, that are (2*(cutoff+1)+1)^3
  // times bigger than the grid of wavevectors. a grid of lattice vectors
  // equal to the grid of wavevectors would be the bare minimum for the
  // interpolation to work, and wouldn't be enough for anything good.

  // now, we loop over the vectors of super cell A, build a vector r.
  for (int i0 = -nx * grid(0); i0 <= nx * grid(0); i0++) {
    for (int i1 = -nx * grid(1); i1 <= nx * grid(1); i1++) {
      for (int i2 = -nx * grid(2); i2 <= nx * grid(2); i2++) {
        // loop over a super cell B of size ((searchSize+1)*2+1)^3
        // bigger than super cell A. We compute the distance between any
        // vector in super cell B w.r.t. the vector "n"

        // We calculate |r-R|^2
        std::vector<double> distances;
        for (int j0 = -nx - 1; j0 <= nx + 1; j0++) {
          for (int j1 = -nx - 1; j1 <= nx + 1; j1++) {
            for (int j2 = -nx - 1; j2 <= nx + 1; j2++) {
              Eigen::Vector3d dist;
              dist(0) = i0 - j0 * grid(0);
              dist(1) = i1 - j1 * grid(1);
              dist(2) = i2 - j2 * grid(2);
              // distances in cartesian space
              dist = directUnitCell * dist;
              double dist2 = dist.norm();
              distances.push_back(dist2);
            }
          }
        }

        // find the minimum distance out of all distances
        double distMin = distances[0]; // [0] is just the first guess
        for (double dist : distances) {
          if (dist < distMin) {
            distMin = dist;
          }
        }

        // the point at midIndex is the reference "n" vector
        // i.e. the one generated by i0=i1=i2=0
        // if it's the minimum vector, than it's in the
        // Wigner Seitz zone of super cell A.
        // Therefore we save this vector.
        unsigned int midIndex = distances.size() / 2;
        if (abs(distances[midIndex] - distMin) < 1.0e-6) {
          // count its degeneracy
          double degeneracy = 0.;
          for (double dist : distances) {
            if (abs(dist - distMin) < 1.0e-6) {
              degeneracy += 1.;
            }
          }
          tmpDegeneracies.push_back(degeneracy);

          // store vector
          Eigen::Vector3d thisVec;
          thisVec << i0, i1, i2;
          tmpVectors.push_back(thisVec);
        }
      }
    }
  }

  // now we store the list of these lattice vectors in the class members
  int numPositionVectors = tmpVectors.size();
  Eigen::VectorXd positionDegeneracies(numPositionVectors);
  Eigen::MatrixXd positionVectors(3, numPositionVectors);
  int originIndex = -1; // to look for R=0 vector
  for (int iR = 0; iR < numPositionVectors; iR++) {
    auto thisVec = tmpVectors[iR];
    // we convert from crystal to cartesian coordinates
    positionVectors.col(iR) = directUnitCell * thisVec;
    positionDegeneracies(iR) = tmpDegeneracies[iR];
    //
    if (thisVec.norm() < 1.0e-6) {
      originIndex = iR;
    }
  }

  // As a convention, we keep the origin vectors at the index iR = 0
  // so we swap it. The other vectors aren't in any particular order.
  double tmp1 = positionDegeneracies(originIndex);
  double tmp2 = positionDegeneracies(0);
  positionDegeneracies(0) = tmp1;
  positionDegeneracies(originIndex) = tmp2;
  Eigen::Vector3d tmpv1 = positionVectors.col(originIndex);
  Eigen::Vector3d tmpv2 = positionVectors.col(0);
  positionVectors.col(0) = tmpv1;
  positionVectors.col(originIndex) = tmpv2;

  // check that we found all vectors
  {
    double tot = positionDegeneracies.sum();
    assert(abs(tot - grid(0) * grid(1) * grid(2)) < 1.0e-6);
    (void)tot;
  }

  return {positionVectors, positionDegeneracies};
}

std::tuple<Eigen::MatrixXd, Eigen::Tensor<double, 3>>
Crystal::buildWignerSeitzVectorsWithShift(const Eigen::Vector3i &grid,
                                          const Eigen::MatrixXd &shift,
                                          const int &supercellFactor) {

  // Note: we expect shift to be in cartesian coordinates!

  if (shift.rows() != 3) {
    Error e("shift should have at least 3 cartesian coordinates");
  }
  int shiftDims = shift.cols();

  int nx = supercellFactor;

  int maxVecs = 2 * nx * grid(0) + 1;
  maxVecs *= 2 * nx * grid(1) + 1;
  maxVecs *= 2 * nx * grid(2) + 1;
  Eigen::MatrixXd tmpNumPoints(shiftDims, shiftDims);
  Eigen::Tensor<double, 3> tmpDegeneraciesAll(shiftDims, shiftDims, maxVecs);
  Eigen::Tensor<double, 4> tmpVectorsAll(shiftDims, shiftDims, maxVecs, 3);

  // what are we doing:
  // the "n" specifies a grid of lattice vectors, that are (2*(cutoff+1)+1)^3
  // times bigger than the grid of wavevectors. a grid of lattice vectors
  // equal to the grid of wavevectors would be the bare minimum for the
  // interpolation to work, and wouldn't be enough for anything good.

  for (int idim = 0; idim < shiftDims; idim++) {
    for (int jdim = 0; jdim < shiftDims; jdim++) {
      // algorithm as before, except the shift

      std::vector<Eigen::Vector3d> tmpVectors;
      std::vector<double> tmpDegeneracies;

      // now, we loop over the vectors of super cell A, build a vector r.
      for (int i0 = -nx * grid(0); i0 <= nx * grid(0); i0++) {
        for (int i1 = -nx * grid(1); i1 <= nx * grid(1); i1++) {
          for (int i2 = -nx * grid(2); i2 <= nx * grid(2); i2++) {
            // loop over a super cell B of size ((searchSize+1)*2+1)^3
            // bigger than super cell A. We compute the distance between any
            // vector in super cell B w.r.t. the vector "n"

            // We calculate |r-R|^2
            std::vector<double> distances;
            for (int j0 = -nx - 1; j0 <= nx + 1; j0++) {
              for (int j1 = -nx - 1; j1 <= nx + 1; j1++) {
                for (int j2 = -nx - 1; j2 <= nx + 1; j2++) {
                  Eigen::Vector3d dist;
                  dist(0) = i0 - j0 * grid(0);
                  dist(1) = i1 - j1 * grid(1);
                  dist(2) = i2 - j2 * grid(2);
                  // distances in cartesian space
                  dist = directUnitCell * dist;

                  dist += shift.col(jdim) - shift.col(idim);

                  double dist2 = dist.norm();
                  distances.push_back(dist2);
                }
              }
            }

            // find the minimum distance out of all distances
            double distMin = distances[0]; // [0] is just the first guess
            for (double dist : distances) {
              if (dist < distMin) {
                distMin = dist;
              }
            }

            // the point at midIndex is the reference "n" vector
            // i.e. the one generated by i0=i1=i2=0
            // if it's the minimum vector, than it's in the
            // Wigner Seitz zone of super cell A.
            // Therefore we save this vector.
            unsigned int midIndex = distances.size() / 2;
            if (abs(distances[midIndex] - distMin) < 1.0e-6) {
              // count its degeneracy
              double degeneracy = 0.;
              for (double dist : distances) {
                if (abs(dist - distMin) < 1.0e-6) {
                  degeneracy += 1.;
                }
              }

              tmpDegeneracies.push_back(degeneracy);

              // store vector
              Eigen::Vector3d thisVec;
              thisVec << i0, i1, i2;
              tmpVectors.push_back(thisVec);
            }
          }
        }
      }

      tmpNumPoints(idim, jdim) = tmpVectors.size();
      for (unsigned int iR = 0; iR < tmpVectors.size(); iR++) {
        tmpDegeneraciesAll(idim, jdim, iR) = tmpDegeneracies[iR];
        for (int i : {0, 1, 2}) {
          tmpVectorsAll(idim, jdim, iR, i) = tmpVectors[iR](i);
        }
      }
    }
  }

  // find the unique set of vectors, so it's a smaller list
  std::vector<Eigen::Vector3d> tmpVectors;
  for (int idim = 0; idim < shiftDims; idim++) {
    for (int jdim = 0; jdim < shiftDims; jdim++) {
      for (int iR = 0; iR < tmpNumPoints(idim, jdim); iR++) {
        Eigen::Vector3d x;
        for (int i : {0, 1, 2}) {
          x(i) = tmpVectorsAll(idim, jdim, iR, i);
        }
        bool found = false;
        if (tmpVectors.empty())
          found = true; // if empty, add vector
        for (const Eigen::Vector3d &y : tmpVectors) {
          if ((x - y).squaredNorm() < 1.0e-6) {
            found = true;
            break;
          }
        }
        if (found) {
          tmpVectors.push_back(x);
        }
      }
    }
  }
  int numVectors = tmpVectors.size();

  // save the degeneracies in a smaller tensor
  Eigen::Tensor<double, 3> degeneracies(shiftDims, shiftDims, numVectors);
  degeneracies.setZero(); // important!
  for (int idim = 0; idim < shiftDims; idim++) {
    for (int jdim = 0; jdim < shiftDims; jdim++) {
      for (int iR = 0; iR < tmpNumPoints(idim, jdim); iR++) {
        Eigen::Vector3d x;
        for (int i : {0, 1, 2}) {
          x(i) = tmpVectorsAll(idim, jdim, iR, i);
        }
        double deg = tmpDegeneraciesAll(idim, jdim, iR);
        int iRNew = 0;
        for (const Eigen::Vector3d &y : tmpVectors) {
          if ((x - y).squaredNorm() < 1.0e-6) {
            degeneracies(idim, jdim, iRNew) = deg;
            break;
          }
          iRNew += 1;
        }
      }
    }
  }

  // convert to cartesian coordinates
  Eigen::MatrixXd bravaisVectors(3, numVectors);
  for (int iR = 0; iR < numVectors; iR++) {
    auto thisVec = tmpVectors[iR];
    // we convert from crystal to cartesian coordinates
    bravaisVectors.col(iR) = directUnitCell * thisVec;
  }

  // check that we found all vectors
  for (int idim = 0; idim < shiftDims; idim++) {
    for (int jdim = 0; jdim < shiftDims; jdim++) {
      double tot = 0.;

      for (int iR = 0; iR < numVectors; iR++) {
        if (degeneracies(idim, jdim, iR) > 0) {
          tot += degeneracies(idim, jdim, iR);
        }
      }
      assert(abs(tot - grid(0) * grid(1) * grid(2)) < 1.0e-6);
    }
  }

  return {bravaisVectors, degeneracies};
}
