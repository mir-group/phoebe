#include "delta_function.h"
#include "constants.h"
#include "context.h"
#include "exceptions.h"
#include "utilities.h"

DeltaFunction::~DeltaFunction() {}

int DeltaFunction::getType() { return id; }

int GaussianDeltaFunction::getType() { return id; }

int AdaptiveGaussianDeltaFunction::getType() { return id; }

int TetrahedronDeltaFunction::getType() { return id; }

// app factory
DeltaFunction *
DeltaFunction::smearingFactory(Context &context,
                               BaseBandStructure &fullBandStructure) {
  auto choice = context.getSmearingMethod();
  if (choice == gaussian) {
    return new GaussianDeltaFunction(context);
  } else if (choice == adaptiveGaussian) {
    return new AdaptiveGaussianDeltaFunction(fullBandStructure);
  } else if (choice == tetrahedron) {
    return new TetrahedronDeltaFunction(fullBandStructure);
  } else {
    Error e("Unrecognized smearing choice");
    return nullptr;
  }
}

GaussianDeltaFunction::GaussianDeltaFunction(Context &context) {
  inverseWidth = 1. / context.getSmearingWidth();
  prefactor = 1. / context.getSmearingWidth() / sqrt(pi);
}

double GaussianDeltaFunction::getSmearing(const double &energy,
                                          const Eigen::Vector3d &velocity) {
  (void)velocity;
  double x = energy * inverseWidth;
  return prefactor * exp(-x * x);
}

double GaussianDeltaFunction::getSmearing(const double &energy,
                                          StateIndex &is) {
  (void)energy;
  (void)is;
  Error e("GaussianDeltaFunction::getSmearing2 not implemented");
  return 1.;
}

AdaptiveGaussianDeltaFunction::AdaptiveGaussianDeltaFunction(
    BaseBandStructure &bandStructure) {
  auto tup = bandStructure.getPoints().getMesh();
  auto mesh = std::get<0>(tup);
  qTensor = bandStructure.getPoints().getCrystal().getReciprocalUnitCell();
  qTensor.row(0) /= mesh(0);
  qTensor.row(1) /= mesh(1);
  qTensor.row(2) /= mesh(2);
}

double
AdaptiveGaussianDeltaFunction::getSmearing(const double &energy,
                                           const Eigen::Vector3d &velocity) {
  if (velocity.norm() == 0. && energy == 0.) {
    // in this case, velocities are parallel, there shouldn't be
    // scattering unless energy is strictly conserved
    return 1.;
  }

  double sigma = 0.;
  for (int i : {0, 1, 2}) {
    sigma += pow(qTensor.row(i).dot(velocity), 2);
  }
  sigma = prefactor * sqrt(sigma / 6.);

  if (sigma == 0.)
    return 0.;

  if (abs(energy) > 2. * sigma)
    return 0.;
  double x = energy / sigma;
  // note: the factor ERF_2 corrects for the cutoff at 2*sigma
  return exp(-x * x) / sqrtPi / sigma / erf2;
}

double AdaptiveGaussianDeltaFunction::getSmearing(const double &energy,
                                                  StateIndex &is) {
  (void)energy;
  (void)is;
  Error e("AdaptiveGaussianDeltaFunction::getSmearing2 not implemented");
  return 1.;
}

TetrahedronDeltaFunction::TetrahedronDeltaFunction(
    BaseBandStructure &fullBandStructure_)
    : fullBandStructure(fullBandStructure_) {
  auto fullPoints = fullBandStructure.getPoints();
  auto tup = fullPoints.getMesh();
  auto grid = std::get<0>(tup);
  auto offset = std::get<1>(tup);
  if (offset.norm() > 0.) {
    Error e("We didnt' implement tetrahedra with offsets", 1);
  }
  if (grid(0) < 1 || grid(1) < 1 || grid(2) < 1) {
    Error e("Tetrahedron method initialized with invalid wavevector grid");
  }

  // shifts of 1 kpoint in the grid in crystal coordinates
  Eigen::Vector3d deltaGrid;
  for (int i : {0,1,2}) {
    deltaGrid(i) = 1. / grid(i);
  }

  // in this tensor, we store the offset to find the vertices of the
  // tetrahedrons to which the current point belongs to
  // (6: number of tetrahedra, 4: number of vertices, 3: cartesian)
  // multiplying this by the vector Delta k of the grid, we can reconstruct
  // all the points of the tetrahedron.

  subcellShift.resize(8,3);
  subcellShift.row(0) << 0., 0., 0.;
  subcellShift.row(1) << 1., 0., 0.;
  subcellShift.row(2) << 0., 1., 0.;
  subcellShift.row(3) << 1., 1., 0.;
  subcellShift.row(4) << 0., 0., 1.;
  subcellShift.row(5) << 1., 0., 1.;
  subcellShift.row(6) << 0., 1., 1.;
  subcellShift.row(7) << 1., 1., 1.;
  for ( int i=0; i<8; i++) {
    for (int j :  {0,1,2}) {
      subcellShift(i, j) *= deltaGrid(j);
    }
  }

  // list of quadruplets identifying the vertices of the tetrahedra.
  // each number corresponds to the one number corresponds to the (interal) coordinate
  vertices.resize(6,4);
  vertices.row(0) << 0, 1, 2, 5;
  vertices.row(1) << 0, 2, 4, 5;
  vertices.row(2) << 2, 4, 5, 6;
  vertices.row(3) << 2, 5, 6, 7;
  vertices.row(4) << 2, 3, 5, 7;
  vertices.row(5) << 1, 2, 3, 5;
}

double TetrahedronDeltaFunction::getDOS(const double &energy) {
  // initialize tetrahedron weight
  double weight = 0.;
  for (long is : fullBandStructure.irrStateIterator()) {
    auto isIndex = StateIndex(is);
    long degeneracy = fullBandStructure.getRotationsStar(isIndex).size();
    weight += getSmearing(energy, isIndex) * degeneracy;
  }
  weight /= fullBandStructure.getNumPoints(true);
  return weight;
}

double TetrahedronDeltaFunction::getSmearing(const double &energy,
                                             StateIndex &is) {
  auto t = fullBandStructure.getIndex(is);
  long ik = std::get<0>(t).get();
  auto ibIndex = std::get<1>(t);

  auto fullPoints = fullBandStructure.getPoints();
  // if the mesh is uniform, each kpoint belongs to 6 tetrahedra
  auto kCoords = fullPoints.getPointCoords(ik, Points::crystalCoords);

  // in this tensor, we store the offset to find the vertices of the
  // tetrahedrons to which the current point belongs to
  // (6: number of tetrahedra, 4: number of vertices, 3: cartesian)
  // multiplying this by the vector Delta k of the grid, we can reconstruct
  // all the points of the tetrahedron.

  Eigen::MatrixXd kVectorsSubcell(8,3);
  for ( int i=0; i<8; i++) {
    Eigen::Vector3d x = subcellShift.row(i);
    kVectorsSubcell.row(i) = kCoords + x;
  }

  Eigen::VectorXi ikSubcellIndices(8);
  ikSubcellIndices(0) = ik;
  for ( int i=1; i<8; i++) {
    ikSubcellIndices(i) = fullPoints.getIndex(kVectorsSubcell.row(i));
  }

  Eigen::VectorXd energies(8);
  for ( int i=1; i<8; i++) {
    long is1 = fullBandStructure.getIndex(WavevectorIndex(ikSubcellIndices(i)),ibIndex);
    StateIndex is1Idx(is1);
    energies(i) = fullBandStructure.getEnergy(is1Idx);
  }

  // initialize tetrahedron weight
  double numTetra = 0.;
  double weight = 0.;
  for (int iTetra=0; iTetra<6; iTetra++) {
    std::vector<double> tmp(4);
    tmp[0] = energies(vertices(iTetra,0));
    tmp[1] = energies(vertices(iTetra,1));
    tmp[2] = energies(vertices(iTetra,2));
    tmp[3] = energies(vertices(iTetra,3));
    std::sort(tmp.begin(), tmp.end());
    double e1 = tmp[0];
    double e2 = tmp[1];
    double e3 = tmp[2];
    double e4 = tmp[3];

    // We follow Eq. B6 of Lambin and Vigneron PRB 29 3430 (1984)

    double cnE = 0.;
    if (e1 <= energy && energy <= e2) {
      if ( e2==e1 || e3==e1 || e4==e1 ) {
        cnE = 0.;
      } else {
        cnE = 3. * (energy-e1)*(energy-e1) / (e2-e1) / (e3-e1) / (e4-e1);
      }
    } else if (e2 <= energy && energy <= e3) {

      cnE = 0.;
      if ( e4==e2 || e3==e2 || e3==e1 ) {
        cnE += 0.;
      } else {
        cnE += (e3-energy)*(energy-e2) / (e4-e2) / (e3-e2) / (e3-e1);
      }
      if ( e4==e1 || e4==e2 || e3==e1 ) {
        cnE += 0.;
      } else {
        cnE += (e4-energy)*(energy-e1) / (e4-e1) / (e4-e2) / (e3-e1);
      }
      cnE *= 3.;

    } else if (e3 <= energy && energy <= e4) {
      if ( e4==e1 || e4==e2 || e4==e3 ) {
        cnE = 0.;
      } else {
        cnE = 3. * (e4-energy)*(e4-energy) / (e4-e1) / (e4-e2) / (e4-e3);
      }
    }

    // exception
    if ((e1 == e2) && (e1 == e3) && (e1 == e4) & (energy == e1)) {
      cnE = 0.25;
    }

    numTetra += 1.;

    weight += cnE;
  } // loop over all tetrahedra

  // Zero out extremely small weights
  if (weight < 1.0e-12) {
    weight = 0.;
  }

  // Normalize by number of tetrahedra and the vertices
  weight /= numTetra;
  return weight;
}

double TetrahedronDeltaFunction::getSmearing(const double &energy,
                                             const Eigen::Vector3d &velocity) {
  (void)energy;
  (void)velocity;
  Error e("TetrahedronDeltaFunction getSmearing1 not implemented");
  return 1.;
}
