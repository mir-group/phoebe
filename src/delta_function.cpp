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

double GaussianDeltaFunction::getSmearing(const double &energy, const long &iq,
                                          const long &ib) {
  (void)energy;
  (void)iq;
  (void)ib;
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
                                                  const long &iq,
                                                  const long &ib) {
  (void)energy;
  (void)iq;
  (void)ib;
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
  if (grid(0) == 1 || grid(1) == 1 || grid(2) == 1) {
    Error e("Tetrahedron method with k-grid dimensionality<3 not "
            "supported");
  }

  // number of grid points (wavevectors)
  long numPoints = fullPoints.getNumPoints();
  long numBands = fullBandStructure.getNumBands();
  // note: the code will fail at numBands if we are using ActiveBandStructure
  // this is intentional, as the implementation of tetrahedra works only
  // assuming a full grid of wavevectors.

  // Number of tetrahedra
  numTetra = 6 * numPoints;
  // Allocate tetrahedron data holders
  tetrahedra = Eigen::MatrixXi::Zero(numTetra, 4);
  qToTetCount = Eigen::VectorXi::Zero(numPoints);
  qToTet = Eigen::MatrixXi::Zero(numPoints, 24);

  // Label the vertices of each tetrahedron in a subcell
  Eigen::MatrixXi verticesLabels(6, 4);
  verticesLabels << 0, 1, 2, 5, 0, 2, 4, 5, 2, 4, 5, 6, 2, 5, 6, 7, 2, 3, 5, 7,
      1, 2, 3, 5;

  // 8 corners of a subcell
  Eigen::MatrixXi subcellCorners(8, 3);

  for (long iq = 0; iq < fullPoints.getNumPoints(); iq++) {
    // point is a vector with coordinates between 0 and 1
    Eigen::Vector3d point =
        fullPoints.getPointCoords(iq, Points::crystalCoords);
    // scale it to integers between 0 and grid size
    point(0) *= grid(0);
    point(1) *= grid(1);
    point(2) *= grid(2);

    int i = int(point(0));
    int j = int(point(1));
    int k = int(point(2));
    int ip1 = mod((i + 1), grid(0));
    int jp1 = mod((j + 1), grid(1));
    int kp1 = mod((k + 1), grid(2));

    subcellCorners << i, j, k, ip1, j, k, i, jp1, k, ip1, jp1, k, i, j, kp1,
        ip1, j, kp1, i, jp1, kp1, ip1, jp1, kp1;

    for (int it = 0; it < 6; it++) {   // over 6 tetrahedra
      for (int iv = 0; iv < 4; iv++) { // over 4 vertices
        // Grab a label
        long aux = verticesLabels(it, iv);
        // Grab a corner of subcell
        point(0) = double(subcellCorners(aux, 0)) / grid(0);
        point(1) = double(subcellCorners(aux, 1)) / grid(1);
        point(2) = double(subcellCorners(aux, 2)) / grid(2);
        // Get combined index of corner
        long aux2 = fullPoints.getIndex(point);
        // Save corner as a tetrahedron vertex
        tetrahedra(iq, iv) = aux2;
        // Save mapping of a wave vector index
        // to the ordered pair (tetrahedron,vertex)
        qToTetCount(aux2) = qToTetCount(aux2) + 1;
        qToTet(aux2, qToTetCount(aux2) - 1) = iq;
      }
    }
  }
  /**
   * Fill all tetrahedra with the eigenvalues.
   *
   * Method for filling the tetrahedra with the eigenvalues for
   * all polarizations. For eigenvalues are sorted along the vertex.
   */

  // Internal variables
  std::vector<double> temp(4);

  // Allocate tetraEigVals
  tetraEigVals = Eigen::Tensor<double, 3>(numTetra, numBands, 4);

  for (long it = 0; it < numTetra; it++) { // over tetrahedra
    // over bands
    for (int ib = 0; ib < numBands; ib++) {
      for (int iv = 0; iv < 4; iv++) { // over vertices
        // Index of wave vector
        long ik = tetrahedra(it, iv);

        // Fill tetrahedron vertex with the band energy
        auto is =
            fullBandStructure.getIndex(WavevectorIndex(ik), BandIndex(ib));
        double energy = fullBandStructure.getEnergy(is);
        tetraEigVals(it, ib, iv) = energy;
        temp[iv] = energy; // save for later
      }

      // sort energies in the vertex
      std::sort(temp.begin(), temp.end());
      // refill tetrahedron vertex
      for (int iv = 0; iv < 4; iv++) {
        tetraEigVals(it, ib, iv) = temp[iv];
      }
    }
  }
}

double TetrahedronDeltaFunction::getDOS(const double &energy) {
  // initialize tetrahedron weight
  double weight = 0.;
  for (int iq = 0; iq < fullBandStructure.getNumPoints(); iq++) {
    for (int ib = 0; ib < fullBandStructure.getNumBands(); ib++) {
      weight += getWeight(energy, iq, ib);
    }
  }
  weight /= fullBandStructure.getNumPoints(true);
  return weight;
}

double TetrahedronDeltaFunction::getSmearing(const double &energy,
                                             const long &iq, const long &ib) {
  // initialize tetrahedron weight
  return getWeight(energy, iq, ib);
}

double TetrahedronDeltaFunction::getWeight(const double &energy, const long &iq,
                                           const long &ib) {

  // initialize tetrahedron weight
  double weight = 0.;

  // loop on the number of tetrahedra in which the wave vector belongs
  for (long i = 0; i < qToTetCount(iq); i++) { // over all tetrahedra
    long it = qToTet(iq, i); // get index of tetrahedron

    // Sorted energies at the 4 vertices
    double e1 = tetraEigVals(it, ib, 0);
    double e2 = tetraEigVals(it, ib, 1);
    double e3 = tetraEigVals(it, ib, 2);
    double e4 = tetraEigVals(it, ib, 3);

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

    weight += cnE;
  } // loop over all tetrahedra

  // Zero out extremely small weights
  if (weight < 1.0e-12) {
    weight = 0.;
  }

  // Normalize by number of tetrahedra and the vertices
  weight /= 6. * 4.;
  return weight;
}

double TetrahedronDeltaFunction::getSmearing(const double &energy,
                                             const Eigen::Vector3d &velocity) {
  (void)energy;
  (void)velocity;
  Error e("TetrahedronDeltaFunction getSmearing1 not implemented");
  return 1.;
}
