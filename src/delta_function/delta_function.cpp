#include "delta_function.h"
#include "exceptions.h"
#include "utilities.h"
#include "constants.h"
#include "context.h"

DeltaFunction::~DeltaFunction() {
}

int DeltaFunction::getType() {
	return id;
}

int GaussianDeltaFunction::getType() {
	return id;
}

int AdaptiveGaussianDeltaFunction::getType() {
	return id;
}

int TetrahedronDeltaFunction::getType() {
	return id;
}

// app factory
DeltaFunction * DeltaFunction::smearingFactory(Context & context,
		FullBandStructure<FullPoints> & fullBandStructure) {
	auto choice = context.getSmearingMethod();
	if ( choice == gaussian ) {
		return new GaussianDeltaFunction(context);
	} else if ( choice == adaptiveGaussian ) {
		return new AdaptiveGaussianDeltaFunction(fullBandStructure);
	} else if ( choice == tetrahedron ) {
		return new TetrahedronDeltaFunction(fullBandStructure);
	} else {
		Error e("Unrecognized smearing choice");
		return nullptr;
	}
}

GaussianDeltaFunction::GaussianDeltaFunction(Context & context) {
	inverseWidth = 1. / context.getSmearingWidth();
	prefactor = 1. / context.getSmearingWidth() / sqrt(pi);
}

double GaussianDeltaFunction::getSmearing(const double & energy,
		const Eigen::Vector3d & velocity) {
	(void) velocity;
	double x = energy * inverseWidth;
	return prefactor * exp( - x*x );
}

double GaussianDeltaFunction::getSmearing(const double & energy,
		const long & iq, const long &ib) {
	(void) energy;
	(void) iq;
	(void) ib;
	Error e("GaussianDeltaFunction::getSmearing2 not implemented");
	return 1.;
}


AdaptiveGaussianDeltaFunction::AdaptiveGaussianDeltaFunction(
		FullBandStructure<FullPoints> & bandStructure) {
	auto [mesh,offset] = bandStructure.getPoints().getMesh();
	qTensor = bandStructure.getPoints().getCrystal().getReciprocalUnitCell();
	qTensor.row(0) /= mesh(0);
	qTensor.row(1) /= mesh(1);
	qTensor.row(2) /= mesh(2);
}

double AdaptiveGaussianDeltaFunction::getSmearing(const double & energy,
		const Eigen::Vector3d & velocity) {

	if ( velocity.norm() == 0. && energy == 0. ) {
		// in this case, velocities are parallel, there shouldn't be
		// scattering unless energy is strictly conserved
		return 1.;
	}

	double sigma = 0.;
    for ( int i : {0,1,2} ) {
        sigma += pow( qTensor.row(i).dot(velocity) , 2);
    }
    sigma = prefactor * sqrt(sigma / 6.);

    if ( sigma == 0. ) return 0.;

	if ( abs(energy) > 2. * sigma ) return 0.;
	double x = energy / sigma;
	// note: the factor ERF_2 corrects for the cutoff at 2*sigma
	return exp( - x*x ) / sqrtPi / sigma / erf2;
}

double AdaptiveGaussianDeltaFunction::getSmearing(const double & energy,
		const long & iq, const long & ib) {
	(void) energy;
	(void) iq;
	(void) ib;
	Error e("AdaptiveGaussianDeltaFunction::getSmearing2 not implemented");
	return 1.;
}

TetrahedronDeltaFunction::TetrahedronDeltaFunction(
		FullBandStructure<FullPoints> & fullBandStructure_) :
		fullBandStructure(fullBandStructure_) {

	auto fullPoints = fullBandStructure.getPoints();
	auto [grid, offset] = fullPoints.getMesh();
	if ( offset.norm() > 0. ) {
		Error e("We didnt' implement tetrahedra with offsets", 1);
	}

	// number of grid points (wavevectors)
	long numPoints = fullPoints.getNumPoints();
	long numBands = fullBandStructure.getNumBands();

	// Number of tetrahedra
	numTetra = 6 * numPoints;
	// Allocate tetrahedron data holders
	tetrahedra = Eigen::MatrixXi::Zero(numTetra,4);
	qToTetCount = Eigen::VectorXi::Zero(numPoints);
	qToTet = Eigen::Tensor<long,3>(numPoints,24,2);
	qToTet.setZero();

	// Label the vertices of each tetrahedron in a subcell
	Eigen::MatrixXi verticesLabels(6,4);
	verticesLabels <<
			0,1,2,5,
			0,2,4,5,
			2,4,5,6,
			2,5,6,7,
			2,3,5,7,
			1,2,3,5;

	// 8 corners of a subcell
	Eigen::MatrixXi subcellCorners(8,3);

	for ( long iq=0; iq<fullPoints.getNumPoints(); iq++ ) {
		// point is a vector with coordinates between 0 and 1
		Eigen::Vector3d point = fullPoints.getPointCoords(iq,
				Points::crystalCoords);
		// scale it to integers between 0 and grid size
		point(0) *= grid(0);
		point(1) *= grid(1);
		point(2) *= grid(2);

		long i = long(point(0));
		long j = long(point(1));
		long k = long(point(2));
		long ip1 = mod((i+1) , grid(0));
		long jp1 = mod((j+1) , grid(1));
		long kp1 = mod((k+1) , grid(2));

		subcellCorners <<
				i,j,k,
				ip1,j,k,
				i,jp1,k,
				ip1,jp1,k,
				i,j,kp1,
				ip1,j,kp1,
				i,jp1,kp1,
				ip1,jp1,kp1;

		for ( long it = 0; it < 6; it++ ) { //over 6 tetrahedra
			for ( long iv = 0; iv < 4; iv++ ) { //over 4 vertices
				// Grab a label
				long aux = verticesLabels(it,iv);
				// Grab a corner of subcell
				point(0) = double(subcellCorners(aux, 0)) / grid(0);
				point(1) = double(subcellCorners(aux, 1)) / grid(1);
				point(2) = double(subcellCorners(aux, 2)) / grid(2);
				// Get combined index of corner
				long aux2 = fullPoints.getIndex(point);
				// Save corner as a tetrahedron vertex
				tetrahedra(iq,iv) = aux2;
				// Save mapping of a wave vector index
				// to the ordered pair (tetrahedron,vertex)
				qToTetCount(aux2) = qToTetCount(aux2) + 1;
				qToTet(aux2, qToTetCount(aux2)-1, 0) = iq;
				qToTet(aux2, qToTetCount(aux2)-1, 1) = iv;
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
	tetraEigVals = Eigen::Tensor<double,3>(numTetra, numBands, 4);

	for ( long it = 0; it < numTetra; it++ ) { //over tetrahedra
		 //over bands
		for ( long ib = 0; ib < numBands; ib++ ) {
			for ( long iv = 0; iv < 4; iv++ ) { //over vertices
				// Index of wave vector
				long ik = tetrahedra(it,iv);

				// Fill tetrahedron vertex with the band energy
				Point point = fullBandStructure.getPoint(ik);
				State state = fullBandStructure.getState(point);
				double energy = state.getEnergy(ib);
				tetraEigVals(it,ib,iv) = energy;
				temp[iv] = energy; //save for later
			}

			//sort energies in the vertex
			std::sort(temp.begin(),temp.end());
			//refill tetrahedron vertex
			for ( long iv = 0; iv < 4; iv++) {
				tetraEigVals(it,ib,iv) = temp[iv];
			}
		}
	}
}

double TetrahedronDeltaFunction::getDOS(const double & energy) {
	// initialize tetrahedron weight
	double weight = 0.;

	for ( long iq=0; iq<fullBandStructure.getNumPoints(); iq++ ) {
		for ( long ib=0; ib<fullBandStructure.getNumBands(); ib++ ) {
			weight += getWeight(energy, iq, ib);
		}
	}
	return weight;
}

double TetrahedronDeltaFunction::getSmearing(const double & energy,
		const long & iq, const long & ib) {
	// initialize tetrahedron weight
	return getWeight(energy, iq, ib);
}

double TetrahedronDeltaFunction::getWeight(const double & energy, const long & iq,
		const long & ib) {

	// initialize tetrahedron weight
	double weight = 0.;

	// Internal variables
	double tmp = 0.0;

	// loop on the number of tetrahedra in which the wave vector belongs
	for ( long i = 0; i < qToTetCount(iq); i++ ) {//over all tetrahedra
		long it = qToTet(iq,i,0); //get index of tetrahedron
		long iv = qToTet(iq,i,1); //get index of vertex

		//Sorted energies at the 4 vertices
		double e1 = tetraEigVals(it,ib,0);
		double e2 = tetraEigVals(it,ib,1);
		double e3 = tetraEigVals(it,ib,2);
		double e4 = tetraEigVals(it,ib,3);

		//Define the shorthands
		// Refer to Lambin and Vigneron prb 29.6 (1984): 3430 to understand
		// what these mean.
		double e1e = e1 - energy;
		double e2e = e2 - energy;
		double e3e = e3 - energy;
		double e4e = e4 - energy;
		double e21 = e2 - e1;
		double e31 = e3 - e1;
		double e41 = e4 - e1;
		double e32 = e3 - e2;
		double e42 = e4 - e2;
		double e43 = e4 - e3;

		//Check the inequalities
		bool c1 = (e1 <= energy) && (energy <= e2);
		bool c2 = (e2 <= energy) && (energy <= e3);
		bool c3 = (e3 <= energy) && (energy <= e4);

		if( !((energy < e1) || (energy > e4)) ){
			if ( iv == 0 ) { //switch over 4 vertices
				if ( c1 ) {
					tmp = (e2e/e21 + e3e/e31 + e4e/e41)*pow(e1e,2)/e41/e31/e21;

					if ( e1 == e2 ) tmp = 0.0;
				} else if ( c2 ) {
					tmp = -0.5*(e3e/pow(e31,2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42
							+ e3e*e1e/e32/e41) + e4e/pow(e41,2)*(e4e*e1e/e42/e31
									+ e4e*e2e/e42/e32 + e3e*e1e/e31/e32));

					if ( e2 == e3 ) {
						tmp = -0.5*(e4e*e1e/e41/e42 + e1e/e41
								+ e4e/pow(e41,2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31));
					}
				} else if ( c3 ) {
					tmp = pow(e4e,3)/pow(e41,2)/e42/e43;

					if ( e3 == e4 ){
						tmp = pow(e4e,2)/pow(e41,2)/e42;
					}
				}
			}  else if ( iv == 1 ) {
				if ( c1 ) {
					tmp = -pow(e1e,3)/pow(e21,2)/e31/e41;

					if ( e1 == e2 ) tmp = 0.0;
				} else if ( c2 ) {
					tmp = -0.5*(e3e/pow(e32,2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41
							+ e3e*e1e/e31/e41)
							+ e4e/pow(e42,2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31
									+ e4e*e2e/e32/e41));

					if ( e2 == e3 ) {
						tmp = -0.5*(e4e/e42/e41
								+ e4e/pow(e42,2)*(e4e*e1e/e41/e31 + 1.0));
					}
				} else if ( c3 ) {
					tmp = pow(e4e,3)/e41/pow(e42,2)/e43;

					if ( e3 == e4 ) tmp = 0.0;
				}
			} else if ( iv == 2 ) {
				if ( c1 ) {
					tmp = -pow(e1e,3)/e21/pow(e31,2)/e41;

					if ( e1 == e2 ) tmp = 0.0;
				} else if ( c2 ) {
					tmp = 0.5*(e2e/pow(e32,2)*(e3e*e2e/e42/e31 + e4e*e2e/e42/e41
							+ e3e*e1e/e31/e41)
							+ e1e/pow(e31,2)*(e3e*e2e/e42/e32 + e4e*e1e/e41/e42
									+ e3e*e1e/e32/e41));

					if ( e2 == e3 ) {
						tmp = 0.5*(e4e/e42/e41 + e1e/e31/e41
								+ e1e/pow(e31,2)*(e4e*e1e/e41/e42 + e1e/e41));
					}
				} else if ( c3 ) {
					tmp = pow(e4e,3)/e41/e42/pow(e43,2);

					if ( e3 == e4 ) tmp = 0.0;
				}
			} else {
				if ( c1 ) {
					tmp = -pow(e1e,3)/e21/e31/pow(e41,2);

					if ( e1 == e2 ) tmp = 0.0;
				} else if ( c2 ) {
					tmp = 0.5*(e2e/pow(e42,2)*(e3e*e2e/e32/e31 + e4e*e1e/e41/e31
							+ e4e*e2e/e32/e41)
							+ e1e/pow(e41,2)*(e4e*e1e/e42/e31 + e4e*e2e/e42/e32
									+ e3e*e1e/e31/e32));

					if ( e2 == e3 ) {
						tmp = 0.5*e1e/pow(e41,2)*(e4e*e1e/e42/e31 + e4e/e42 + e1e/e31);
					}
				} else if ( c3 ) {
					tmp = -(e3e/e43 + e2e/e42 + e1e/e41)*pow(e4e,2)/e41/e42/e43;

					if ( e3 == e4 ) tmp = 0.0;
				}
			} // switch over 4 vertices

			if( (e1 == e2) && (e1 == e3) && (e1 == e4) & (energy == e1) ) tmp = 0.25;

			weight += tmp;
		} //!((energy < e1) || (energy > e4))
	} //over all tetrahedra

	// Zero out extremely small weights
	if ( weight < 1.0e-12 ) weight = 0.;

	// Normalize by number of tetrahedra
	weight /= double(numTetra);
	return weight;
}

double TetrahedronDeltaFunction::getSmearing(const double & energy,
		const Eigen::Vector3d & velocity) {
	(void) energy;
	(void) velocity;
	Error e("TetrahedronDeltaFunction getSmearing1 not implemented");
	return 1.;
}
