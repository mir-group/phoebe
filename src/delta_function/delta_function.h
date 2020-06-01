#ifndef DELTAF_H
#define DELTAF_H

#include "eigen.h"
#include "points.h"
#include "bandstructure.h"
#include "context.h"

class DeltaFunction {
public:
	virtual ~DeltaFunction();
	// here a smearing factory
	static DeltaFunction * smearingFactory(Context & context,
			FullBandStructure<FullPoints> & fullBandStructure);
	static const int gaussian = 0;
	static const int adaptiveGaussian = 1;
	static const int tetrahedron = 2;
	const int id = -1;

	virtual double getSmearing(const double & energy,
			const Eigen::Vector3d & velocity=Eigen::Vector3d::Zero()) = 0;
	virtual double getSmearing(const double & energy, const long & iq,
			const long &ib) = 0;
};

class GaussianDeltaFunction : public DeltaFunction {
public:
	GaussianDeltaFunction(Context & context); // context to get amplitude
	virtual double getSmearing(const double & energy,
			const Eigen::Vector3d & velocity=Eigen::Vector3d::Zero());
	virtual double getSmearing(const double & energy, const long & iq,
			const long &ib);
	const int id = 0;
protected:
	double inverseWidth;
	double prefactor;
};

class AdaptiveGaussianDeltaFunction : public DeltaFunction {
public:
	AdaptiveGaussianDeltaFunction(FullBandStructure<FullPoints>&bandStructure);
	virtual double getSmearing(const double & energy,
			const Eigen::Vector3d & velocity=Eigen::Vector3d::Zero());
	virtual double getSmearing(const double & energy, const long & iq,
			const long &ib);
	const int id = 1;
protected:
	const double smearingCutoff = 1.0e-8;
	const double prefactor = 1.;
	Eigen::Matrix3d qTensor;
};

/**
 * Class for approximating the Delta function with the tetrahedron method
 */
class TetrahedronDeltaFunction : public DeltaFunction {
public:
	const int id = 2;

	/**
	 * Form all tetrahedra for 3D wave vector mesh.
	 *
	 * Method for creating and enumerating all the tetrahedra
	 * for a given 3D mesh of wave vectors following Fig. 5 of
	 * Bloechl, Jepsen and Andersen prb 49.23 (1994): 16223.
	 *
	 * @param[in] grid: the mesh points along the three lattice vectors.
	 *
	 */
	TetrahedronDeltaFunction(FullBandStructure<FullPoints>&fullBandStructure_);

	/**
	 * Calculate tetrehedron weight.
	 *
	 * Method for calculating the tetrahedron weight (normalized by the number of
	 * tetrahedra) for given wave vector and polarization following Lambin and
	 * Vigneron prb 29.6 (1984): 3430.
	 *
	 * @param[in] energy Energy of mode.
	 */
	double getDOS(const double & energy);

	/**
	 * Calculate tetrehedron weight.
	 *
	 * Method for calculating the tetrahedron weight (normalized by the number of
	 * tetrahedra) for given wave vector and polarization following Lambin and
	 * Vigneron prb 29.6 (1984): 3430.
	 *
	 * @param[in] energy Energy of mode.
	 * @param[in] State: state at which the tetrahedron is computed.
	 * @returns The tetrahedron weight.
	 *
	 */
	virtual double getSmearing(const double & energy, const long & iq,
			const long & ib);
	virtual double getSmearing(const double & energy,
			const Eigen::Vector3d & velocity=Eigen::Vector3d::Zero());
protected:
	FullBandStructure<FullPoints> & fullBandStructure;

	/** Number of tetrahedra. */
	long numTetra;
	/** Holder for the indices of the vertices of of each tetrahedron. */
	Eigen::MatrixXi tetrahedra;
	/** Count of how many tetrahedra wave vector belongs to. */
	Eigen::VectorXi qToTetCount;
	/** Mapping of a wave vector to a tetrahedron. */
	Eigen::Tensor<long,3> qToTet;
	/** Holder for the eigenvalues. */
	Eigen::Tensor<double,3> tetraEigVals;

	double getWeight(const double & energy, const long & iq, const long & ib);

};

#endif
