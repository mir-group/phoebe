#ifndef DELTAF_H
#define DELTAF_H

#include "eigen.h"
#include "points.h"
#include "bandstructure.h"
#include "context.h"

/** Base class for the approximations to the Delta Function.
 * Currently used for density of states calculation or for transport
 * calculations (scattering rates require a dirac-delta).
 */
class DeltaFunction {
public:
	/** Destructor for the DeltaFunction method.
	 * Needed to destruct the pointer returned by the factory method.
	 */
	virtual ~DeltaFunction();

	/** A smearing factory for loading a specific approximation to the
	 * dirac-delta. Returns a pointer. Therefore, it should be deleted
	 * afterwards.
	 * @param context: object with the user input, containing the choice of
	 * smearing.
	 * @param fullBandStructure: the object containing the bandstructure
	 * computed on a mesh of wavevectors in the Brillouin zone.
	 */
	static DeltaFunction * smearingFactory(Context & context,
			FullBandStructure<FullPoints> & fullBandStructure);

	// These are integers identifying the smearing selection of the user.
	static const int gaussian = 0;
	static const int adaptiveGaussian = 1;
	static const int tetrahedron = 2;

	/** Method for obtaining a smearing value (approximating \delta).
	 * @param energy: the energy difference of the dirac delta
	 * @param velocity optional: a vector of velocity (used for adaptive
	 * smearing schemes).
	 */
	virtual double getSmearing(const double & energy,
			const Eigen::Vector3d & velocity=Eigen::Vector3d::Zero()) = 0;

	/** Method for obtaining a smearing value (approximating \delta).
	 * Overloads the previous one, and is used by the tetrahedron method.
	 * @param energy: the energy difference of the dirac delta.
	 * @param iq: wavevector index.
	 * @param ib: band index.
	 */
	virtual double getSmearing(const double & energy, const long & iq,
			const long &ib) = 0;

	/** Method to identify which kind of smearing is being used.
	 * Returns a int value between gaussian, adaptiveGaussian, tetrahedron.
	 */
	virtual int getType();
private:
	int id = -1;
};

/** Gaussian smearing scheme. Simply replace the dirac delta with a gaussian.
 * The width of such gaussian is specified by the user.
 */
class GaussianDeltaFunction : public DeltaFunction {
public:
	/** Constructor of the gaussian smearing scheme.
	 * @param context: object with user input.
	 */
	GaussianDeltaFunction(Context & context); // context to get amplitude

	/** Method to obtain the value of smearing.
	 * @param energy: the energy difference.
	 * @param[optional] velocity: ignored parameter.
	 * @return smearing: the approximation to the dirac-delta.
	 */
	virtual double getSmearing(const double & energy,
			const Eigen::Vector3d & velocity=Eigen::Vector3d::Zero());

	/** phantom method that should not be used. Will throw an error.
	 */
	virtual double getSmearing(const double & energy, const long & iq,
			const long &ib);

	/** returns an integer identifying this class as AdaptiveGaussian
	 * @return int: id.
	 */
	virtual int getType();
protected:
	int id = DeltaFunction::gaussian;
	double inverseWidth;
	double prefactor;
};

/** Adaptive smearing scheme. Uses the group velocity to adjust the width
 * of a gaussian approximation to the dirac-delta.
 * See Ref doi: 10.1103/PhysRevB.75.195121
 */
class AdaptiveGaussianDeltaFunction : public DeltaFunction {
public:
	/** Constructor of the adaptive gaussian smearing scheme.
	 * @param bandStructure: mainly to see what mesh of points and crystal is
	 * being used, and to prepare a suitable scaling of velocities.
	 */
	AdaptiveGaussianDeltaFunction(FullBandStructure<FullPoints>&bandStructure);

	/** Method to obtain the value of smearing.
	 * @param energy: the energy difference.
	 * @param velocity; the velocity (difference) of the quasiparticles.
	 * @return smearing: the approximation to the dirac-delta
	 */
	virtual double getSmearing(const double & energy,
			const Eigen::Vector3d & velocity=Eigen::Vector3d::Zero());

	/** phantom method that should not be used. Will throw an error.
	 */
	virtual double getSmearing(const double & energy, const long & iq,
			const long &ib);

	/** returns an integer identifying this class as AdaptiveGaussian
	 * @return int: id.
	 */
	virtual int getType();
protected:
	int id = DeltaFunction::adaptiveGaussian;

	// we use a cutoff to treat specially cases where the smearing approaches 0
	const double smearingCutoff = 1.0e-8;

	const double prefactor = 1.; // could be exposed to the user
	Eigen::Matrix3d qTensor; // to normalize by the number of wavevectors.

	// since we truncate the smearing above 2sigma, we use erf to correct
	// results
	const double erf2 = 0.9953222650189527;
};

/** Class for approximating the Delta function with the tetrahedron method
 */
class TetrahedronDeltaFunction : public DeltaFunction {
public:
	virtual int getType();

	/** Constructor method.
	 * Forms all tetrahedra for the 3D wave vector mesh.
	 * Lower dimensionality is not currently supported.
	 *
	 * Method for creating and enumerating all the tetrahedra
	 * for a given 3D mesh of wave vectors following Fig. 5 of
	 * Bloechl, Jepsen and Andersen prb 49.23 (1994): 16223.
	 *
	 * @param[in] bandStructure: the bandstructure on which the dirac delta
	 * will be computed.
	 */
	TetrahedronDeltaFunction(FullBandStructure<FullPoints>&fullBandStructure_);

	/** Calculate the total tetrehedron weight for all states at given energy.
	 *
	 * Method for calculating the tetrahedron weight for given wave vector and
	 * polarization following Lambin and Vigneron prb 29.6 (1984): 3430.
	 *
	 * @param[in] energy Energy of mode.
	 * @return dos: density of states for the given value of energy
	 */
	double getDOS(const double & energy);

	/** Calculate tetrehedron weight.
	 *
	 * Method for calculating the tetrahedron weight for given wave vector and
	 * polarization following Lambin and Vigneron prb 29.6 (1984): 3430.
	 *
	 * @param[in] energy Energy of mode.
	 * @param[in] State: state at which the tetrahedron is computed.
	 * @returns weight: the tetrahedron weight approximating the dirac delta.
	 */
	virtual double getSmearing(const double & energy, const long & iq,
			const long & ib);

	/** overload abstract base method, but simply raises an error if it's
	 * called. One should really use the other getSmearing method.
	 */
	virtual double getSmearing(const double & energy,
			const Eigen::Vector3d & velocity=Eigen::Vector3d::Zero());
protected:
	FullBandStructure<FullPoints> & fullBandStructure;
	int id = DeltaFunction::tetrahedron;

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
