#ifndef PHONONVISCOSITY_H
#define PHONONVISCOSITY_H

#include "observable.h"
#include "drift.h"
#include "ph_scattering.h"

/** Object to compute and store the phonon viscosity.
 * TODO: test if electron viscosity works with this class as well.
 */
class PhononViscosity : public Observable {
public:
	/** Object constructor. Simply stores references to the inputs
	 * @param statisticsSweep: object with info on the temperature and chemical
	 * potential loops
	 * @param crystal: object with the crystal structure
	 * @param FullBandStructure: object with the quasiparticle energies and
	 * velocities computed on a mesh of wavevectors.
	 */
	PhononViscosity(StatisticsSweep & statisticsSweep_, Crystal & crystal_,
			FullBandStructure<FullPoints> & bandStructure_);

	/** copy constructor
	 */
	PhononViscosity(const PhononViscosity & that);

	/** copy assignment
	 */
	PhononViscosity & operator = (const PhononViscosity & that);

	/** Compute the viscosity within the relaxation time approximation
	 * Stores it internally.
	 * @param n: the absolute phonon population.
	 */
	virtual void calcRTA(VectorBTE & n);

	/** Prints the viscosity to screen for the user.
	 */
	void print();

	/** Computes the viscosity from the scattering matrix eigenvectors.
	 * Following Simoncelli PRX 2020. Stores it internally.
	 * @param vector0: VectorBTE object with the energy-conserving eigenvector.
	 * @param relTimes: the VectorBTE object with relaxon relaxation times.
	 * @oaram sMatrix: the scattering matrix.
	 * @oaram eigenvectors: the eigenvectors of the scattering matrix above.
	 */
	void calcFromRelaxons(Vector0 & vector0, VectorBTE & relTimes,
			PhScatteringMatrix & sMatrix, Eigen::MatrixXd & eigenvectors);
protected:
	virtual int whichType();
	FullBandStructure<FullPoints> & bandStructure;
};

#endif
