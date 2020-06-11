#ifndef STATISTICS_H
#define STATISTICS_H

#include "eigen.h"

/** This object decides whether we are dealing with phonons or electrons.
 * Additionally, it contains the functionality to compute equilibrium
 * occupation numbers, i.e. Fermi-Dirac or Bose-Einstein distributions, their
 * derivatives w.r.t. temperature or energy, and factors N(N+1) or F(1-F).
 */
class Statistics {
public:
	static const int electron = -1;
	static const int phonon   = -2;

	/** Constructor of Statistics. Saves which kind of particle we are dealing
	 * with.
	 * @param statistics: an integer equal to Statistics::electron or
	 * Statistics::phonon.
	 */
	Statistics(int statistics_);

	/** Copy constructor
	 */
	Statistics(const Statistics & obj);

	/** Copy assignment operator
	 */
	Statistics & operator=(const Statistics & obj);

	/** Returns either a Bose--Einstein or a Fermi--Dirac distribution,
	 * depending on the value of "statistics".
	 * @param energy: value of quasiparticle energy.
	 * @param temperature: value of temperature.
	 * @param chemicalPotential: 0 by default, set a value for electrons.
	 * @return n: either the Bose or Fermi population.
	 */
	double getPopulation(const double & energy, const double & temperature,
			const double & chemicalPotential=0.);

	/** Returns dn/dT, with T temperature, and n being either a Bose--Einstein
	 * or a Fermi--Dirac distribution, depending on the value of "statistics"
	 * @param energy: value of quasiparticle energy.
	 * @param temperature: value of temperature.
	 * @param chemicalPotential: 0 by default, set a value for electrons.
	 * @return dndT: derivative wrt temperature of the equilibrium distribution
	 */
	double getDndt(const double & energy, const double & temperature,
			const double & chemicalPotential=0.);

	/** Returns dn/dE, with E the particle energy, and n being either a
	 * Bose or Fermi distribution, depending on the value of "statistics".
	 * @param energy: value of quasiparticle energy.
	 * @param temperature: value of temperature.
	 * @param chemicalPotential: 0 by default, set a value for electrons.
	 * @return dndE: derivative wrt QP energy of the equilibrium distribution.
	 */
	double getDnde(const double & energy, const double & temperature,
			const double & chemicalPotential=0.);

	/** Returns bose(bose+1) for bosons, fermi(1-fermi) for fermions.
	 * @param energy: value of quasiparticle energy.
	 * @param temperature: value of temperature.
	 * @param chemicalPotential: 0 by default, set a value for electrons.
	 * @return popPopPm1: bose(bose+1) for bosons, fermi(1-fermi) for fermions.
	 */
	double getPopPopPm1(const double & energy,
			const double & temperature, const double & chemicalPotential=0.);

	/** Method to check if the particle is a fermion.
	 */
	bool isFermi();

	/** Method to check if the particle is a boson.
	 */
	bool isBose();

	/** Method to check if the particle is an electron.
	 */
	bool isElectron();

	/** Method to check if the particle is a phonon.
	 */
	bool isPhonon();
private:
	static const int fermi = 0;
	static const int bose = 1;

	// identifies whether the instance is a fermion of a boson
	int statistics;
	// identifies whether we have an electron or a phonon
	// Note: we keep particle != statistics because we might add particles
	int particle;
};

#endif
