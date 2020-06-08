#ifndef STATISTICS_H
#define STATISTICS_H

#include "eigen.h"

class Statistics {
public:
	static const int electron = -1;
	static const int phonon   = -2;

	double getPopulation(const double & energy, const double & temperature,
			const double & chemicalPotential=0.);

	double getDndt(const double & energy, const double & temperature,
			const double & chemicalPotential=0.);
	double getDnde(const double & energy, const double & temperature,
			const double & chemicalPotential=0.);

	Eigen::VectorXd getDndt(const Eigen::VectorXd & energies,
			const double & temperature, const double & chemicalPotential=0.);
	Eigen::VectorXd getDnde(const Eigen::VectorXd & energies,
			const double & temperature, const double & chemicalPotential=0.);

	/** Returns bose(bose+1) for bosons, fermi(1-fermi) for fermions
	 *
	 */
	double getPopPopPm1(const double & energy,
			const double & temperature, const double & chemicalPotential=0.);

//	Statistics() = default; // default constructor, to be used only for
	// temporary object initializations

	Statistics(int statistics_);
	Statistics(const Statistics & obj);// copy constructor
	Statistics & operator=(const Statistics & obj); // copy assignment operator
	~Statistics();

	bool isFermi();
	bool isBose();
	bool isElectron();
	bool isPhonon();
private:
	static const int fermi = 0;
	static const int bose = 1;

	int statistics;
	int particle;

	Eigen::VectorXd temperatures;
	Eigen::VectorXd chemicalPotentials;
	Eigen::VectorXd dopings;
};

#endif
