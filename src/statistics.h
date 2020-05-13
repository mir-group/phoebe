#ifndef STATISTICS_H
#define STATISTICS_H

#include "eigen.h"

class Statistics {
public:
	static const int fermi = 0;
	static const int bose = 1;

	double getPopulation(double energy, double temperature,
			double chemicalPotential=0.);

	double getDndt(double energy, double temperature,
			double chemicalPotential=0.);
	double getDnde(double energy, double temperature,
			double chemicalPotential=0.);

	Eigen::MatrixXd getDndt(Eigen::VectorXd energies,
			Eigen::VectorXd temperatures, double chemicalPotential=0.);
	Eigen::MatrixXd getDnde(Eigen::VectorXd energies,
			Eigen::VectorXd temperatures, double chemicalPotential=0.);

//	Statistics() = default; // default constructor, to be used only for
	// temporary object initializations

	Statistics(int statistics_);

	Statistics(const Statistics &obj);// copy constructor
	Statistics & operator=(const Statistics & obj); // copy assignment operator

	bool isFermi();
	bool isBose();
private:
	int statistics;
};

#endif
