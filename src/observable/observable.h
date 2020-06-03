#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include "eigen.h"
#include "context.h"
#include "vector_bte.h"

class Observable {
public:
	 // constructor
	Observable(Context & context_, Crystal & crystal_);
	// copy constructor
	Observable(const Observable & that);
	// copy assignment
	Observable & operator = (const Observable & that);
	// diff operator overload
	Observable operator - (const Observable & that);

	Eigen::VectorXd getNorm();
protected:
	Context & context;
	Crystal & crystal;

	long numChemPots;
	long numTemps;
	long dimensionality = 3; // to fix the dimensionality of quantities (if vectors)
	long numCalcs;

	int type;
	const int isScalar = 0;
	const int isVector = 1;
	const int is2Tensor = 2;
	const int is4Tensor = 3;

	Eigen::VectorXd scalar; // e.g. specific heat
	Eigen::MatrixXd vectord; // e.g. sound speed
	Eigen::Tensor<double,3> tensordxd; // e.g. conductivity
	Eigen::Tensor<double,5> tensordxdxdxd; // e.g. viscosity

	long glob2Loc(const long & imu, const long & it);
	std::tuple<long,long> loc2Glob(const long & i);

//	void calcFromPopulation(VectorBTE & n);
};

#endif
