#ifndef OBSERVABLE_H
#define OBSERVABLE_H

#include "eigen.h"
#include "context.h"
#include "vector_bte.h"

class Observable {
public:
	 // constructor
	Observable(StatisticsSweep & statisticsSweep_, Crystal & crystal_);
	// copy constructor
	Observable(const Observable & that);
	// copy assignment
	Observable & operator = (const Observable & that);
	// diff operator overload
	Observable operator - (const Observable & that);

	Eigen::VectorXd getNorm();
protected:
	StatisticsSweep & statisticsSweep;
	Crystal & crystal;

	long numChemPots;
	long numTemps;
	long dimensionality = 3;
	long numCalcs;

	virtual int whichType();
	const int isScalar = 0;
	const int isVector = 1;
	const int is2Tensor = 2;
	const int is4Tensor = 3;

	Eigen::VectorXd scalar; // e.g. specific heat
	Eigen::MatrixXd vectord; // e.g. sound speed
	Eigen::Tensor<double,3> tensordxd; // e.g. conductivity
	Eigen::Tensor<double,5> tensordxdxdxd; // e.g. viscosity

	long glob2Loc(const ChemPotIndex & imu, const TempIndex & it);
	std::tuple<ChemPotIndex,TempIndex> loc2Glob(const long & i);

	void baseOperatorMinus(Observable & newObservable, const Observable &that);

//	void calcFromPopulation(VectorBTE & n);
};

#endif
