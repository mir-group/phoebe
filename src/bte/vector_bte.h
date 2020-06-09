#ifndef VECTORBTE_H
#define VECTORBTE_H

#include "eigen.h"
#include "context.h"
#include "active_bandstructure.h"

class VectorBTE {
public:
	 // constructor
	VectorBTE(StatisticsSweep & statisticsSweep_,
			FullBandStructure<FullPoints> & bandStructure_,
			const long & dimensionality_=3);
	// copy constructor
	VectorBTE(const VectorBTE & that);
	// copy assignment
	VectorBTE & operator = (const VectorBTE & that);
	// scalar product
	Eigen::VectorXd dot(const VectorBTE & that);
	// product operator overload
	VectorBTE operator * (VectorBTE & that);
	// product operator overload
	VectorBTE operator * (const double & scalar);
	// product operator overload
	VectorBTE operator * (const Eigen::VectorXd & vector);
	// sum operator overload
	VectorBTE operator + (VectorBTE & that);
	// diff operator overload
	VectorBTE operator - (VectorBTE & that);
	// diff operator overload
	VectorBTE operator - ();
	// division operator overload
	VectorBTE operator / (VectorBTE & that);

	VectorBTE sqrt();
	VectorBTE reciprocal();

	void setConst(const double & constant);

	void canonical2Population();
	void population2Canonical();

//protected:
	StatisticsSweep & statisticsSweep;
	FullBandStructure<FullPoints> & bandStructure;
	long numCalcs;
	long numStates;
	long numChemPots;
	long numTemps;
	long dimensionality;
	long glob2Loc(const long & imu, const long & it, const long & idim);
	std::tuple<long,long,long> loc2Glob(const long & i);
	Eigen::MatrixXd data;
	friend class ScatteringMatrix; // this is also to remember that
	// if we change the index order of this class, we should check the
	// ScatteringMatrix implementations: they are high efficiency methods
	// so they need low-level access to the raw buffer

	VectorBTE baseOperator(VectorBTE & that, const int & operatorType);
	const int operatorSums = 0;
	const int operatorDivs = 1;
	const int operatorProd = 2;
	const int operatorDiff = 3;

	std::vector<long> excludeIndeces;
};

#endif
