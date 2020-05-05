#ifndef VECTORBTE_H
#define VECTORBTE_H

#include "eigen.h"
#include "context.h"
#include "bandstructure.h"

//class Observable;

// base class
class VectorBTE {
public:
	 // constructor
	VectorBTE(Context & context_, ActiveBandStructure & activeBandStructure_);
	// copy constructor
	VectorBTE(const VectorBTE & that);
	// copy assignment
	VectorBTE & operator = (const VectorBTE & that);
	// product operator overload
	Eigen::VectorXd operator * (const VectorBTE & that);
	// product operator overload
	VectorBTE operator * (const double & scalar);
	// product operator overload
	VectorBTE operator * (const Eigen::VectorXd & vector);
	// sum operator overload
	VectorBTE operator + (const VectorBTE & that);
	// diff operator overload
	VectorBTE operator - (const VectorBTE & that);
	// diff operator overload
	VectorBTE operator - ();
	// division operator overload
	VectorBTE operator / (const VectorBTE & that);

	VectorBTE sqrt();
//protected:
	Context & context;
	ActiveBandStructure & activeBandStructure;
	long numChemPots;
	long numTemps;
	long dimensionality;
	long numCalcs;
	long numStates;
	long glob2Loc(long & imu, long & it, long & idim);
	std::tuple<long,long,long> loc2Glob(long & i);
	Eigen::MatrixXd data;
//	friend class Observable;
};

#endif
