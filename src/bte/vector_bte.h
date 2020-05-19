#ifndef VECTORBTE_H
#define VECTORBTE_H

#include "eigen.h"
#include "context.h"
#include "active_bandstructure.h"

class VectorBTE {
public:
	 // constructor
	VectorBTE(Context & context_,
			FullBandStructure<FullPoints> & bandStructure_,
			const long & dimensionality_=0);
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

	void setConst(const double & constant);

//protected:
	Context & context;
	FullBandStructure<FullPoints> & bandStructure;
	long numChemPots;
	long numTemps;
	long numCalcs;
	long numStates;
	long dimensionality;
	long glob2Loc(long & imu, long & it, long & idim);
	std::tuple<long,long,long> loc2Glob(long & i);
	Eigen::MatrixXd data;
	friend class ScatteringMatrix; // this is also to remember that
	// if we change the index order of this class, we should check the
	// ScatteringMatrix implementations: they are high efficiency methods
	// so they need low-level access to the raw buffer
};

#endif
