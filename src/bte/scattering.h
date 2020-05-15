#ifndef SCATTERING_H
#define SCATTERING_H

#include "context.h"
#include "vector_bte.h"

class ScatteringMatrix {
public:
	ScatteringMatrix(Context & context_,
			FullBandStructure & activeBandStructure_,
			StatisticsSweep & statisticsSweep_);
	ScatteringMatrix(const ScatteringMatrix & that); // copy constructor
	ScatteringMatrix & operator=(const ScatteringMatrix & that);//assignment op

	VectorBTE diagonal();
	VectorBTE offDiagonalDot(VectorBTE & popOld);
	VectorBTE dot(VectorBTE & popRTA);

	// note: CGScalign only affects the results of dot() and diagonal()
	void setCGScaling();
	void unsetCGScaling();

//	std::tuple<Eigen::VectorXd,Eigen::MatrixXd> diagonalize();
 private:
	Context & context;
	FullBandStructure & innerBandStructure;
	FullBandStructure & outerBandStructure;
	StatisticsSweep & statisticsSweep;

	 // constant relaxation time approximation -> the matrix is just a scalar
	bool constantRTA = true;
	bool highMemory = true;
	bool hasCGScaling = false;

	VectorBTE internalDiagonal;
	std::vector<Eigen::MatrixXd> theMatrix;
	long numStates;
	long numPoints;
	long numCalcs;

	long deltaFunctionSelection;

	// pure virtual function
	virtual VectorBTE builderManager(VectorBTE * inPopulation) = 0;
};

#endif
