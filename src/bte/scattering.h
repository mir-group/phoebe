#ifndef SCATTERING_H
#define SCATTERING_H

#include "context.h"
#include "vector_bte.h"

class ScatteringMatrix {
public:
	ScatteringMatrix(Context & context_, StatisticsSweep & statisticsSweep_
			FullBandStructure & innerBandStructure_,
			FullBandStructure * outerBandStructure_ = nullptr);
	ScatteringMatrix(const ScatteringMatrix & that); // copy constructor
	ScatteringMatrix & operator=(const ScatteringMatrix & that);//assignment op

	VectorBTE diagonal();
	VectorBTE offDiagonalDot(VectorBTE & inPopulation);
	VectorBTE dot(VectorBTE & inPopulation);

	// note: CGScalign only affects the results of dot() and diagonal()
	void setCGScaling();
	void unsetCGScaling();

//	std::tuple<Eigen::VectorXd,Eigen::MatrixXd> diagonalize();
 private:
	Context & context;
	StatisticsSweep & statisticsSweep;
	FullBandStructure & innerBandStructure;
	FullBandStructure * outerBandStructure;

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
	// needs an implementation in every subclass
	virtual void builder(Eigen::MatrixXd * matrix, VectorBTE * linewidth,
			Vector3BTE * inPopulation, VectorBTE * outPopulation) = 0;
};

#endif
