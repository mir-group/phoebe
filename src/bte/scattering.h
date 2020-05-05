#ifndef SCATTERING_H
#define SCATTERING_H

#include "context.h"
#include "vector_bte.h"

class ScatteringMatrix {
public:
	ScatteringMatrix(Context & context_);
	VectorBTE diagonal();
	VectorBTE offDiagonalDot(VectorBTE & popOld);
	VectorBTE dot(VectorBTE & popRTA);
	void setCGScaling();

//	std::tuple<Eigen::VectorXd,Eigen::MatrixXd> diagonalize();
 private:
	Context & context;
};


#endif
