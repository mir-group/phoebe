#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

#include "eigen.h"
#include "points.h"
#include "bandstructure.h"

/**
 * Class for everything related to the tetrahedron method.
 */
class Tetrahedra {
public:
	/**
	 * Form all tetrahedra for 3D wave vector mesh.
	 *
	 * Method for creating and enumerating all the tetrahedra
	 * for a given 3D mesh of wave vectors following Fig. 5 of
	 * Bloechl, Jepsen and Andersen prb 49.23 (1994): 16223.
	 *
	 * @param[in] grid: the mesh points along the three lattice vectors.
	 *
	 */
	Tetrahedra(FullPoints & fullPoints_, FullBandStructure & fullBandStructure_);

	/**
	 * Calculate tetrehedron weight.
	 *
	 * Method for calculating the tetrahedron weight (normalized by the number of
	 * tetrahedra) for given wave vector and polarization following Lambin and
	 * Vigneron prb 29.6 (1984): 3430.
	 *
	 * @param[in] energy Energy of mode.
	 */
	double getDOS(const double & energy);

private:
	FullPoints fullPoints;
	FullBandStructure fullBandStructure;

	/** Number of tetrahedra. */
	long numTetra;
	/** Holder for the indices of the vertices of of each tetrahedron. */
	Eigen::MatrixXi tetrahedra;
	/** Count of how many tetrahedra wave vector belongs to. */
	Eigen::VectorXi qToTetCount;
	/** Mapping of a wave vector to a tetrahedron. */
	Eigen::Tensor<long,3> qToTet;
	/** Holder for the eigenvalues. */
	Eigen::Tensor<double,3> tetraEigVals;

	/**
	 * Calculate tetrehedron weight.
	 *
	 * Method for calculating the tetrahedron weight (normalized by the number of
	 * tetrahedra) for given wave vector and polarization following Lambin and
	 * Vigneron prb 29.6 (1984): 3430.
	 *
	 * @param[in] energy Energy of mode.
	 * @param[in] ib Band index.
	 * @param[in] iq Muxed index of wave vector.
	 * @param[in] tetra All the data related to the analytic tetrahedron method.
	 * @returns The tetrahedron weight.
	 *
	 */
	double getWeight(const double & energy, const long & iq, const long & ib);
};
