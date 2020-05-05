#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

#include "eigen.h"

/**
 * Data structure to hold all data related to the tetrahedron method.
 */
struct TetraData{
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
};
  
/**
 * Form all tetrahedra for 3D wave vector mesh.
 * 
 * Method for creating and enumerating all the tetrahedra
 * for a given 3D mesh of wave vectors following Fig. 5 of 
 * Bloechl, Jepsen and Andersen prb 49.23 (1994): 16223.
 *
 * @param[in] grid The number of mesh points along the three lattice vectors.
 * @param[out] tetra All the data related to the analytic tetrahedron method.
 *  
 */
void formTets(const Eigen::Vector3i & grid, TetraData & tetra);

/**
 * Fill all tetrahedra with the eigenvalues.
 *
 * Method for filling the tetrahedra with the eigenvalues for
 * all polarizations. For eigenvalues are sorted along the vertex.
 *
 * @param[in] numBands Number of polarizations.
 * @param[in] energy Energy spectrum.
 * @param[out] tetra All the data related to the analytic tetrahedron method.
 *
 */
void fillTetsEigs(const long & numBands, const Eigen::MatrixXd & energy,
		TetraData & tetra);

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
double fillTetsWeights(const double & energy, const long & ib, const long & iq,
		const TetraData & tetra);
