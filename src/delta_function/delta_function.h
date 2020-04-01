#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

using namespace std;

struct TetraData{
  //The rank-2 tensor holding the indices of the vertices of of each tetrahedron.
  Eigen::MatrixXi tetrahedra;
  //Count of how many tetrahedra wave vector belongs to.
  Eigen::VectorXi qToTetCount;
  //Mapping of a wave vector to a tetrahedron.
  Eigen::Tensor<int,3> qToTet;
};
  
/**
 * Form all tetrahedra for 3D wave vector mesh.
 * 
 * Method for creating and enumerating all the tetrahedra
 * for a given 3D mesh of wave vectors following Lambin &
 * Vigneron Physical Review B 29.6 (1984): 3430.
 *
 * @param[in] grid The number of mesh points along the three lattice vectors.
 * @param[out] tetra All the data related to the analytic tetrahedron method.
 *  
 */
void formTets(int grid[3], TetraData &tetra);
