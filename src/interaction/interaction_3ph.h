#ifndef PHINTERACTION_H
#define PHINTERACTION_H

#include <complex>
#include "eigen.h"
#include "state.h"
#include "crystal.h"

// * Class to calculate the coupling of a 3-phonon interaction given three
// * phonon modes.
class Interaction3Ph{
private:
	Crystal & crystal;
	long numTriplets;
	Eigen::Tensor<double,4> ifc3Tensor;
	Eigen::Tensor<double,3> cellPositions;
	Eigen::Tensor<long,2> displacedAtoms;

	std::tuple<Eigen::Tensor<std::complex<double>,3>,
			Eigen::Tensor<std::complex<double>,3>>
			calcCouplingSquared(State<FullPoints> & state1,
					State<FullPoints> & state2,
					State<FullPoints> & state3Plus,
					State<FullPoints> & state3Mins);
public:

	Interaction3Ph(Crystal & crystal_,
			long & numTriplets_,
			Eigen::Tensor<double,4> & ifc3Tensor_,
			Eigen::Tensor<double,3> & cellPositions_,
			Eigen::Tensor<long,2> & displacedAtoms_); // default constructor
	Interaction3Ph(const Interaction3Ph & that); // copy constructor
	Interaction3Ph & operator=(const Interaction3Ph & that);// assignment op

	// this is an interface: we can compute it on the fly or read the cache.
	std::tuple<Eigen::Tensor<std::complex<double>,3>,
			Eigen::Tensor<std::complex<double>,3>>
			getCouplingSquared(State<FullPoints> & state1,
					State<FullPoints> & state2,
					State<FullPoints> & state3Plus,
					State<FullPoints> & state3Mins);

	/**
	 * Calculate single three-phonon matrix element (V^{+}/V^{-1}).
	 *
	 * Method for calculating the matrix element of a single, plus
	 * or minus, three-phonon scattering process.
	 *
	 * @param[in] interactingPhonons phononTriplet data structure
	 *  containing the three phonon modes (in the full Brillouin zone)
	 *  taking part in an interaction event.
	 * @param[in] phononEigenvectors Eigenvectors of all phonons in the
	 *  full Brillouin zone.
	 * @param[in] numTriplets Number of triplets considered in the supercell
	 *  for the third order force constants calculation.
	 * @param[in] ifc3Tensor Third order force constants.
	 * @param[in] cellPositions Cartesian coordinates of the 2nd and 3rd unitcells.
	 * @param[in] displacedAtoms Index of the displaced atom for every triplet
	 *  of unitcells.
	 * @param[in] crystal Contains of crystal information
	 * @param[in] procType Character '+' or '-' to choose the +/- type process.
	 * @return The complex matrix element.
	 */

//	double calculateSingleV(const PhononTriplet &interactingPhonons, const Eigen::MatrixXd &q,
//			const int numTriplets, const Eigen::Tensor<double,4> &ifc3Tensor,
//			const Eigen::Tensor<double,3> &cellPositions,
//			const Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo,
//			const char procType);

//  void calculateAllVminus(const int *grid, const PhononMode &mode,
//			    const Eigen::MatrixXd &qFBZ,
//			    const Eigen::Tensor<complex<double>,3> &ev, const int numTriplets,
//			    const Eigen::Tensor<double,4> &ifc3Tensor,
//			    const Eigen::Tensor<double,3> &cellPositions,
//			    const Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo);
//
//  void calculateAllW(const double T,const int *grid, const PhononMode &mode,
//		     const Eigen::MatrixXi &indexMesh, const CrystalInfo &crysInfo,
//		     const Eigen::MatrixXd omega, const TetraData tetra);
};

#endif
